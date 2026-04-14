"""
Conserved Gene Neighborhood (Synteny Block) Detection Module.

This module implements a spacedust-inspired algorithm for discovering conserved
gene neighborhoods across multiple genomes, built on top of PGAP2's existing
pangenome partitioning results.

The algorithm works in four steps:
  1. Build a positional index: for each strain/contig, order genes by position
     and record their cluster membership.
  2. For each genome pair, extract shared homologous hits (genes belonging to
     the same GeneCluster that appear in both genomes).
  3. Agglomerative clustering: merge nearby co-linear hits into syntenic blocks,
     scored by a combination of spatial clustering significance and gene-order
     conservation significance (following spacedust's statistical framework).
  4. Aggregate blocks across all genome pairs to identify pan-level conserved
     gene neighborhoods.

Reference:
  Zhang, R., Mirdita, M., & Söding, J. (2024). De novo discovery of conserved
  gene clusters in microbial genomes with Spacedust. bioRxiv.
"""

import math
import itertools
import networkx as nx

from tqdm import tqdm
from loguru import logger
from collections import defaultdict, Counter
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Set, Optional
from multiprocessing import get_context

from pgap2.lib.pangenome import Pangenome
from pgap2.lib.tree import Tree
from pgap2.utils.supply import tqdm_


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class GenePosition:
    """Positional record for one gene in a genome."""
    gene_id: str          # e.g. "0:1:23"
    strain_index: int
    contig_index: int
    position_index: int   # ordinal position on the contig (0-based)
    strand: str           # '+' or '-' or ''
    cluster_id: int       # root cluster id from tree.leaf_root


@dataclass
class SyntenyHit:
    """A shared homologous hit between a query and target genome."""
    cluster_id: int
    q_pos: int            # gene position index in query genome's contig
    t_pos: int            # gene position index in target genome's contig
    q_strand: str
    t_strand: str
    q_gene_id: str
    t_gene_id: str


@dataclass
class SyntenyBlock:
    """A conserved gene neighborhood / syntenic block."""
    block_id: int
    hits: List[SyntenyHit] = field(default_factory=list)
    pvalue_order: float = 1.0     # gene-order conservation P-value
    # genome pairs where this block pattern is found
    genome_pairs: List[Tuple[int, int]] = field(default_factory=list)

    @property
    def num_genes(self) -> int:
        return len(self.hits)

    @property
    def cluster_ids(self) -> Set[int]:
        return {h.cluster_id for h in self.hits}

    @property
    def q_span(self) -> int:
        positions = [h.q_pos for h in self.hits]
        return max(positions) - min(positions) + 1 if positions else 0

    @property
    def t_span(self) -> int:
        positions = [h.t_pos for h in self.hits]
        return max(positions) - min(positions) + 1 if positions else 0


@dataclass
class ConservedNeighborhood:
    """A pan-level conserved gene neighborhood aggregated across genome pairs."""
    neighborhood_id: int
    cluster_set: frozenset         # template set of cluster_ids
    blocks: List[SyntenyBlock] = field(default_factory=list)
    num_genomes: int = 0           # number of distinct genomes where this is found
    representative_block: Optional[SyntenyBlock] = None
    # Per-genome membership: strain_index -> list of GenePosition in the neighborhood
    genome_members: Dict[int, List['GenePosition']
                         ] = field(default_factory=dict)
    pvalue: float = 1.0            # combined neighborhood-level P-value (Fisher's method)

    @property
    def num_clusters(self) -> int:
        return len(self.cluster_set)

    def genome_completeness(self, strain_idx: int) -> float:
        """Fraction of the template clusters present syntenically in a genome."""
        gps = self.genome_members.get(strain_idx, [])
        present = {gp.cluster_id for gp in gps}
        return len(present) / len(self.cluster_set) if self.cluster_set else 0.0


# ---------------------------------------------------------------------------
# Statistical functions (following spacedust's framework)
# ---------------------------------------------------------------------------

def _log_gamma(x: float) -> float:
    """Log-gamma function using Lanczos approximation (Pugh 2004)."""
    if x <= 0:
        return 0.0
    return math.lgamma(x)


def _find_conserved_pairs(hits: List[SyntenyHit]) -> int:
    """
    Count the number of conserved gene pairs (consecutive hits maintaining
    strand-consistent relative order), following spacedust's logic.

    For each consecutive pair of hits (sorted by query position):
      - isSameStrand  = (h1.q_strand == h1.t_strand)
      - isSameStrand2 = (h2.q_strand == h2.t_strand)
      - isSameOrder   = (h2.t_pos > h1.t_pos)
    A pair is conserved iff both strand conditions match the order:
      (isSameStrand == isSameOrder) and (isSameStrand2 == isSameOrder)

    Falls back to the strand-unaware forward/reverse counting when strand
    information is unavailable (empty strings).
    """
    if len(hits) < 2:
        return 0

    # Sort by query position
    sorted_hits = sorted(hits, key=lambda h: h.q_pos)

    # Check if strand info is available
    has_strand = any(h.q_strand and h.t_strand for h in sorted_hits)

    if has_strand:
        # Strand-aware counting (spacedust-style)
        m = 0
        for i in range(len(sorted_hits) - 1):
            h1 = sorted_hits[i]
            h2 = sorted_hits[i + 1]
            is_same_order = h2.t_pos > h1.t_pos
            is_same_strand = (h1.q_strand == h1.t_strand)
            is_same_strand2 = (h2.q_strand == h2.t_strand)
            if (is_same_strand == is_same_order) and (is_same_strand2 == is_same_order):
                m += 1
        return m
    else:
        # Fallback: strand-unaware (original behaviour)
        forward = 0
        reverse = 0
        for i in range(len(sorted_hits) - 1):
            h1 = sorted_hits[i]
            h2 = sorted_hits[i + 1]
            if h2.t_pos > h1.t_pos:
                forward += 1
            elif h2.t_pos < h1.t_pos:
                reverse += 1
        return max(forward, reverse)


def _ordering_pvalue(k: int, m: int) -> float:
    """
    Compute log P-value for gene order conservation.

    Tests whether m conserved pairs out of k hits is more than expected
    by chance.

    log P_order = log(1 - m/k) - m*log(2) - log(Gamma(m+1))
    (when m > 0 and k > 0)
    """
    if k <= 1 or m <= 0:
        return 0.0
    ratio = m / k
    if ratio >= 1.0:
        ratio = 0.999  # cap to avoid log(0)
    log_p = math.log(max(1.0 - ratio, 1e-300)) - m * \
        math.log(2) - _log_gamma(m + 1)
    return log_p


def _chi2_survival_even_df(x: float, k: int) -> float:
    """
    Survival function P(χ²(2k) ≥ x) for chi-squared with even degrees
    of freedom (df = 2k), using the exact closed-form formula:

        P(χ²(2k) ≥ x) = e^{-x/2} · Σ_{i=0}^{k-1} (x/2)^i / i!

    Computed in log-space to avoid 0 × ∞ = NaN when x is very large.
    """
    if x <= 0:
        return 1.0
    if k <= 0:
        return 1.0

    half_x = x / 2.0

    # When the Poisson mean (half_x) greatly exceeds k, every term
    # e^{-half_x} · (half_x)^i / i!  is negligible → P ≈ 0.
    if half_x > k + 700:
        return 0.0

    # Compute each log-term: log_term_i = i·ln(half_x) − lgamma(i+1)
    # Then result = exp(−half_x) · Σ exp(log_term_i)
    #             = exp(−half_x + logsumexp(log_terms))
    log_hx = math.log(half_x) if half_x > 0 else float('-inf')
    max_lt = float('-inf')
    log_terms = []
    for i in range(k):
        lt = i * log_hx - _log_gamma(i + 1)
        log_terms.append(lt)
        if lt > max_lt:
            max_lt = lt

    # logsumexp
    acc = 0.0
    for lt in log_terms:
        acc += math.exp(lt - max_lt)
    log_sum = max_lt + math.log(acc)

    log_result = -half_x + log_sum
    if log_result < -745:
        return 0.0
    return min(math.exp(log_result), 1.0)


def _fisher_combined_pvalue(pvalues: List[float]) -> float:
    """
    Fisher's combined probability test (meta-analysis of P-values).

    Given k independent P-values p₁, p₂, ..., pₖ, Fisher's statistic is:
        χ² = -2 Σ ln(pᵢ)
    Under H₀, χ² ~ χ²(2k).

    Returns the combined P-value: P(χ²(2k) ≥ observed).

    This is used to aggregate pairwise syntenic block P-values into a
    single neighborhood-level P-value, answering: "Is this conserved gene
    neighborhood collectively more significant than expected by chance
    across all genome pairs where it was observed?"
    """
    if not pvalues:
        return 1.0

    k = len(pvalues)
    chi2_stat = 0.0
    for p in pvalues:
        # Clamp to avoid log(0)
        p_clamped = max(p, 1e-300)
        chi2_stat += -2.0 * math.log(p_clamped)

    return _chi2_survival_even_df(chi2_stat, k)


def _order_pvalue_for_block(hits: List[SyntenyHit]) -> float:
    """
    Compute the gene-order conservation P-value for a set of hits.

    Returns: pvalue_order (float in [0, 1])
    """
    k = len(hits)
    if k < 2:
        return 1.0

    conserved_pairs = _find_conserved_pairs(hits)
    log_p_order = _ordering_pvalue(k, conserved_pairs)

    # Clamp: log(p) must be ≤ 0 (p ≤ 1)
    log_p_order = min(log_p_order, 0.0)
    p_order = math.exp(max(log_p_order, -500))

    return p_order


# ---------------------------------------------------------------------------
# Step 1: Build positional index
# ---------------------------------------------------------------------------

def build_position_index(
    pg: Pangenome,
    tree: Tree,
    G: nx.Graph
) -> Dict[int, Dict[int, List[GenePosition]]]:
    """
    Build a positional index: strain -> contig -> sorted list of GenePosition.

    For each gene in the pangenome, record its cluster membership (the root
    of the connected component in tree.leaf_root via the merged network G).
    """
    # Build a mapping from gene_id -> cluster_id (using G's node membership)
    gene_to_cluster = {}
    for node_id, node_data in G.nodes(data=True):
        # Each node in G has 'members' = set of gene_ids
        # and 'repre_nodes' = list of representative node names
        members = node_data.get('members', set())
        for gene_id in members:
            gene_to_cluster[gene_id] = node_id

    # Build position index per strain per contig
    position_index: Dict[int, Dict[int, List[GenePosition]]] = defaultdict(
        lambda: defaultdict(list)
    )

    for strain_index, strain in pg.strain_dict.items():
        gene_num = strain.gene_num  # list of gene counts per contig
        for contig_index, gene_count in enumerate(gene_num):
            for gene_index in range(gene_count):
                gene_id = f"{strain_index}:{contig_index}:{gene_index}"
                cluster_id = gene_to_cluster.get(gene_id, None)
                if cluster_id is None:
                    continue
                annot_entry = pg.annot.get(gene_id, {}) if pg.annot else {}
                strand = annot_entry.get('strand', '')
                gp = GenePosition(
                    gene_id=gene_id,
                    strain_index=strain_index,
                    contig_index=contig_index,
                    position_index=gene_index,
                    strand=strand,
                    cluster_id=cluster_id
                )
                position_index[strain_index][contig_index].append(gp)

    # Sort each contig's genes by position_index (should already be sorted)
    for strain_idx in position_index:
        for contig_idx in position_index[strain_idx]:
            position_index[strain_idx][contig_idx].sort(
                key=lambda gp: gp.position_index
            )

    return {s: dict(contigs) for s, contigs in position_index.items()}


# ---------------------------------------------------------------------------
# Step 2: Extract shared hits per genome pair
# ---------------------------------------------------------------------------

def _build_cluster_to_genes(
    position_index: Dict[int, Dict[int, List[GenePosition]]]
) -> Dict[int, Dict]:
    """
    Build a mapping: cluster_id -> {strain_index: [(contig, pos_index, strand, gene_id), ...]}
    """
    cluster_genes: Dict[int, Dict[int, list]
                        ] = defaultdict(lambda: defaultdict(list))
    for strain_idx, contigs in position_index.items():
        for contig_idx, genes in contigs.items():
            for gp in genes:
                cluster_genes[gp.cluster_id][strain_idx].append(
                    (gp.contig_index, gp.position_index, gp.strand, gp.gene_id)
                )
    return dict(cluster_genes)


def _build_strain_to_clusters(
    cluster_genes: Dict[int, Dict[int, list]]
) -> Dict[int, Set[int]]:
    """
    Build an inverted index: strain_index -> set of cluster_ids present in that strain.
    Used for fast set-intersection to find shared clusters between a genome pair.
    """
    strain_clusters: Dict[int, Set[int]] = defaultdict(set)
    for cluster_id, strain_map in cluster_genes.items():
        for strain_idx in strain_map:
            strain_clusters[strain_idx].add(cluster_id)
    return dict(strain_clusters)


def extract_hits_for_pair(
    strain_q: int,
    strain_t: int,
    cluster_genes: Dict,
    shared_clusters: Set[int],
) -> Dict[Tuple[int, int], List[SyntenyHit]]:
    """
    For a pair of genomes, extract all shared hits grouped by (contig_q, contig_t).
    Only iterates over pre-computed shared clusters (set intersection),
    not the full cluster set.

    Returns: dict of (contig_q, contig_t) -> list of SyntenyHit
    """
    hits_by_contig_pair: Dict[Tuple[int, int],
                              List[SyntenyHit]] = defaultdict(list)

    for cluster_id in shared_clusters:
        strain_map = cluster_genes[cluster_id]
        q_genes = strain_map[strain_q]
        t_genes = strain_map[strain_t]
        # For each q-t gene combination in this cluster
        for q_contig, q_pos, q_strand, q_gid in q_genes:
            for t_contig, t_pos, t_strand, t_gid in t_genes:
                hit = SyntenyHit(
                    cluster_id=cluster_id,
                    q_pos=q_pos, t_pos=t_pos,
                    q_strand=q_strand, t_strand=t_strand,
                    q_gene_id=q_gid, t_gene_id=t_gid
                )
                hits_by_contig_pair[(q_contig, t_contig)].append(hit)

    return dict(hits_by_contig_pair)


# ---------------------------------------------------------------------------
# Step 3: Agglomerative clustering to discover syntenic blocks
# ---------------------------------------------------------------------------

def find_syntenic_blocks(
    hits: List[SyntenyHit],
    max_gap: int = 3,
    min_genes: int = 2,
    pval_threshold: float = 0.01,
    use_order: bool = True,
) -> List[SyntenyBlock]:
    """
    Union-Find sweep to discover syntenic blocks — O(n · max_gap) instead of O(n²).

    Algorithm:
      1. Sort hits by query position.
      2. For each hit, check a small forward window of subsequent hits (bounded
         by max_gap in query space).  If a pair is also within max_gap in
         target space, union the two hits.
      3. Collect connected components, optionally filter by P_order, and return.

    When use_order=True, blocks are filtered by gene-order P-value ≤ pval_threshold.
    When use_order=False, blocks are accepted purely based on min_genes.
    """
    if len(hits) < min_genes:
        return []

    # --- Union-Find with path compression ---
    parent = list(range(len(hits)))
    rank = [0] * len(hits)

    def _find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def _union(x: int, y: int) -> None:
        rx, ry = _find(x), _find(y)
        if rx == ry:
            return
        if rank[rx] < rank[ry]:
            rx, ry = ry, rx
        parent[ry] = rx
        if rank[rx] == rank[ry]:
            rank[rx] += 1

    # Sort by query position
    hits_sorted = sorted(hits, key=lambda h: h.q_pos)
    n = len(hits_sorted)

    # Sweep: for each hit, look ahead at subsequent hits within max_gap
    # in query space. Because hits are sorted by q_pos, we can break early.
    for i in range(n):
        qi = hits_sorted[i].q_pos
        for j in range(i + 1, n):
            qj = hits_sorted[j].q_pos
            if qj - qi > max_gap + 1:  # gap = qj - qi - 1 > max_gap
                break
            # Check target-side gap
            t_gap = abs(hits_sorted[j].t_pos - hits_sorted[i].t_pos) - 1
            if t_gap <= max_gap:
                _union(i, j)

    # Collect groups
    groups: Dict[int, List[SyntenyHit]] = defaultdict(list)
    for i in range(n):
        groups[_find(i)].append(hits_sorted[i])

    # Filter
    results: List[SyntenyBlock] = []
    block_counter = 0
    for group_hits in groups.values():
        if len(group_hits) < min_genes:
            continue
        p_order = _order_pvalue_for_block(group_hits) if use_order else 1.0
        if use_order and p_order > pval_threshold:
            continue
        block = SyntenyBlock(
            block_id=block_counter,
            hits=group_hits,
            pvalue_order=p_order,
        )
        results.append(block)
        block_counter += 1

    return results


# ---------------------------------------------------------------------------
# Step 4: Aggregate across genome pairs -> pan-level neighborhoods
# ---------------------------------------------------------------------------

def aggregate_neighborhoods(
    all_blocks: List[Tuple[int, int, SyntenyBlock]],
    position_index: Dict[int, Dict[int, List[GenePosition]]],
    min_genomes: int = 2,
    min_genes: int = 2,
    max_gap: int = 3,
    resolution: float = 1.0,
    use_order: bool = True,
) -> List[ConservedNeighborhood]:
    """
    Graph-based aggregation of syntenic blocks into conserved neighborhoods.

    Instead of requiring exact cluster-composition identity (frozenset match),
    this algorithm uses community detection on a cluster co-occurrence graph
    to naturally tolerate partial overlap and support transitive merging.

    Algorithm
    ---------
    A. **Build cluster co-occurrence graph** from all raw syntenic blocks.
       Nodes = cluster IDs; an edge (A, B) is weighted by the number of
       distinct genomes where A and B co-occur in any syntenic block.
    B. **Prune** edges observed in fewer than *min_genomes* genomes and
       remove resulting isolate nodes.
    C. **Community detection** (Louvain, weight-aware) partitions the graph
       into communities — each community defines a *neighborhood template*.
    D. **Per-genome validation**: for each template, scan every candidate
       genome's contigs (via *position_index*) to find the best syntenic
       instance (allowing partial membership, i.e. missing a few clusters).
    E. **Filter**: keep only neighborhoods present in ≥ *min_genomes* genomes.

    This approach:
      • tolerates missing genes (genome with 4/5 template clusters still
        counts, with completeness = 0.80);
      • supports transitive merging (blocks A–B and B–C with overlapping
        cluster sets are naturally merged into one community);
      • leverages graph-theoretic community structure instead of exact
        set identity.

    Args:
        all_blocks:     list of (strain_q, strain_t, SyntenyBlock)
        position_index: strain -> contig -> sorted list of GenePosition
        min_genomes:    minimum genomes a neighborhood must span
        min_genes:      minimum clusters in a neighborhood template
        max_gap:        maximum positional gap for syntenic grouping

    Returns:
        List of ConservedNeighborhood objects, sorted by num_genomes desc.
    """
    import time as _time

    if not all_blocks:
        return []

    t0 = _time.perf_counter()

    # ---- A: Build cluster co-occurrence graph from raw blocks ----------
    # Track per-edge: genome set (for min_genomes filter) and
    # edge weight: block count (each block contributes 1.0).
    edge_genomes: Dict[Tuple[int, int], Set[int]] = {}
    edge_weight: Dict[Tuple[int, int], float] = {}
    for sq, st, block in all_blocks:
        cids = sorted(block.cluster_ids)
        nc = len(cids)
        for i in range(nc):
            a = cids[i]
            for j in range(i + 1, nc):
                key = (a, cids[j])
                gs = edge_genomes.get(key)
                if gs is None:
                    edge_genomes[key] = {sq, st}
                    edge_weight[key] = 1.0
                else:
                    gs.add(sq)
                    gs.add(st)
                    edge_weight[key] += 1.0

    t1 = _time.perf_counter()
    logger.debug(f"  Step A (co-occurrence dict): {t1 - t0:.3f}s, "
                 f"{len(edge_genomes)} edges from {len(all_blocks)} blocks")

    # ---- B: Prune weak edges & build nx.Graph in one pass ---------------
    N = nx.Graph()
    for (a, b), gs in edge_genomes.items():
        if len(gs) >= min_genomes:
            N.add_edge(a, b, weight=edge_weight[(a, b)])
    del edge_genomes, edge_weight  # free memory

    N.remove_nodes_from(list(nx.isolates(N)))
    if len(N) == 0:
        return []

    t2 = _time.perf_counter()
    logger.debug(f"  Step B (prune + build graph): {t2 - t1:.3f}s, "
                 f"{N.number_of_nodes()} nodes, {N.number_of_edges()} edges")

    # ---- C: Community detection (Louvain) ------------------------------
    try:
        communities = nx.community.louvain_communities(
            N, weight='weight', resolution=resolution, seed=42)
    except AttributeError:
        communities = list(nx.connected_components(N))

    communities = [c for c in communities if len(c) >= min_genes]

    t3 = _time.perf_counter()
    logger.info(
        f"Co-occurrence graph: {N.number_of_nodes()} nodes, "
        f"{N.number_of_edges()} edges → {len(communities)} communities "
        f"({t3 - t0:.2f}s)")
    del N  # no longer needed

    if not communities:
        return []

    # ---- Pre-index: (strain, contig, cluster) → sorted positions -------
    # This replaces the naive "filter all genes" approach in Step D.
    # Index: strain → cluster → [(contig, pos, GenePosition), ...]
    strain_cluster_positions: Dict[int, Dict[int, List[Tuple[int, int, GenePosition]]]] = \
        defaultdict(lambda: defaultdict(list))
    cluster_to_strains: Dict[int, Set[int]] = defaultdict(set)
    for strain_idx, contigs in position_index.items():
        for contig_idx, genes in contigs.items():
            for gp in genes:
                strain_cluster_positions[strain_idx][gp.cluster_id].append(
                    (contig_idx, gp.position_index, gp))
                cluster_to_strains[gp.cluster_id].add(strain_idx)

    t4 = _time.perf_counter()
    logger.debug(f"  Pre-index: {t4 - t3:.3f}s")

    # ---- Pre-build: community → best block (skip full copy) ------------
    cluster_to_community: Dict[int, int] = {}
    for comm_idx, comm_set in enumerate(communities):
        for cid in comm_set:
            cluster_to_community[cid] = comm_idx

    # For each community, track the best block (lowest P_order)
    # and collect per-block P-values for Fisher's combined test.
    comm_best_block: Dict[int, SyntenyBlock] = {}
    comm_best_block_p: Dict[int, float] = {}  # p_order of best block
    comm_block_count: Dict[int, int] = defaultdict(int)
    comm_block_pvals: Dict[int, List[float]] = defaultdict(list)
    for sq, st, block in all_blocks:
        comm_votes: Dict[int, int] = Counter()
        for cid in block.cluster_ids:
            ci = cluster_to_community.get(cid)
            if ci is not None:
                comm_votes[ci] += 1
        block_p = block.pvalue_order
        for ci, count in comm_votes.items():
            if count >= min_genes:
                comm_block_count[ci] += 1
                if use_order:
                    comm_block_pvals[ci].append(block_p)
                prev_p = comm_best_block_p.get(ci)
                if prev_p is None or block_p < prev_p:
                    comm_best_block[ci] = block
                    comm_best_block_p[ci] = block_p

    t5 = _time.perf_counter()
    logger.debug(f"  Pre-build community index: {t5 - t4:.3f}s")

    # ---- D: Per-genome membership validation ---------------------------
    neighborhoods: List[ConservedNeighborhood] = []
    neigh_id = 0

    for comm_idx, template_set in enumerate(communities):
        template = frozenset(template_set)

        # Fast pre-filter: genomes with ≥ min_genes template clusters
        strain_counts: Dict[int, int] = Counter()
        for cid in template:
            for sid in cluster_to_strains.get(cid, set()):
                strain_counts[sid] += 1
        candidate_strains = {
            s for s, c in strain_counts.items() if c >= min_genes}

        # For each candidate genome, find its best syntenic instance.
        # Use the pre-built (strain, cluster) → positions index for
        # O(|template| × hits_per_cluster) instead of O(|all genes|).
        genome_members: Dict[int, List[GenePosition]] = {}

        for strain_idx in candidate_strains:
            scp = strain_cluster_positions.get(strain_idx, {})

            # Gather all relevant (contig, pos, gp) tuples in one pass
            by_contig: Dict[int, List[GenePosition]] = defaultdict(list)
            for cid in template:
                for contig_idx, pos_idx, gp in scp.get(cid, []):
                    by_contig[contig_idx].append(gp)

            best_instance: List[GenePosition] = []
            best_n_clusters = 0

            for contig_idx, relevant in by_contig.items():
                if len(relevant) < min_genes:
                    continue

                # Gap-based sweep to group nearby template genes
                relevant.sort(key=lambda gp: gp.position_index)
                current_group: List[GenePosition] = [relevant[0]]

                for k in range(1, len(relevant)):
                    gap = (relevant[k].position_index
                           - relevant[k - 1].position_index - 1)
                    if gap <= max_gap:
                        current_group.append(relevant[k])
                    else:
                        n_clust = len({g.cluster_id
                                       for g in current_group})
                        if (len(current_group) >= min_genes
                                and n_clust > best_n_clusters):
                            best_instance = list(current_group)
                            best_n_clusters = n_clust
                        current_group = [relevant[k]]

                # Flush last group
                n_clust = len({g.cluster_id for g in current_group})
                if (len(current_group) >= min_genes
                        and n_clust > best_n_clusters):
                    best_instance = list(current_group)
                    best_n_clusters = n_clust

            if best_n_clusters >= min_genes:
                genome_members[strain_idx] = best_instance

        # ---- E: Filter by min_genomes ----------------------------------
        if len(genome_members) < min_genomes:
            continue

        # Use tracked best block — no mass copying
        best_block = comm_best_block.get(comm_idx)

        # Compute neighborhood-level P-value via Fisher's combined test
        # on all pairwise block P_order values assigned to this community.
        # When use_order=False, no P-values are available — set to 1.0.
        if use_order:
            block_pvals = comm_block_pvals.get(comm_idx, [])
            combined_pval = _fisher_combined_pvalue(block_pvals)
        else:
            combined_pval = 1.0

        neigh = ConservedNeighborhood(
            neighborhood_id=neigh_id,
            cluster_set=template,
            blocks=[],     # raw block list omitted for performance
            num_genomes=len(genome_members),
            representative_block=best_block,
            genome_members=genome_members,
            pvalue=combined_pval,
        )
        neighborhoods.append(neigh)
        neigh_id += 1

    t6 = _time.perf_counter()
    logger.debug(f"  Step D+E (validation + assembly): {t6 - t5:.3f}s, "
                 f"{len(neighborhoods)} neighborhoods from {len(communities)} communities")

    # Sort by num_genomes desc, then by combined P-value asc
    neighborhoods.sort(
        key=lambda n: (
            -n.num_genomes,
            n.pvalue))

    logger.debug(f"  aggregate_neighborhoods total: {t6 - t0:.3f}s")
    return neighborhoods


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Worker for multiprocessing
# ---------------------------------------------------------------------------

# Module-level globals set by the pool initializer
_MP_CLUSTER_GENES = None
_MP_STRAIN_CLUSTERS = None
_MP_MAX_GAP = 3
_MP_MIN_GENES = 2
_MP_PVAL = 0.01
_MP_USE_ORDER = True


def _init_worker(cluster_genes, strain_clusters,
                 max_gap, min_genes, pval_threshold, use_order):
    """Initializer for multiprocessing workers — avoids pickling large dicts."""
    global _MP_CLUSTER_GENES, _MP_STRAIN_CLUSTERS
    global _MP_MAX_GAP, _MP_MIN_GENES, _MP_PVAL, _MP_USE_ORDER
    _MP_CLUSTER_GENES = cluster_genes
    _MP_STRAIN_CLUSTERS = strain_clusters
    _MP_MAX_GAP = max_gap
    _MP_MIN_GENES = min_genes
    _MP_PVAL = pval_threshold
    _MP_USE_ORDER = use_order


def _process_pair(pair: Tuple[int, int]) -> List[Tuple[int, int, SyntenyBlock]]:
    """Process a single genome pair — designed for multiprocessing.Pool.map."""
    sq, st = pair
    shared = _MP_STRAIN_CLUSTERS[sq] & _MP_STRAIN_CLUSTERS[st]
    if len(shared) < _MP_MIN_GENES:
        return []

    hits_by_contig = extract_hits_for_pair(
        sq, st, _MP_CLUSTER_GENES, shared)

    result: List[Tuple[int, int, SyntenyBlock]] = []
    for hits in hits_by_contig.values():
        blocks = find_syntenic_blocks(
            hits,
            max_gap=_MP_MAX_GAP,
            min_genes=_MP_MIN_GENES,
            pval_threshold=_MP_PVAL,
            use_order=_MP_USE_ORDER,
        )
        for block in blocks:
            result.append((sq, st, block))
    return result


def discover_neighborhoods(
    G: nx.Graph,
    pg: Pangenome,
    tree: Tree,
    max_gap: int = 3,
    min_genes: int = 2,
    pval_threshold: float = 0.01,
    min_genomes: int = 2,
    resolution: float = 1.0,
    use_order: bool = True,
    disable: bool = False,
    threads: int = 1,
) -> List[ConservedNeighborhood]:
    """
    Full pipeline: discover conserved gene neighborhoods from a partitioned pangenome.

    Args:
        G: The pangenome network graph (output of partition).
        pg: The Pangenome object.
        tree: The identity tree.
        max_gap: Maximum gap (in gene positions) allowed between clusters to merge.
        min_genes: Minimum number of genes to define a syntenic block.
        pval_threshold: P-value threshold for reporting a block.
        min_genomes: Minimum number of genomes a neighborhood must span.
        disable: Disable progress bar.
        threads: Number of parallel worker processes.

    Returns:
        List of ConservedNeighborhood objects.
    """
    logger.info("Building gene position index...")
    position_index = build_position_index(pg, tree, G)
    logger.info(f"Position index built for {len(position_index)} strains")

    logger.info("Building cluster-to-gene mapping...")
    cluster_genes = _build_cluster_to_genes(position_index)
    logger.info(
        f"Found {len(cluster_genes)} gene clusters with positional data")

    logger.info("Building strain-to-cluster index...")
    strain_clusters = _build_strain_to_clusters(cluster_genes)

    strain_indices = sorted(position_index.keys())
    pairs = list(itertools.combinations(strain_indices, 2))
    num_pairs = len(pairs)
    logger.info(f"Scanning {num_pairs} genome pairs for syntenic blocks...")

    all_blocks: List[Tuple[int, int, SyntenyBlock]] = []

    if threads > 1:
        # --- Parallel mode ---
        ctx = get_context('fork')
        with ctx.Pool(
            processes=threads,
            initializer=_init_worker,
            initargs=(cluster_genes, strain_clusters,
                      max_gap, min_genes, pval_threshold, use_order),
        ) as pool:
            for pair_blocks in tqdm(
                pool.imap_unordered(_process_pair, pairs, chunksize=64),
                total=num_pairs,
                unit=" pair",
                desc=tqdm_.step(-1),
                disable=disable,
            ):
                all_blocks.extend(pair_blocks)
    else:
        # --- Serial mode (avoids fork overhead for small jobs) ---
        _init_worker(cluster_genes, strain_clusters,
                     max_gap, min_genes, pval_threshold, use_order)
        for pair in tqdm(
            pairs,
            total=num_pairs,
            unit=" pair",
            desc=tqdm_.step(-1),
            disable=disable,
        ):
            all_blocks.extend(_process_pair(pair))

    logger.info(
        f"Found {len(all_blocks)} raw syntenic blocks across all genome pairs")

    logger.info("Aggregating into pan-level conserved neighborhoods...")
    neighborhoods = aggregate_neighborhoods(
        all_blocks,
        position_index=position_index,
        min_genomes=min_genomes,
        min_genes=min_genes,
        max_gap=max_gap,
        resolution=resolution,
        use_order=use_order,
    )
    logger.info(
        f"Identified {len(neighborhoods)} conserved gene neighborhoods")

    return neighborhoods


# ---------------------------------------------------------------------------
# Output functions
# ---------------------------------------------------------------------------

def write_neighborhood_tsv(
    neighborhoods: List[ConservedNeighborhood],
    pg: Pangenome,
    G: nx.Graph,
    outdir: str,
    prefix: str = "pgap2.neighborhood"
):
    """
    Write the neighborhood results to a TSV file.

    Output columns:
      neighborhood_id, num_clusters, num_genomes, pvalue, pval_order
    """
    outfile = f"{outdir}/{prefix}.synteny_blocks.tsv"
    with open(outfile, 'w') as fh:
        # Header
        fh.write("#neighborhood_id\tnum_clusters\tnum_genomes\t"
                 "pvalue\tpval_order\n")

        for neigh in neighborhoods:
            rep = neigh.representative_block
            if rep is None:
                continue

            fh.write(f"neigh_{neigh.neighborhood_id}\t"
                     f"{neigh.num_clusters}\t"
                     f"{neigh.num_genomes}\t"
                     f"{neigh.pvalue:.2e}\t"
                     f"{rep.pvalue_order:.2e}\n")

    logger.info(f"Neighborhood results written to {outfile}")
    return outfile


def write_neighborhood_detail_tsv(
    neighborhoods: List[ConservedNeighborhood],
    pg: Pangenome,
    G: nx.Graph,
    outdir: str,
    prefix: str = "pgap2.neighborhood"
):
    """
    Write per-genome membership details for each neighborhood.

    Each row represents one genome's instance of a neighborhood, showing
    which clusters are present, the completeness fraction, and gene IDs.
    """
    outfile = f"{outdir}/{prefix}.synteny_detail.tsv"
    with open(outfile, 'w') as fh:
        fh.write("#neighborhood_id\tgenome\tnum_present\tnum_total\t"
                 "completeness\tgene_ids\tcluster_ids\n")

        for neigh in neighborhoods:
            for strain_idx in sorted(neigh.genome_members.keys()):
                gp_list = neigh.genome_members[strain_idx]
                strain_name = (
                    pg.strain_dict[strain_idx].strain_name
                    if strain_idx in pg.strain_dict
                    else str(strain_idx))
                present_clusters = {gp.cluster_id for gp in gp_list}
                completeness = neigh.genome_completeness(strain_idx)

                gene_ids = ";".join(
                    gp.gene_id for gp in sorted(
                        gp_list, key=lambda g: g.position_index))
                cluster_ids = ",".join(
                    str(c) for c in sorted(
                        present_clusters, key=lambda x: str(x)))

                fh.write(f"neigh_{neigh.neighborhood_id}\t"
                         f"{strain_name}\t"
                         f"{len(present_clusters)}\t"
                         f"{neigh.num_clusters}\t"
                         f"{completeness:.2f}\t"
                         f"{gene_ids}\t"
                         f"{cluster_ids}\n")

    logger.info(f"Neighborhood membership detail written to {outfile}")
    return outfile


def write_neighborhood_gml(
    neighborhoods: List[ConservedNeighborhood],
    outdir: str,
    prefix: str = "pgap2.neighborhood"
):
    """
    Write a gene neighborhood network in GML format.

    Nodes = gene clusters (labelled with the neighborhood they belong to),
    edges = co-occurrence in a neighborhood (weighted by number of genomes).
    """
    H = nx.Graph()

    for neigh in neighborhoods:
        cluster_list = sorted(neigh.cluster_set, key=lambda x: str(x))
        for cid in cluster_list:
            nid = str(cid)
            if not H.has_node(nid):
                H.add_node(nid, neighborhood=neigh.neighborhood_id)
            else:
                # Node may belong to multiple neighborhoods — record all
                prev = H.nodes[nid].get('neighborhood', '')
                H.nodes[nid]['neighborhood'] = (
                    f"{prev},{neigh.neighborhood_id}"
                    if prev != '' else neigh.neighborhood_id)
        for i in range(len(cluster_list)):
            for j in range(i + 1, len(cluster_list)):
                a, b = str(cluster_list[i]), str(cluster_list[j])
                if H.has_edge(a, b):
                    H[a][b]['weight'] += neigh.num_genomes
                else:
                    H.add_edge(a, b, weight=neigh.num_genomes)

    outfile = f"{outdir}/{prefix}.synteny_network.gml"
    nx.write_gml(H, outfile)
    logger.info(f"Neighborhood network written to {outfile}")
    return outfile
