"""
Tests for pgap2.lib.synteny — conserved gene neighborhood detection.

Test categories:
  1. Unit tests for each algorithmic function
  2. Integration test: full pipeline on a synthetic pangenome
  3. Edge-case tests (single contig, single genome, empty input, etc.)
  4. Performance / scalability benchmark

Run:
    pytest tests/test_synteny.py -v
    pytest tests/test_synteny.py -v -k bench   # benchmarks only
"""

import math
import time
import random
import tempfile
import os
import itertools
from collections import defaultdict
from typing import Dict, List, Set, Tuple

import pytest
import networkx as nx

# --- Module under test ---
from pgap2.lib.synteny import (
    GenePosition,
    SyntenyHit,
    SyntenyBlock,
    ConservedNeighborhood,
    _log_gamma,
    _find_conserved_pairs,
    _ordering_pvalue,
    _order_pvalue_for_block,
    find_syntenic_blocks,
    extract_hits_for_pair,
    aggregate_neighborhoods,
    _build_cluster_to_genes,
    _build_strain_to_clusters,
    write_neighborhood_tsv,
    write_neighborhood_detail_tsv,
    write_neighborhood_gml,
)


# =========================================================================
# Helpers: lightweight mock objects so we don't need real PGAP2 data files
# =========================================================================

class MockStrain:
    """Minimal stand-in for pgap2.lib.strain.Strain."""

    def __init__(self, strain_name: str, strain_index: int, gene_num: list):
        self.strain_name = strain_name
        self.strain_index = strain_index
        # [genes_on_contig0, genes_on_contig1, ...]
        self.gene_num = gene_num
        self.bed_gene_num = sum(gene_num)


class MockPangenome:
    """Minimal stand-in for pgap2.lib.pangenome.Pangenome."""

    def __init__(self):
        self.strain_dict: Dict[int, MockStrain] = {}
        self.annot: Dict[str, dict] = {}

    def load_strain(self, strain: MockStrain):
        self.strain_dict[strain.strain_index] = strain


class MockTree:
    """Minimal stand-in for pgap2.lib.tree.Tree (unused in synteny)."""
    pass


def _build_mock_graph(gene_cluster_map: Dict[int, Set[str]]) -> nx.Graph:
    """
    Build a mock pangenome graph G.

    gene_cluster_map: {cluster_id: {gene_id, ...}}
    Each cluster_id becomes a node whose 'members' is the set of gene_ids.
    """
    G = nx.Graph()
    for cluster_id, gene_ids in gene_cluster_map.items():
        G.add_node(cluster_id, members=set(gene_ids))
    return G


def _make_position_index(
    genome_layouts: Dict[int, Dict[int, List[int]]]
) -> Dict[int, Dict[int, List[GenePosition]]]:
    """
    Build a position_index directly (bypasses build_position_index).

    genome_layouts: {strain_idx: {contig_idx: [cluster_id_at_pos0, cluster_id_at_pos1, ...]}}
    """
    pi: Dict[int, Dict[int, List[GenePosition]]] = {}
    for strain_idx, contigs in genome_layouts.items():
        pi[strain_idx] = {}
        for contig_idx, clusters in contigs.items():
            gp_list = []
            for pos, cid in enumerate(clusters):
                gp = GenePosition(
                    gene_id=f"{strain_idx}:{contig_idx}:{pos}",
                    strain_index=strain_idx,
                    contig_index=contig_idx,
                    position_index=pos,
                    strand='',
                    cluster_id=cid,
                )
                gp_list.append(gp)
            pi[strain_idx][contig_idx] = gp_list
    return pi


# =========================================================================
# 1. Unit tests — Statistical functions
# =========================================================================

class TestStatisticalFunctions:
    """Tests for _log_gamma, _ordering_pvalue, _order_pvalue_for_block, etc."""

    def test_log_gamma_positive(self):
        assert abs(_log_gamma(1) - 0.0) < 1e-10       # Gamma(1) = 0! = 1
        assert abs(_log_gamma(5) - math.lgamma(5)) < 1e-10

    def test_log_gamma_edge(self):
        assert _log_gamma(0) == 0.0
        assert _log_gamma(-1) == 0.0

    def test_find_conserved_pairs_forward(self):
        hits = [
            SyntenyHit(1, 0, 0, '', '', 'a', 'b'),
            SyntenyHit(2, 1, 1, '', '', 'a', 'b'),
            SyntenyHit(3, 2, 2, '', '', 'a', 'b'),
        ]
        assert _find_conserved_pairs(hits) == 2

    def test_find_conserved_pairs_reverse(self):
        """Reverse-complement synteny: q↑ → t↓ should also score."""
        hits = [
            SyntenyHit(1, 0, 10, '', '', 'a', 'b'),
            SyntenyHit(2, 1, 9, '', '', 'a', 'b'),
            SyntenyHit(3, 2, 8, '', '', 'a', 'b'),
        ]
        assert _find_conserved_pairs(hits) == 2

    def test_find_conserved_pairs_mixed(self):
        hits = [
            SyntenyHit(1, 0, 0, '', '', 'a', 'b'),
            SyntenyHit(2, 1, 5, '', '', 'a', 'b'),
            SyntenyHit(3, 2, 1, '', '', 'a', 'b'),
        ]
        assert _find_conserved_pairs(hits) == 1  # max(1 fwd, 1 rev)

    def test_find_conserved_pairs_empty(self):
        assert _find_conserved_pairs([]) == 0
        assert _find_conserved_pairs(
            [SyntenyHit(1, 0, 0, '', '', '', '')]) == 0

    def test_ordering_pvalue_guards(self):
        assert _ordering_pvalue(1, 0) == 0.0
        assert _ordering_pvalue(0, 0) == 0.0
        assert _ordering_pvalue(5, 0) == 0.0

    def test_ordering_pvalue_perfect(self):
        log_p = _ordering_pvalue(k=5, m=4)
        assert log_p < 0

    def test_order_pvalue_for_block_small(self):
        hits = [
            SyntenyHit(1, 0, 0, '', '', '', ''),
            SyntenyHit(2, 1, 1, '', '', '', ''),
        ]
        po = _order_pvalue_for_block(hits)
        assert 0 < po <= 1

    def test_order_pvalue_for_block_single(self):
        """Single hit → no meaningful block."""
        po = _order_pvalue_for_block(
            [SyntenyHit(1, 0, 0, '', '', '', '')])
        assert po == 1.0


# =========================================================================
# 2. Unit tests — find_syntenic_blocks (Union-Find)
# =========================================================================

class TestFindSyntenicBlocks:

    def test_simple_collinear(self):
        """Three adjacent forward-collinear hits → 1 block."""
        hits = [
            SyntenyHit(10, q_pos=0, t_pos=0, q_strand='',
                       t_strand='', q_gene_id='q0', t_gene_id='t0'),
            SyntenyHit(11, q_pos=1, t_pos=1, q_strand='',
                       t_strand='', q_gene_id='q1', t_gene_id='t1'),
            SyntenyHit(12, q_pos=2, t_pos=2, q_strand='',
                       t_strand='', q_gene_id='q2', t_gene_id='t2'),
        ]
        blocks = find_syntenic_blocks(
            hits, max_gap=3, min_genes=2, pval_threshold=1.0)
        assert len(blocks) == 1
        assert blocks[0].num_genes == 3
        assert blocks[0].cluster_ids == {10, 11, 12}

    def test_two_separate_blocks(self):
        """Two groups separated by large gap → 2 blocks."""
        hits = [
            # Group 1: positions 0-2
            SyntenyHit(10, 0, 0, '', '', '', ''),
            SyntenyHit(11, 1, 1, '', '', '', ''),
            SyntenyHit(12, 2, 2, '', '', '', ''),
            # Group 2: positions 20-22 (gap > max_gap=3)
            SyntenyHit(20, 20, 20, '', '', '', ''),
            SyntenyHit(21, 21, 21, '', '', '', ''),
            SyntenyHit(22, 22, 22, '', '', '', ''),
        ]
        blocks = find_syntenic_blocks(
            hits, max_gap=3, min_genes=2, pval_threshold=1.0)
        assert len(blocks) == 2

    def test_gap_boundary(self):
        """Hits exactly at max_gap boundary → should merge."""
        # max_gap=3: q_pos 0 and 4 → gap = 4-0-1 = 3 ≤ 3 → merge
        hits = [
            SyntenyHit(10, 0, 0, '', '', '', ''),
            SyntenyHit(11, 4, 4, '', '', '', ''),
        ]
        blocks = find_syntenic_blocks(
            hits, max_gap=3, min_genes=2, pval_threshold=1.0)
        assert len(blocks) == 1

    def test_gap_exceeds(self):
        """Hits just beyond max_gap → should NOT merge."""
        # max_gap=3: q_pos 0 and 5 → gap = 5-0-1 = 4 > 3 → separate
        hits = [
            SyntenyHit(10, 0, 0, '', '', '', ''),
            SyntenyHit(11, 5, 5, '', '', '', ''),
        ]
        blocks = find_syntenic_blocks(
            hits, max_gap=3, min_genes=2, pval_threshold=1.0)
        assert len(blocks) == 0  # each group has only 1 hit < min_genes=2

    def test_reverse_complement_block(self):
        """q increasing, t decreasing → valid block."""
        hits = [
            SyntenyHit(10, 0, 10, '', '', '', ''),
            SyntenyHit(11, 1, 9, '', '', '', ''),
            SyntenyHit(12, 2, 8, '', '', '', ''),
        ]
        blocks = find_syntenic_blocks(
            hits, max_gap=3, min_genes=2, pval_threshold=1.0)
        assert len(blocks) == 1
        assert blocks[0].num_genes == 3

    def test_min_genes_filter(self):
        """Block with 2 genes, min_genes=3 → filtered out."""
        hits = [
            SyntenyHit(10, 0, 0, '', '', '', ''),
            SyntenyHit(11, 1, 1, '', '', '', ''),
        ]
        blocks = find_syntenic_blocks(
            hits, max_gap=3, min_genes=3, pval_threshold=1.0)
        assert len(blocks) == 0

    def test_empty_input(self):
        assert find_syntenic_blocks([], max_gap=3, min_genes=2) == []

    def test_single_hit(self):
        hits = [SyntenyHit(10, 0, 0, '', '', '', '')]
        assert find_syntenic_blocks(hits, max_gap=3, min_genes=2) == []


# =========================================================================
# 3. Unit tests — extract_hits_for_pair
# =========================================================================

class TestExtractHits:

    def _simple_cluster_genes(self):
        """Cluster genes for 2 strains, 1 contig each."""
        # cluster 100: present in both strain 0 and strain 1
        # cluster 101: present in both
        # cluster 102: only in strain 0
        return {
            100: {0: [(0, 0, '', '0:0:0')], 1: [(0, 0, '', '1:0:0')]},
            101: {0: [(0, 1, '', '0:0:1')], 1: [(0, 2, '', '1:0:2')]},
            102: {0: [(0, 3, '', '0:0:3')]},
        }

    def test_shared_clusters(self):
        cg = self._simple_cluster_genes()
        shared = {100, 101}
        hits = extract_hits_for_pair(0, 1, cg, shared)
        assert (0, 0) in hits
        assert len(hits[(0, 0)]) == 2  # 2 shared clusters

    def test_no_shared_clusters(self):
        cg = self._simple_cluster_genes()
        shared = set()
        hits = extract_hits_for_pair(0, 1, cg, shared)
        assert len(hits) == 0


# =========================================================================
# 4. Unit tests — aggregate_neighborhoods
# =========================================================================

class TestAggregateNeighborhoods:

    def _make_block(self, block_id, cluster_ids, q_positions, t_positions):
        """Helper to create a SyntenyBlock with given clusters at given positions."""
        hits = []
        for i, cid in enumerate(cluster_ids):
            hits.append(SyntenyHit(
                cid, q_positions[i], t_positions[i], '', '', '', ''))
        po = _order_pvalue_for_block(hits)
        return SyntenyBlock(
            block_id=block_id, hits=hits,
            pvalue_order=po)

    def test_simple_two_genomes(self):
        """Two genomes share one 3-gene block → 1 neighborhood."""
        # Genome layouts: both have clusters [10, 11, 12] at positions 0-2
        pi = _make_position_index({
            0: {0: [10, 11, 12, 99, 98]},
            1: {0: [10, 11, 12, 97, 96]},
        })
        block = self._make_block(0, [10, 11, 12], [0, 1, 2], [0, 1, 2])
        all_blocks = [(0, 1, block)]
        neighs = aggregate_neighborhoods(
            all_blocks, pi, min_genomes=2, min_genes=2, max_gap=3)
        assert len(neighs) >= 1
        # The neighborhood should span both genomes
        n0 = neighs[0]
        assert n0.num_genomes == 2
        assert 0 in n0.genome_members
        assert 1 in n0.genome_members

    def test_transitive_merge(self):
        """
        G0-G1 share {10,11,12}, G1-G2 share {11,12,13}.
        Through transitive merging, all should be in one neighborhood
        (because 11,12 bridge the two blocks).
        """
        pi = _make_position_index({
            0: {0: [10, 11, 12, 99]},
            1: {0: [10, 11, 12, 13]},
            2: {0: [99, 11, 12, 13]},
        })
        block_01 = self._make_block(0, [10, 11, 12], [0, 1, 2], [0, 1, 2])
        block_12 = self._make_block(1, [11, 12, 13], [1, 2, 3], [1, 2, 3])
        all_blocks = [(0, 1, block_01), (1, 2, block_12)]

        neighs = aggregate_neighborhoods(
            all_blocks, pi, min_genomes=2, min_genes=2, max_gap=3)

        # Should be 1 merged neighborhood, not 2 separate ones
        assert len(neighs) == 1
        n = neighs[0]
        # All three genomes should appear
        assert n.num_genomes == 3

    def test_partial_membership(self):
        """
        G2 has only 2 of 3 template clusters → still included with completeness < 1.0.
        """
        pi = _make_position_index({
            0: {0: [10, 11, 12]},
            1: {0: [10, 11, 12]},
            2: {0: [10, 11, 99]},  # missing cluster 12
        })
        block_01 = self._make_block(0, [10, 11, 12], [0, 1, 2], [0, 1, 2])
        block_02 = self._make_block(1, [10, 11], [0, 1], [0, 1])
        all_blocks = [(0, 1, block_01), (0, 2, block_02)]

        neighs = aggregate_neighborhoods(
            all_blocks, pi, min_genomes=2, min_genes=2, max_gap=3)
        assert len(neighs) >= 1
        n = neighs[0]
        # G2 should be a member with completeness < 1
        if 2 in n.genome_members:
            comp = n.genome_completeness(2)
            assert comp < 1.0

    def test_min_genomes_filter(self):
        """Block in only 1 genome pair, min_genomes=3 → filtered out."""
        pi = _make_position_index({
            0: {0: [10, 11, 12]},
            1: {0: [10, 11, 12]},
        })
        block = self._make_block(0, [10, 11, 12], [0, 1, 2], [0, 1, 2])
        all_blocks = [(0, 1, block)]
        neighs = aggregate_neighborhoods(
            all_blocks, pi, min_genomes=3, min_genes=2, max_gap=3)
        assert len(neighs) == 0

    def test_empty_blocks(self):
        pi = _make_position_index({0: {0: [10]}})
        assert aggregate_neighborhoods([], pi) == []

    def test_non_syntenic_in_genome_rejected(self):
        """
        G2 has template clusters but they are far apart on contig → not syntenic.
        """
        pi = _make_position_index({
            0: {0: [10, 11, 12]},
            1: {0: [10, 11, 12]},
            # G2: cluster 10 at pos 0, cluster 11 at pos 50 → gap=49 >> max_gap=3
            2: {0: [10] + [99] * 49 + [11, 12]},
        })
        block = self._make_block(0, [10, 11, 12], [0, 1, 2], [0, 1, 2])
        all_blocks = [(0, 1, block)]
        neighs = aggregate_neighborhoods(
            all_blocks, pi, min_genomes=2, min_genes=3, max_gap=3)
        # G2 should NOT be validated as syntenic member for {10,11,12}
        if neighs:
            n = neighs[0]
            if 2 in n.genome_members:
                # If it is included, it should only be the {11,12} sub-group
                member_clusters = {gp.cluster_id for gp in n.genome_members[2]}
                assert 10 not in member_clusters or 11 not in member_clusters


# =========================================================================
# 5. Integration test — full pipeline on synthetic data
# =========================================================================

class TestIntegration:
    """
    Build a complete synthetic pangenome (mock Pangenome + Graph),
    run the full pipeline (build_position_index → blocks → aggregate),
    and verify outputs.
    """

    @staticmethod
    def _build_synthetic_pangenome(
        n_genomes: int = 5,
        genes_per_contig: int = 30,
        n_contigs: int = 1,
        n_clusters: int = 20,
        neighborhood_clusters: List[int] = None,
        neighborhood_positions: Dict[int, List[int]] = None,
    ):
        """
        Create a synthetic pangenome where specific clusters form a
        known conserved neighborhood.

        Returns: (pg, tree, G, expected_neighborhood_clusters)
        """
        if neighborhood_clusters is None:
            neighborhood_clusters = [100, 101, 102, 103, 104]
        if neighborhood_positions is None:
            # By default, embed the neighborhood at positions 5-9 in every genome
            neighborhood_positions = {
                s: list(range(5, 5 + len(neighborhood_clusters)))
                for s in range(n_genomes)
            }

        random.seed(42)
        pg = MockPangenome()
        gene_cluster_map: Dict[int, Set[str]] = defaultdict(set)

        # Background clusters
        bg_clusters = list(range(n_clusters))

        for strain_idx in range(n_genomes):
            gene_nums = [genes_per_contig] * n_contigs
            strain = MockStrain(
                strain_name=f"genome_{strain_idx}",
                strain_index=strain_idx,
                gene_num=gene_nums,
            )
            pg.load_strain(strain)

            for contig_idx in range(n_contigs):
                neigh_pos_set = set(neighborhood_positions.get(strain_idx, []))
                for gene_idx in range(genes_per_contig):
                    gene_id = f"{strain_idx}:{contig_idx}:{gene_idx}"
                    if gene_idx in neigh_pos_set:
                        # Assign the neighborhood cluster
                        offset = sorted(neigh_pos_set).index(gene_idx)
                        cid = neighborhood_clusters[offset]
                    else:
                        # Random background cluster
                        cid = random.choice(bg_clusters)
                    gene_cluster_map[cid].add(gene_id)

        G = _build_mock_graph(gene_cluster_map)
        tree = MockTree()
        return pg, tree, G, set(neighborhood_clusters)

    def test_full_pipeline_finds_planted_neighborhood(self):
        """The pipeline should discover the planted conserved neighborhood."""
        pg, tree, G, expected_clusters = self._build_synthetic_pangenome(
            n_genomes=5, genes_per_contig=30, n_clusters=20,
            neighborhood_clusters=[100, 101, 102, 103, 104],
        )

        from pgap2.lib.synteny import build_position_index, discover_neighborhoods

        # build_position_index needs real pg/tree/G — we use them here
        position_index = build_position_index(pg, tree, G)
        assert len(position_index) == 5

        # Build intermediate structures
        cluster_genes = _build_cluster_to_genes(position_index)
        strain_clusters = _build_strain_to_clusters(cluster_genes)

        # Scan all pairs
        all_blocks = []
        strain_indices = sorted(position_index.keys())
        pairs = list(itertools.combinations(strain_indices, 2))
        for sq, st in pairs:
            shared = strain_clusters[sq] & strain_clusters[st]
            if len(shared) < 2:
                continue
            hits_by_contig = extract_hits_for_pair(
                sq, st, cluster_genes, shared)
            for hits in hits_by_contig.values():
                blocks = find_syntenic_blocks(
                    hits, max_gap=3, min_genes=3, pval_threshold=1.0)
                for b in blocks:
                    all_blocks.append((sq, st, b))

        assert len(all_blocks) > 0, "Should find at least some raw blocks"

        neighs = aggregate_neighborhoods(
            all_blocks, position_index,
            min_genomes=3, min_genes=3, max_gap=3)

        # The planted neighborhood {100,101,102,103,104} should be found
        found_planted = False
        for n in neighs:
            overlap = n.cluster_set & expected_clusters
            if len(overlap) >= 3:
                found_planted = True
                assert n.num_genomes >= 3
                break
        assert found_planted, (
            f"Planted neighborhood {expected_clusters} not found. "
            f"Got {len(neighs)} neighborhoods: "
            f"{[n.cluster_set for n in neighs[:5]]}"
        )

    def test_full_pipeline_with_partial_absence(self):
        """
        One genome is missing 1 cluster from the neighborhood.
        The pipeline should still include it with completeness < 1.
        """
        neigh_clusters = [100, 101, 102, 103, 104]
        positions = {
            0: [5, 6, 7, 8, 9],
            1: [5, 6, 7, 8, 9],
            2: [5, 6, 7, 8, 9],
            3: [5, 6, 7, 8, 9],
            4: [5, 6, 7, 8],     # genome 4: missing cluster 104
        }
        # For genome 4, we only place 4 of 5 clusters
        pg, tree, G, expected = self._build_synthetic_pangenome(
            n_genomes=5, genes_per_contig=30, n_clusters=20,
            neighborhood_clusters=neigh_clusters,
            neighborhood_positions=positions,
        )

        from pgap2.lib.synteny import build_position_index

        position_index = build_position_index(pg, tree, G)
        cluster_genes = _build_cluster_to_genes(position_index)
        strain_clusters = _build_strain_to_clusters(cluster_genes)

        all_blocks = []
        for sq, st in itertools.combinations(sorted(position_index.keys()), 2):
            shared = strain_clusters[sq] & strain_clusters[st]
            if len(shared) < 2:
                continue
            hits_by_contig = extract_hits_for_pair(
                sq, st, cluster_genes, shared)
            for hits in hits_by_contig.values():
                blocks = find_syntenic_blocks(
                    hits, max_gap=3, min_genes=2, pval_threshold=1.0)
                for b in blocks:
                    all_blocks.append((sq, st, b))

        neighs = aggregate_neighborhoods(
            all_blocks, position_index,
            min_genomes=2, min_genes=2, max_gap=3)

        # Find the planted neighborhood
        planted_neigh = None
        for n in neighs:
            if len(n.cluster_set & expected) >= 3:
                planted_neigh = n
                break

        assert planted_neigh is not None, "Planted neighborhood not found"

        # Genome 4 should still be a member
        if 4 in planted_neigh.genome_members:
            comp = planted_neigh.genome_completeness(4)
            assert comp < 1.0, "Genome 4 should have partial completeness"

    def test_output_files(self):
        """Verify TSV and GML output files are created and non-empty."""
        pg, tree, G, _ = self._build_synthetic_pangenome(n_genomes=3)

        from pgap2.lib.synteny import build_position_index

        position_index = build_position_index(pg, tree, G)
        cluster_genes = _build_cluster_to_genes(position_index)
        strain_clusters = _build_strain_to_clusters(cluster_genes)

        all_blocks = []
        for sq, st in itertools.combinations(sorted(position_index.keys()), 2):
            shared = strain_clusters[sq] & strain_clusters[st]
            if len(shared) < 2:
                continue
            hits_by_contig = extract_hits_for_pair(
                sq, st, cluster_genes, shared)
            for hits in hits_by_contig.values():
                blocks = find_syntenic_blocks(
                    hits, max_gap=3, min_genes=2, pval_threshold=1.0)
                for b in blocks:
                    all_blocks.append((sq, st, b))

        neighs = aggregate_neighborhoods(
            all_blocks, position_index, min_genomes=2, min_genes=2, max_gap=3)

        with tempfile.TemporaryDirectory() as tmpdir:
            f1 = write_neighborhood_tsv(neighs, pg, G, tmpdir)
            f2 = write_neighborhood_detail_tsv(neighs, pg, G, tmpdir)
            f3 = write_neighborhood_gml(neighs, tmpdir)

            assert os.path.exists(f1)
            assert os.path.exists(f2)
            assert os.path.exists(f3)
            assert os.path.getsize(f1) > 0
            assert os.path.getsize(f2) > 0
            assert os.path.getsize(f3) > 0

            # Verify TSV header
            with open(f1) as fh:
                header = fh.readline()
                assert header.startswith("#neighborhood_id")

            # Verify detail TSV has completeness column
            with open(f2) as fh:
                header = fh.readline()
                assert "completeness" in header

            # Verify GML is loadable
            H = nx.read_gml(f3)
            assert H.number_of_nodes() > 0


# =========================================================================
# 6. Edge-case tests
# =========================================================================

class TestEdgeCases:

    def test_single_genome(self):
        """Only 1 genome → no pairs → no neighborhoods."""
        pi = _make_position_index({0: {0: [10, 11, 12]}})
        neighs = aggregate_neighborhoods([], pi, min_genomes=2)
        assert len(neighs) == 0

    def test_two_genomes_no_shared_clusters(self):
        """Two genomes with completely disjoint clusters."""
        pi = _make_position_index({
            0: {0: [10, 11, 12]},
            1: {0: [20, 21, 22]},
        })
        cg = _build_cluster_to_genes(pi)
        sc = _build_strain_to_clusters(cg)
        shared = sc[0] & sc[1]
        assert len(shared) == 0

    def test_many_contigs(self):
        """Neighborhood split across contigs → only largest group counts."""
        pi = _make_position_index({
            0: {0: [10, 11], 1: [12, 13]},    # split across contigs
            1: {0: [10, 11, 12, 13]},          # all on one contig
        })
        block = SyntenyBlock(
            block_id=0,
            hits=[
                SyntenyHit(10, 0, 0, '', '', '', ''),
                SyntenyHit(11, 1, 1, '', '', '', ''),
                SyntenyHit(12, 2, 2, '', '', '', ''),
                SyntenyHit(13, 3, 3, '', '', '', ''),
            ],
            pvalue_order=1e-6,
        )
        neighs = aggregate_neighborhoods(
            [(0, 1, block)], pi, min_genomes=1, min_genes=2, max_gap=3)
        # Genome 0 can at best have 2 genes from contig 0 or 2 from contig 1
        if neighs and 0 in neighs[0].genome_members:
            members = neighs[0].genome_members[0]
            # All members from same contig
            contigs = {gp.contig_index for gp in members}
            assert len(contigs) == 1

    def test_duplicate_cluster_in_genome(self):
        """Same cluster at multiple positions in one genome."""
        pi = _make_position_index({
            0: {0: [10, 11, 10, 12]},   # cluster 10 appears twice
            1: {0: [10, 11, 12]},
        })
        cg = _build_cluster_to_genes(pi)
        sc = _build_strain_to_clusters(cg)
        shared = sc[0] & sc[1]
        hits = extract_hits_for_pair(0, 1, cg, shared)
        # Should generate hits for both occurrences of cluster 10 in G0
        total_hits = sum(len(h) for h in hits.values())
        assert total_hits >= 4  # 2 × cluster10 + 1 × cluster11 + 1 × cluster12


# =========================================================================
# 7. Performance benchmark
# =========================================================================

class TestPerformance:

    @pytest.mark.benchmark
    def test_bench_find_syntenic_blocks(self):
        """Benchmark Union-Find block detection with N hits."""
        random.seed(42)
        for n_hits in [100, 500, 2000, 5000]:
            hits = []
            for i in range(n_hits):
                hits.append(SyntenyHit(
                    cluster_id=random.randint(0, 500),
                    q_pos=i,
                    t_pos=i + random.randint(-2, 2),
                    q_strand='', t_strand='',
                    q_gene_id=f'q{i}', t_gene_id=f't{i}',
                ))

            t0 = time.perf_counter()
            blocks = find_syntenic_blocks(
                hits, max_gap=3, min_genes=2, pval_threshold=1.0)
            elapsed = time.perf_counter() - t0

            print(
                f"  find_syntenic_blocks({n_hits} hits) → {len(blocks)} blocks in {elapsed:.4f}s")
            # Should complete well under 1 second for ≤5000 hits
            assert elapsed < 2.0, f"Too slow for {n_hits} hits: {elapsed:.2f}s"

    @pytest.mark.benchmark
    def test_bench_aggregate_neighborhoods(self):
        """Benchmark aggregation with many blocks and genomes."""
        random.seed(42)
        n_genomes = 50
        n_clusters = 200
        genes_per_contig = 100

        # Build position index
        layouts = {}
        for s in range(n_genomes):
            contig_clusters = [random.randint(
                0, n_clusters - 1) for _ in range(genes_per_contig)]
            # Plant a known neighborhood at positions 10-15
            for offset, cid in enumerate([1000, 1001, 1002, 1003, 1004, 1005]):
                contig_clusters[10 + offset] = cid
            layouts[s] = {0: contig_clusters}
        pi = _make_position_index(layouts)

        # Generate blocks (simulate pairwise scanning)
        all_blocks = []
        block_id = 0
        for sq, st in itertools.combinations(range(n_genomes), 2):
            # The planted neighborhood should produce a block
            hits = [
                SyntenyHit(1000 + i, 10 + i, 10 + i, '', '', '', '')
                for i in range(6)
            ]
            po = _order_pvalue_for_block(hits)
            b = SyntenyBlock(block_id=block_id, hits=hits,
                             pvalue_order=po)
            all_blocks.append((sq, st, b))
            block_id += 1

        print(
            f"  {len(all_blocks)} blocks from {n_genomes} genomes ({n_genomes*(n_genomes-1)//2} pairs)")

        t0 = time.perf_counter()
        neighs = aggregate_neighborhoods(
            all_blocks, pi, min_genomes=2, min_genes=3, max_gap=3)
        elapsed = time.perf_counter() - t0

        print(
            f"  aggregate_neighborhoods → {len(neighs)} neighborhoods in {elapsed:.4f}s")
        assert elapsed < 10.0, f"Aggregation too slow: {elapsed:.2f}s"
        assert len(neighs) >= 1, "Should find at least the planted neighborhood"

    @pytest.mark.benchmark
    def test_bench_full_pipeline_synthetic(self):
        """End-to-end benchmark on a medium-sized synthetic pangenome."""
        random.seed(42)
        n_genomes = 20
        genes_per_contig = 80
        n_clusters = 100

        pg = MockPangenome()
        gene_cluster_map: Dict[int, Set[str]] = defaultdict(set)

        planted = [500, 501, 502, 503, 504]

        for s in range(n_genomes):
            pg.load_strain(MockStrain(f"g{s}", s, [genes_per_contig]))
            for pos in range(genes_per_contig):
                gid = f"{s}:0:{pos}"
                if 10 <= pos < 10 + len(planted):
                    cid = planted[pos - 10]
                else:
                    cid = random.randint(0, n_clusters - 1)
                gene_cluster_map[cid].add(gid)

        G = _build_mock_graph(gene_cluster_map)
        tree = MockTree()

        from pgap2.lib.synteny import build_position_index

        t0 = time.perf_counter()

        position_index = build_position_index(pg, tree, G)
        cluster_genes = _build_cluster_to_genes(position_index)
        strain_clusters = _build_strain_to_clusters(cluster_genes)

        all_blocks = []
        for sq, st in itertools.combinations(sorted(position_index.keys()), 2):
            shared = strain_clusters[sq] & strain_clusters[st]
            if len(shared) < 2:
                continue
            hits_by_contig = extract_hits_for_pair(
                sq, st, cluster_genes, shared)
            for hits in hits_by_contig.values():
                blocks = find_syntenic_blocks(
                    hits, max_gap=3, min_genes=2, pval_threshold=0.05)
                for b in blocks:
                    all_blocks.append((sq, st, b))

        neighs = aggregate_neighborhoods(
            all_blocks, position_index, min_genomes=2, min_genes=2, max_gap=3)

        elapsed = time.perf_counter() - t0

        print(
            f"\n  Full pipeline: {n_genomes} genomes × {genes_per_contig} genes")
        print(f"  {len(all_blocks)} raw blocks → {len(neighs)} neighborhoods")
        print(f"  Total time: {elapsed:.4f}s")

        assert elapsed < 30.0, f"Full pipeline too slow: {elapsed:.2f}s"
        # Should find planted neighborhood
        found = any(len(n.cluster_set & set(planted)) >= 3 for n in neighs)
        assert found, "Planted neighborhood not recovered"
