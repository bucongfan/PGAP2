"""
Pangenome graph post-processing: annotated GML generation + reference layout.

``pgap2 post graph -i <dir> -o <dir> [--ref <strain>]``

Pipeline
--------
1. **Annotated GML** — Read raw ``pgap2.partition.map.gml`` together with
   ``total.involved_annot.tsv`` and ``pgap2.partition.gene_content.detail.tsv``
   to produce ``pgap2.partition.map_annotated.gml`` with rich per-node
   attributes (gene_name, product, group, genomeIDs, geneIDs, …).

2. **Reference-based layout** — BFS + min-cut algorithm (see
   :mod:`pgap2.postprocess.reference_layout`) to detect and sever long-range
   edges.  Produces ``pgap2.partition.map_annotated_layout.gml`` (cut edges
   removed) and ``pgap2.partition.map_annotated_complete.gml`` (all edges with
   ``is_cut`` attribute).

Outputs
-------
- ``pgap2.partition.map_annotated.gml``           — enriched pangenome graph
- ``pgap2.partition.map_annotated_layout.gml``     — long-range edges removed
- ``pgap2.partition.map_annotated_complete.gml``   — all edges, ``is_cut`` flag
- ``postprocess_annotate.log``                     — log file

Usage
-----
    pgap2 post graph -i <partition_output_dir> [-o <out_dir>] [--ref <strain>]
"""

import os
import re
import json
import math
import argparse
import collections
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
from loguru import logger
from argparse import ArgumentParser, _SubParsersAction

from pgap2.utils.supply import set_verbosity_level


def parse_gml(gml_path: str) -> Tuple[Dict[int, str], List[Tuple[int, int]]]:
    """Parse a PGAP2-flavour GML file.

    Returns
    -------
    nodes : dict[int, str]
        Mapping of node-id → label string.
    edges : list of (source_id, target_id)
    """
    nodes: Dict[int, str] = {}
    edges: List[Tuple[int, int]] = []

    with open(gml_path, "r", encoding="utf-8") as fh:
        in_node = False
        in_edge = False
        cur_id: Optional[int] = None
        cur_label: Optional[str] = None
        cur_src: Optional[int] = None
        cur_tgt: Optional[int] = None

        for raw in fh:
            line = raw.strip()
            if line == "node [":
                in_node, cur_id, cur_label = True, None, None
            elif line == "edge [":
                in_edge, cur_src, cur_tgt = True, None, None
            elif line == "]":
                if in_node and cur_id is not None:
                    nodes[cur_id] = cur_label if cur_label else str(cur_id)
                elif in_edge and cur_src is not None and cur_tgt is not None:
                    edges.append((cur_src, cur_tgt))
                in_node = in_edge = False
            elif in_node:
                if line.startswith("id "):
                    cur_id = int(line[3:].strip())
                elif line.startswith("label "):
                    cur_label = line[6:].strip().strip('"')
            elif in_edge:
                if line.startswith("source "):
                    cur_src = int(line[7:].strip())
                elif line.startswith("target "):
                    cur_tgt = int(line[7:].strip())

    return nodes, edges


# ---------------------------------------------------------------------------
# Annotation loader
# ---------------------------------------------------------------------------

def load_annot(annot_path: str) -> Dict[str, dict]:
    """Load ``total.involved_annot.tsv`` into a dict keyed by Gene_index.

    Each value is a dict with keys: strain, contig, start, end, strand,
    length, gene_id, gene_name, product_name.
    """
    loc_re = re.compile(r"\[(\d+):(\d+)\]\(([+-])\)")
    annot: Dict[str, dict] = {}

    with open(annot_path, "r", encoding="utf-8") as fh:
        header = fh.readline()  # skip header
        for raw in fh:
            parts = raw.split("\t")
            if len(parts) < 8:
                continue
            gene_index = parts[0]
            strain = parts[1]
            contig = parts[2]
            loc_str = parts[3]
            length = int(parts[4]) if parts[4].isdigit() else 0
            gene_id = parts[5]
            gene_name = parts[6]
            product_name = parts[7]

            m = loc_re.match(loc_str)
            start, end, strand = (int(m.group(1)), int(
                m.group(2)), m.group(3)) if m else (0, 0, "+")

            annot[gene_index] = {
                "strain": strain,
                "contig": contig,
                "start": start,
                "end": end,
                "strand": strand,
                "length": length,
                "gene_id": gene_id,
                "gene_name": gene_name,
                "product_name": product_name,
            }
    return annot


# ---------------------------------------------------------------------------
# Build strain_idx → strain_name mapping from annotation
# ---------------------------------------------------------------------------

def build_strain_map(annot: Dict[str, dict]) -> Dict[int, str]:
    """Return {strain_idx: strain_name} derived from Gene_index keys.

    Gene_index format: ``strain_idx:copy_idx:gene_order``
    """
    strain_map: Dict[int, str] = {}
    for gene_index, info in annot.items():
        parts = gene_index.split(":")
        if len(parts) >= 3:
            s_idx = int(parts[0])
            if s_idx not in strain_map:
                strain_map[s_idx] = info["strain"]
    return strain_map


# ---------------------------------------------------------------------------
# Expand compound labels
# ---------------------------------------------------------------------------

def expand_label(label: str) -> List[str]:
    """Split a possibly compound label ``a:b:c_d:e:f`` into individual
    gene-index strings ``['a:b:c', 'd:e:f']``.
    """
    if "_" not in label:
        return [label]
    # Each gene index has the form <int>:<int>:<int>.
    # Compound labels join them with ``_``.
    parts = label.split("_")
    result = []
    buf: List[str] = []
    for p in parts:
        buf.append(p)
        joined = "_".join(buf)
        # Check if it looks like a valid gene index (3 colon-separated ints)
        segs = joined.split(":")
        if len(segs) == 3 and all(s.lstrip("-").isdigit() for s in segs):
            result.append(joined)
            buf = []
    # If there's leftover (shouldn't happen), just join it back
    if buf:
        result.append("_".join(buf))
    return result


# ---------------------------------------------------------------------------
# Core layout algorithm
# ---------------------------------------------------------------------------

def compute_layout(
    gml_path: str,
    annot_path: str,
    ref_strain: Optional[str] = None,
    threshold: Optional[int] = None,
    max_iterations: int = 50,
) -> dict:
    """Compute reference-based linear layout for a pangenome graph.

    Parameters
    ----------
    gml_path : str
        Path to ``pgap2.partition.map.gml``.
    annot_path : str
        Path to ``total.involved_annot.tsv``.
    ref_strain : str or None
        Name of the reference strain.  If *None*, the strain with the most
        genes is used.
    threshold : int or None
        Distance threshold for classifying long-range edges.  If *None*,
        ``max(20, int(genome_length * 0.05))`` is used.
    max_iterations : int
        Max iterations for Laplacian smoothing of non-ref node coordinates.

    Returns
    -------
    dict
        JSON-serialisable result with keys ``nodes``, ``edges``, ``stats``,
        ``ref_genome``, ``available_refs``.
    """
    # 1. Load data
    logger.info("Loading GML …")
    nodes_raw, edges_raw = parse_gml(gml_path)
    logger.info(f"  {len(nodes_raw)} nodes, {len(edges_raw)} edges")

    logger.info("Loading annotations …")
    annot = load_annot(annot_path)
    logger.info(f"  {len(annot)} annotation entries")

    strain_map = build_strain_map(annot)
    available_refs = sorted(set(strain_map.values()))
    logger.info(f"  {len(available_refs)} strains available")

    # 2. Pick reference strain
    if ref_strain is None:
        # Count genes per strain
        strain_gene_count: Dict[str, int] = collections.Counter()
        for info in annot.values():
            strain_gene_count[info["strain"]] += 1
        ref_strain = strain_gene_count.most_common(1)[0][0]
        logger.info(f"  Auto-selected reference: {ref_strain}")
    else:
        if ref_strain not in available_refs:
            raise ValueError(
                f"Reference strain '{ref_strain}' not found. "
                f"Available: {available_refs}"
            )
        logger.info(f"  Using reference: {ref_strain}")

    # Find the strain_idx for the reference
    ref_strain_idx: Optional[int] = None
    for idx, name in strain_map.items():
        if name == ref_strain:
            ref_strain_idx = idx
            break

    # 3. Build adjacency
    adj: Dict[int, Set[int]] = collections.defaultdict(set)
    for s, t in edges_raw:
        adj[s].add(t)
        adj[t].add(s)

    # 4. Identify reference nodes and assign ref_pos from gene_order
    #    A node belongs to the reference if *any* of its gene indices has
    #    strain_idx == ref_strain_idx.
    ref_node_pos: Dict[int, float] = {}  # node_id → ref_pos (gene_order)
    node_gene_indices: Dict[int, List[str]] = {}  # node_id → [gene_index, ...]

    for nid, label in nodes_raw.items():
        indices = expand_label(label)
        node_gene_indices[nid] = indices
        for gi in indices:
            parts = gi.split(":")
            if len(parts) == 3 and int(parts[0]) == ref_strain_idx:
                gene_order = int(parts[2])
                # If multiple ref genes in same node, take mean position
                if nid in ref_node_pos:
                    ref_node_pos[nid] = (ref_node_pos[nid] + gene_order) / 2.0
                else:
                    ref_node_pos[nid] = float(gene_order)

    genome_length = max(ref_node_pos.values()) + 1 if ref_node_pos else 1
    logger.info(
        f"  Reference nodes: {len(ref_node_pos)}, genome length: {int(genome_length)}")

    # 5. Adaptive threshold
    if threshold is None:
        threshold = max(20, int(genome_length * 0.05))
    logger.info(f"  Long-range threshold: {threshold}")

    # 6. Propagate coordinates to non-reference nodes via Laplacian smoothing
    #    pos[nid] for ref nodes is fixed; for non-ref nodes, iteratively set to
    #    the mean of neighbour positions.
    pos: Dict[int, float] = dict(ref_node_pos)  # start with ref positions

    # Initialize non-ref nodes: mean of ref-neighbours, or NaN
    non_ref_ids = [nid for nid in nodes_raw if nid not in ref_node_pos]

    for nid in non_ref_ids:
        nbr_positions = [ref_node_pos[nb]
                         for nb in adj[nid] if nb in ref_node_pos]
        if nbr_positions:
            pos[nid] = float(np.mean(nbr_positions))

    # Iterative smoothing
    for iteration in range(max_iterations):
        max_delta = 0.0
        for nid in non_ref_ids:
            nbr_positions = [pos[nb] for nb in adj[nid] if nb in pos]
            if not nbr_positions:
                continue
            new_pos = float(np.mean(nbr_positions))
            if nid in pos:
                max_delta = max(max_delta, abs(new_pos - pos[nid]))
            pos[nid] = new_pos
        if max_delta < 0.5:
            logger.info(f"  Laplacian converged at iteration {iteration + 1}")
            break

    # Nodes that still have no position (isolated from reference): assign -1
    for nid in nodes_raw:
        if nid not in pos:
            pos[nid] = -1.0

    # 7. Count genomes per node
    node_num_genomes: Dict[int, int] = {}
    for nid, indices in node_gene_indices.items():
        strains = set()
        for gi in indices:
            if gi in annot:
                strains.add(annot[gi]["strain"])
            else:
                parts = gi.split(":")
                if len(parts) == 3:
                    s_idx = int(parts[0])
                    if s_idx in strain_map:
                        strains.add(strain_map[s_idx])
        node_num_genomes[nid] = len(strains)

    # 8. Classify edges
    backbone_edges = []
    local_edges = []
    longrange_edges = []
    unplaced_edges = []

    for s, t in edges_raw:
        sp = pos.get(s, -1)
        tp = pos.get(t, -1)

        if sp < 0 or tp < 0:
            unplaced_edges.append((s, t, -1.0))
            continue

        # Circular distance
        raw_dist = abs(sp - tp)
        dist = min(raw_dist, genome_length -
                   raw_dist) if genome_length > 1 else raw_dist

        is_s_ref = s in ref_node_pos
        is_t_ref = t in ref_node_pos

        if is_s_ref and is_t_ref and dist <= 1.5:
            backbone_edges.append((s, t, dist))
        elif dist <= threshold:
            local_edges.append((s, t, dist))
        else:
            longrange_edges.append((s, t, dist))

    logger.info(
        f"  Edges: {len(backbone_edges)} backbone, {len(local_edges)} local, "
        f"{len(longrange_edges)} long-range, {len(unplaced_edges)} unplaced"
    )

    # 9. Build annotation lookup for each node (compact)
    def _clean_annot_str(s: str) -> str:
        """Strip surrounding list brackets/quotes, e.g. \"['dnaA']\" → \"dnaA\"."""
        s = s.strip()
        if s.startswith("[") and s.endswith("]"):
            s = s[1:-1].strip()
        if s.startswith("'") and s.endswith("'"):
            s = s[1:-1].strip()
        if s.startswith('"') and s.endswith('"'):
            s = s[1:-1].strip()
        return s

    def node_annotation(nid: int) -> str:
        """Return a short annotation string for the node."""
        for gi in node_gene_indices.get(nid, []):
            if gi in annot:
                info = annot[gi]
                gn = _clean_annot_str(info["gene_name"])
                pn = _clean_annot_str(info["product_name"])
                if gn and gn != "-":
                    return gn
                if pn and pn != "-":
                    return pn[:60]
        return ""

    # 10. Compute circular 2D coordinates
    #     Map ref_pos → angle around a circle so the graph looks like a
    #     circular bacterial genome (similar to Panaroo reference_based_layout.py).
    R_BASE = 500.0          # base radius for the circle
    R_REF_OFFSET = 0.0      # ref nodes sit exactly on the circle
    R_ACC_OFFSET = 30.0     # accessory nodes sit slightly inside/outside
    TWO_PI = 2.0 * math.pi

    node_x: Dict[int, float] = {}
    node_y: Dict[int, float] = {}

    for nid in nodes_raw:
        p = pos.get(nid, -1.0)
        if p < 0:
            # Unplaced node – put in the center
            node_x[nid] = 0.0
            node_y[nid] = 0.0
            continue

        angle = TWO_PI * (p / genome_length) if genome_length > 1 else 0.0
        is_ref = nid in ref_node_pos
        # Accessory nodes get a small radial offset based on their id
        # to avoid overlapping; use a deterministic jitter
        if is_ref:
            r = R_BASE
        else:
            # Offset based on degree and num_genomes for visual separation
            jitter = ((nid * 7) % 31 - 15) / 15.0  # ∈ [-1, 1]
            r = R_BASE + R_ACC_OFFSET * jitter

        node_x[nid] = round(r * math.cos(angle), 2)
        node_y[nid] = round(r * math.sin(angle), 2)

    logger.info("  Circular coordinates computed")

    # 11. Assemble output
    #     Long-range edges are REMOVED from the output (like Panaroo) so
    #     the circular layout is clean.
    out_nodes = []
    for nid in sorted(nodes_raw.keys()):
        out_nodes.append({
            "id": nid,
            "label": nodes_raw[nid],
            "ref_pos": round(pos.get(nid, -1.0), 2),
            "is_ref": nid in ref_node_pos,
            "num_genomes": node_num_genomes.get(nid, 0),
            "degree": len(adj[nid]),
            "annotation": node_annotation(nid),
            "x": node_x[nid],
            "y": node_y[nid],
        })

    def fmt_edges(edge_list, etype):
        return [
            {"source": s, "target": t, "type": etype, "distance": round(d, 2)}
            for s, t, d in edge_list
        ]

    # Only keep backbone + local (+ unplaced); longrange are cut
    out_edges = (
        fmt_edges(backbone_edges, "backbone")
        + fmt_edges(local_edges, "local")
        + fmt_edges(unplaced_edges, "unplaced")
    )

    stats = {
        "total_nodes": len(nodes_raw),
        "total_edges": len(out_edges),
        "ref_nodes": len(ref_node_pos),
        "non_ref_nodes": len(non_ref_ids),
        "backbone_edges": len(backbone_edges),
        "local_edges": len(local_edges),
        "longrange_edges": len(longrange_edges),   # how many were cut
        "unplaced_edges": len(unplaced_edges),
        "genome_length": int(genome_length),
        "threshold": threshold,
    }

    # Long-range edges are recorded separately (for reference/download)
    cut_edges = [
        {"source": s, "target": t, "distance": round(d, 2)}
        for s, t, d in longrange_edges
    ]

    return {
        "nodes": out_nodes,
        "edges": out_edges,
        "stats": stats,
        "ref_genome": ref_strain,
        "available_refs": available_refs,
        "cut_edges": cut_edges,
    }


def write_gml(
    result: dict,
    original_nodes: Dict[int, str],
    out_path: str,
):
    """Write a processed GML with long-range edges removed and
    reference coordinates / circular x,y embedded in node attributes.

    This GML can be opened directly in Cytoscape Desktop.
    """
    node_lookup = {n["id"]: n for n in result["nodes"]}

    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("graph [\n  directed 0\n")

        for nid in sorted(original_nodes.keys()):
            info = node_lookup.get(nid, {})
            label = original_nodes[nid]
            fh.write("  node [\n")
            fh.write(f"    id {nid}\n")
            fh.write(f'    label "{label}"\n')
            fh.write(f"    ref_pos {info.get('ref_pos', -1.0):.2f}\n")
            fh.write(f"    is_ref {1 if info.get('is_ref') else 0}\n")
            fh.write(f"    num_genomes {info.get('num_genomes', 0)}\n")
            annot_str = info.get('annotation', '').replace('"', "'")
            fh.write(f'    annotation "{annot_str}"\n')
            fh.write(f"    x {info.get('x', 0.0):.2f}\n")
            fh.write(f"    y {info.get('y', 0.0):.2f}\n")
            fh.write("  ]\n")

        # Only write kept edges (backbone + local + unplaced)
        for e in result["edges"]:
            fh.write("  edge [\n")
            fh.write(f"    source {e['source']}\n")
            fh.write(f"    target {e['target']}\n")
            fh.write(f'    type "{e["type"]}"\n')
            fh.write(f"    distance {e['distance']:.2f}\n")
            fh.write("  ]\n")

        fh.write("]\n")

    logger.info(f"  GML written to {out_path}")


def _parse_detail_tsv(detail_path: str) -> Tuple[List[str], List[dict]]:
    """Parse ``pgap2.partition.gene_content.detail.tsv``.

    Returns
    -------
    strain_names : list[str]
        Ordered strain names from the header (e.g. ``['EC01', 'EC02', …]``).
    clusters : list[dict]
        One dict per cluster with keys:
        ``clust_id``, ``gene_name``, ``product``, ``group``, ``repre_gene``,
        ``min``, ``mean``, ``var``, ``uni``, ``involved_strain``,
        ``para_strain``, ``involved_gene``, ``para_gene``,
        ``per_strain_genes`` (list[str] – comma-separated gene indices).
    """
    clusters: List[dict] = []
    strain_names: List[str] = []

    with open(detail_path, "r", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n")
        # Header: #Clust\tgene_name\tproduct\tgroup\trepre_gene\tmin\tmean\tvar\tuni
        #         \tinvolved_strain\tpara_strain\tinvolved_gene\tpara_gene\tEC01,EC02,...
        hdr_parts = header.split("\t")
        if len(hdr_parts) >= 14:
            strain_names = hdr_parts[13].split(",")

        for raw in fh:
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 14:
                continue
            per_strain = parts[13].split(",")   # gene indices per strain
            clusters.append({
                "clust_id": parts[0],           # e.g. "clust_0"
                "gene_name": parts[1],
                "product": parts[2],
                "group": parts[3],
                "repre_gene": parts[4],
                "min": parts[5],
                "mean": parts[6],
                "var": parts[7],
                "uni": parts[8],
                "involved_strain": int(parts[9]) if parts[9].isdigit() else 0,
                "para_strain": int(parts[10]) if parts[10].isdigit() else 0,
                "involved_gene": int(parts[11]) if parts[11].isdigit() else 0,
                "para_gene": int(parts[12]) if parts[12].isdigit() else 0,
                "per_strain_genes": per_strain,
            })

    return strain_names, clusters


def _build_gene_to_cluster_map(clusters: List[dict]) -> Dict[str, int]:
    """Build a mapping from individual gene_index → cluster list index.

    Each cluster's ``per_strain_genes`` may contain compound entries
    (``a:b:c;d:e:f`` for paralogs within the same strain).  We split on
    ``;`` and register every single gene_index.
    """
    gene2clust: Dict[str, int] = {}
    for idx, cl in enumerate(clusters):
        for cell in cl["per_strain_genes"]:
            if not cell:
                continue
            for gene_index in cell.split(";"):
                gene_index = gene_index.strip()
                if gene_index:
                    gene2clust[gene_index] = idx
    return gene2clust


def _gml_escape(s: str) -> str:
    """Escape a string for GML attribute values."""
    return s.replace("\\", "\\\\").replace('"', '\\"')


def _clean_bracket_str(s: str) -> str:
    """Clean Python list-repr strings that may be multi-valued and
    semicolon-separated.

    Examples::

        "['dnaA']"                       → "dnaA"
        "[]"                             → ""
        "['rsmI'];"                      → "rsmI"
        "['a'];['b']"                    → "a; b"
        "['SAM-dep'];[\"16S rRNA ...\"]" → "SAM-dep; 16S rRNA ..."
    """
    s = s.strip()
    if not s:
        return ""

    # Split on ';' that separates bracket groups, but NOT on ';' inside
    # quoted strings.  A simple heuristic: split on '];' or '];[' boundaries.
    import re
    # Split on '];' optionally followed by '[' to separate groups
    parts = re.split(r"\]\s*;\s*\[?", s)

    cleaned: list = []
    for part in parts:
        p = part.strip()
        # Remove leading '[' or trailing ']'
        if p.startswith("["):
            p = p[1:]
        if p.endswith("]"):
            p = p[:-1]
        # Remove trailing ';'
        p = p.rstrip(";").strip()
        if not p:
            continue
        # Remove surrounding quotes
        if (p.startswith("'") and p.endswith("'")) or \
           (p.startswith('"') and p.endswith('"')):
            p = p[1:-1]
        p = p.strip()
        if p:
            cleaned.append(p)
    return "; ".join(cleaned)


def generate_annotated_gml(
    indir: str,
    outdir: str,
):
    """Read PGAP2 output files and write a richly annotated GML.

    The output GML has:
    - Nodes renamed to ``clust_<id>``
    - Rich annotation attributes on every node
    - Edge attributes with ``size`` (# shared strains)

    Parameters
    ----------
    indir : str
        Directory containing PGAP2 partition outputs.
    outdir : str
        Directory where ``pgap2.partition.map_annotated.gml`` will be written.
    """
    gml_path = os.path.join(indir, "pgap2.partition.map.gml")
    detail_path = os.path.join(
        indir, "pgap2.partition.gene_content.detail.tsv")
    annot_path = os.path.join(indir, "total.involved_annot.tsv")

    for fp, desc in [
        (gml_path, "GML file"),
        (detail_path, "Detail TSV"),
        (annot_path, "Annotation file"),
    ]:
        if not os.path.exists(fp):
            raise FileNotFoundError(f"{desc} not found: {fp}")

    # 1. Parse the raw GML
    logger.info("Loading GML …")
    nodes_raw, edges_raw = parse_gml(gml_path)
    logger.info(f"  {len(nodes_raw)} nodes, {len(edges_raw)} edges")

    # 2. Parse detail.tsv (cluster metadata + per-strain gene lists)
    logger.info("Loading detail.tsv …")
    strain_names, clusters = _parse_detail_tsv(detail_path)
    logger.info(f"  {len(clusters)} clusters, {len(strain_names)} strains")

    # 3. Load annotation for gene_id lookup
    annot = load_annot(annot_path)

    # 4. Map each GML node → cluster
    #    PGAP2's partition.py iterates G.nodes() and writes one cluster per
    #    node (clust_0, clust_1, …).  The subsequent nx.write_gml preserves
    #    the same ordering, so GML node id=i corresponds to clust_i.
    n_mapped = min(len(nodes_raw), len(clusters))
    n_unmapped = max(0, len(nodes_raw) - len(clusters))
    logger.info(
        f"  Direct id→cluster mapping: {n_mapped} mapped, "
        f"{n_unmapped} unmapped"
    )

    # 5. Build per-node enrichment
    #    For each cluster, collect:
    #      - genomeIDs: semicolon-separated strain names that have this gene
    #      - geneIDs:   semicolon-separated gene_index strings for all copies
    #      - members:   semicolon-separated gene_id (e.g. 1_1_1) per copy
    #      - size:      total number of gene copies (= involved_gene)
    #      - paralogs:  1 if para_strain > 0, else 0

    def _collect_cluster_attrs(cl: dict) -> dict:
        """Derive enrichment attributes for one cluster."""
        genome_ids = []
        gene_ids = []
        member_ids = []

        for i, cell in enumerate(cl["per_strain_genes"]):
            if not cell:
                continue
            strain_name = strain_names[i] if i < len(strain_names) else str(i)
            for gi in cell.split(";"):
                gi = gi.strip()
                if not gi:
                    continue
                gene_ids.append(gi)
                genome_ids.append(strain_name)
                # Look up the Gene_ID (e.g. "1_1_1") from annotation
                if gi in annot:
                    member_ids.append(annot[gi]["gene_id"])
                else:
                    member_ids.append(gi)

        return {
            "genomeIDs": ";".join(genome_ids),
            "geneIDs": ";".join(gene_ids),
            "members": ";".join(member_ids),
            "size": cl["involved_gene"],
        }

    cluster_attrs_cache: Dict[int, dict] = {}
    for idx, cl in enumerate(clusters):
        cluster_attrs_cache[idx] = _collect_cluster_attrs(cl)

    # 6. Build per-node strain sets for edge capacity computation
    node_strain_set: Dict[int, Set[str]] = {}
    for nid in nodes_raw:
        cidx = nid  # direct mapping: node id = cluster index
        if cidx < len(clusters):
            cl = clusters[cidx]
            strains = set()
            for i, cell in enumerate(cl["per_strain_genes"]):
                if cell:
                    strain_name = strain_names[i] if i < len(
                        strain_names) else str(i)
                    strains.add(strain_name)
            node_strain_set[nid] = strains

    # 7. Write annotated GML
    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, "pgap2.partition.map_annotated.gml")

    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("graph [\n  directed 0\n")

        for nid in sorted(nodes_raw.keys()):
            cidx = nid  # direct mapping
            if cidx < len(clusters):
                cl = clusters[cidx]
                ca = cluster_attrs_cache[cidx]
                clust_label = cl["clust_id"]        # e.g. "clust_0"
                gene_name = _gml_escape(_clean_bracket_str(cl["gene_name"]))
                product = _gml_escape(_clean_bracket_str(cl["product"]))
                group = cl["group"]
                repre_gene = cl["repre_gene"]
                size = ca["size"]
                involved_strain = cl["involved_strain"]
                para_strain = cl["para_strain"]
                involved_gene = cl["involved_gene"]
                para_gene = cl["para_gene"]
                paralogs = 1 if para_strain > 0 else 0
                genome_ids = _gml_escape(ca["genomeIDs"])
                gene_ids = _gml_escape(ca["geneIDs"])
                members = _gml_escape(ca["members"])
                min_val = cl["min"]
                mean_val = cl["mean"]
                var_val = cl["var"]
                uni_val = cl["uni"]
            else:
                # Unmapped node — keep original label, minimal attributes
                clust_label = nodes_raw[nid]
                gene_name = ""
                product = ""
                group = "Unknown"
                repre_gene = nodes_raw[nid]
                size = 0
                involved_strain = 0
                para_strain = 0
                involved_gene = 0
                para_gene = 0
                paralogs = 0
                genome_ids = ""
                gene_ids = ""
                members = ""
                min_val = "0"
                mean_val = "0"
                var_val = "0"
                uni_val = "0"

            fh.write("  node [\n")
            fh.write(f"    id {nid}\n")
            fh.write(f'    label "{_gml_escape(clust_label)}"\n')
            fh.write(f'    gene_name "{gene_name}"\n')
            fh.write(f'    product "{product}"\n')
            fh.write(f'    group "{group}"\n')
            fh.write(f'    repre_gene "{repre_gene}"\n')
            fh.write(f"    size {size}\n")
            fh.write(f"    involved_strain {involved_strain}\n")
            fh.write(f"    para_strain {para_strain}\n")
            fh.write(f"    involved_gene {involved_gene}\n")
            fh.write(f"    para_gene {para_gene}\n")
            fh.write(f"    paralogs {paralogs}\n")
            fh.write(f"    min {min_val}\n")
            fh.write(f"    mean {mean_val}\n")
            fh.write(f"    var {var_val}\n")
            fh.write(f"    uni {uni_val}\n")
            fh.write(f'    genomeIDs "{genome_ids}"\n')
            fh.write(f'    geneIDs "{gene_ids}"\n')
            fh.write(f'    members "{members}"\n')
            fh.write("  ]\n")

        for s, t in edges_raw:
            # Edge size = number of shared strains between source and target
            s_strains = node_strain_set.get(s, set())
            t_strains = node_strain_set.get(t, set())
            edge_size = len(s_strains & t_strains)

            fh.write("  edge [\n")
            fh.write(f"    source {s}\n")
            fh.write(f"    target {t}\n")
            fh.write(f"    size {edge_size}\n")
            fh.write("  ]\n")

        fh.write("]\n")

    logger.info(f"  Annotated GML written to {out_path}")
    logger.info(
        f"  Summary: {len(nodes_raw)} nodes, {len(edges_raw)} edges, "
        f"{len(clusters)} clusters, {len(strain_names)} strains"
    )
    logger.success("Done")


_GROUP_COLOR = {
    "Strict_core": "#409EFF",
    "Core": "#67C23A",
    "Soft_core": "#E6A23C",
    "Shell": "#F56C6C",
    "Cloud": "#9099A4",
}


def compute_igraph_layout(
    layout_gml_path: str,
    output_json_path: str,
    algorithm: str = "drl",
) -> str:
    """Pre-compute 2-D coordinates for every node using igraph.

    Reads the ``_layout.gml`` file (long-range edges already removed)
    and writes a compact JSON file that the frontend can render with
    D3 Canvas.

    Parameters
    ----------
    layout_gml_path : str
        Path to ``pgap2.partition.map_annotated_layout.gml``.
    output_json_path : str
        Destination JSON file.
    algorithm : str
        igraph layout algorithm.  ``"drl"`` (DrL / distributed recursive
        layout) is recommended for 10 000+ nodes; ``"fr"`` (Fruchterman–
        Reingold) works better for smaller graphs.

    Returns
    -------
    str
        Path to the written JSON file.
    """
    import igraph as ig

    logger.info(f"Computing igraph layout (algorithm={algorithm}) …")

    # ---- 1. Parse GML into igraph ----
    # We parse manually for speed and to preserve per-node attributes.
    node_attrs: Dict[int, dict] = {}       # gml_id → {label, gene_name, …}
    edge_list: List[Tuple[int, int]] = []
    edge_attrs: List[dict] = []

    with open(layout_gml_path, "r", encoding="utf-8") as fh:
        in_node = in_edge = False
        cur: dict = {}
        for raw in fh:
            line = raw.strip()
            if line == "node [":
                in_node, cur = True, {}
            elif line == "edge [":
                in_edge, cur = True, {}
            elif line == "]":
                if in_node:
                    nid = cur.get("id")
                    if nid is not None:
                        node_attrs[nid] = cur
                elif in_edge:
                    s, t = cur.get("source"), cur.get("target")
                    if s is not None and t is not None:
                        edge_list.append((s, t))
                        edge_attrs.append(cur)
                in_node = in_edge = False
            elif in_node or in_edge:
                # Parse key-value
                m = re.match(r'(\w+)\s+(.*)', line)
                if m:
                    key, val = m.group(1), m.group(2)
                    # Strip quotes
                    if val.startswith('"') and val.endswith('"'):
                        val = val[1:-1]
                    else:
                        # Try int / float
                        try:
                            val = int(val)
                        except ValueError:
                            try:
                                val = float(val)
                            except ValueError:
                                pass
                    cur[key] = val

    # Build gml_id → igraph_idx mapping
    sorted_ids = sorted(node_attrs.keys())
    id_to_idx = {gml_id: idx for idx, gml_id in enumerate(sorted_ids)}
    n_nodes = len(sorted_ids)

    logger.info(f"  Nodes: {n_nodes}, Edges: {len(edge_list)}")

    # ---- 2. Build igraph Graph ----
    g = ig.Graph(n=n_nodes, directed=False)
    ig_edges = []
    ig_weights = []
    for (s, t), ea in zip(edge_list, edge_attrs):
        si, ti = id_to_idx.get(s), id_to_idx.get(t)
        if si is not None and ti is not None:
            ig_edges.append((si, ti))
            ig_weights.append(ea.get("size", 1))
    g.add_edges(ig_edges)
    g.es["weight"] = ig_weights

    # ---- 3. Compute layout ----
    if algorithm == "drl":
        lo = g.layout_drl(weights="weight")
    elif algorithm == "fr":
        lo = g.layout_fruchterman_reingold(weights="weight", niter=500)
    elif algorithm == "kk":
        lo = g.layout_kamada_kawai()
    elif algorithm == "graphopt":
        lo = g.layout_graphopt(niter=500)
    else:
        lo = g.layout_drl(weights="weight")

    coords = lo.coords  # list of (x, y) tuples

    # ---- 4. Normalise to [0, width] × [0, height] canvas ----
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    range_x = max_x - min_x if max_x > min_x else 1.0
    range_y = max_y - min_y if max_y > min_y else 1.0

    CANVAS = 4000.0  # virtual canvas size
    PAD = 50.0

    def norm_x(v):
        return round(PAD + (v - min_x) / range_x * (CANVAS - 2 * PAD), 2)

    def norm_y(v):
        return round(PAD + (v - min_y) / range_y * (CANVAS - 2 * PAD), 2)

    # ---- 5. Assemble JSON ----
    json_nodes = []
    for idx, gml_id in enumerate(sorted_ids):
        na = node_attrs[gml_id]
        json_nodes.append({
            "id": gml_id,
            "label": na.get("label", f"clust_{gml_id}"),
            "gene_name": na.get("gene_name", ""),
            "product": na.get("product", ""),
            "group": na.get("group", ""),
            "size": na.get("size", 1),
            "involved_strain": na.get("involved_strain", 0),
            "color": _GROUP_COLOR.get(na.get("group", ""), "#909399"),
            "x": norm_x(coords[idx][0]),
            "y": norm_y(coords[idx][1]),
        })

    json_edges = []
    for (s, t), ea in zip(edge_list, edge_attrs):
        json_edges.append({
            "source": s,
            "target": t,
            "size": ea.get("size", 1),
        })

    # Group statistics
    group_counts: Dict[str, int] = collections.Counter()
    for n in json_nodes:
        group_counts[n["group"]] += 1

    result = {
        "canvas": {"width": CANVAS, "height": CANVAS},
        "total_nodes": n_nodes,
        "total_edges": len(json_edges),
        "groups": dict(group_counts),
        "algorithm": algorithm,
        "nodes": json_nodes,
        "edges": json_edges,
    }

    with open(output_json_path, "w", encoding="utf-8") as fh:
        json.dump(result, fh)

    size_mb = os.path.getsize(output_json_path) / (1024 * 1024)
    logger.info(
        f"  Layout JSON written to {output_json_path} ({size_mb:.1f} MB)")
    return output_json_path


def main(
    indir: str,
    outdir: str,
    ref: Optional[str] = None,
    distance_threshold: int = 100,
    compute_layout_json: bool = False,
    layout_algorithm: str = "drl",
):
    """Unified graph post-processing pipeline.

    1. Generate annotated GML  → ``pgap2.partition.map_annotated.gml``
    2. Reference-based layout  → ``pgap2.partition.map_annotated_layout.gml``
                                 ``pgap2.partition.map_annotated_complete.gml``
    3. (optional) igraph pre-computed layout → ``…_layout.json``
    """
    from pgap2.postprocess.reference_layout import layout, resolve_ref, read_gml

    os.makedirs(outdir, exist_ok=True)

    # Step 1 – annotated GML
    logger.info("=== Step 1: Generate annotated GML ===")
    generate_annotated_gml(indir=indir, outdir=outdir)

    annotated_gml = os.path.join(outdir, "pgap2.partition.map_annotated.gml")
    if not os.path.exists(annotated_gml):
        raise FileNotFoundError(
            f"Annotated GML was not produced: {annotated_gml}"
        )

    # Step 2 – resolve reference strain
    logger.info("=== Step 2: Reference-based layout ===")
    G_tmp = read_gml(annotated_gml)
    ref_name = resolve_ref(G_tmp, ref)
    del G_tmp

    # Step 3 – BFS + min-cut layout
    output_prefix = os.path.join(outdir, "pgap2.partition.map_annotated")
    G_layout, cut_edges = layout(
        gml_path=annotated_gml,
        ref_name=ref_name,
        output_prefix=output_prefix,
        distance_threshold=distance_threshold,
    )

    logger.info(
        f"Pipeline complete. Nodes: {G_layout.number_of_nodes()}, "
        f"Edges kept: {G_layout.number_of_edges()}, "
        f"Edges cut: {len(cut_edges)}"
    )

    # Step 4 – (optional) igraph pre-computed layout → JSON
    if compute_layout_json:
        logger.info("=== Step 4: igraph pre-computed layout ===")
        layout_gml = output_prefix + "_layout.gml"
        layout_json = output_prefix + "_layout.json"
        compute_igraph_layout(
            layout_gml_path=layout_gml,
            output_json_path=layout_json,
            algorithm=layout_algorithm,
        )

    logger.success("Done")


def launch(args: argparse.Namespace):
    outdir = os.path.abspath(
        args.outdir) if args.outdir else os.path.abspath(args.indir)
    os.makedirs(outdir, exist_ok=True)
    indir = os.path.abspath(args.indir)
    assert os.path.exists(indir), f"Input directory does not exist: {indir}"
    main(
        indir=indir,
        outdir=outdir,
        ref=getattr(args, "ref", None),
        distance_threshold=getattr(args, "distance_threshold", 100),
        compute_layout_json=getattr(args, "compute_layout", False),
        layout_algorithm=getattr(args, "layout_algorithm", "drl"),
    )


def postprocess_portal(args):
    outdir = args.outdir if args.outdir else args.indir
    set_verbosity_level(outdir, args.verbose,
                        args.debug, "postprocess_annotate")
    launch(args)


def graph_cmd(subparser: _SubParsersAction):
    p: ArgumentParser = subparser.add_parser(
        "graph",
        help="Generate annotated pangenome graph with reference-based layout",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--indir", "-i", required=True,
        help="Input directory (output of partition step, containing GML, "
             "detail.tsv, and annot.tsv)",
    )
    p.add_argument(
        "--outdir", "-o", required=False, default=None,
        help="Output directory (default: same as input directory)",
    )
    p.add_argument(
        "--ref", "-r", required=False, default=None,
        help="Reference strain name (default: auto-select strain with most nodes)",
    )
    p.add_argument(
        "--distance-threshold", required=False, default=100, type=int,
        help="Gene-order distance threshold for long-range edge detection",
    )
    p.add_argument(
        "--compute-layout", required=False, action="store_true", default=False,
        help="Pre-compute 2-D coordinates with igraph and output a JSON file "
             "for D3 Canvas rendering (requires python-igraph)",
    )
    p.add_argument(
        "--layout-algorithm", required=False, default="drl",
        choices=["drl", "fr", "kk", "graphopt"],
        help="igraph layout algorithm (only used with --compute-layout). "
             "drl=DrL (best for 10k+ nodes), fr=Fruchterman-Reingold, "
             "kk=Kamada-Kawai, graphopt=GraphOpt",
    )
    p.add_argument(
        "--verbose", "-V", required=False, action="store_true", default=False,
        help="Verbose output",
    )
    p.add_argument(
        "--debug", required=False, action="store_true", default=False,
        help="Debug mode",
    )
    p.set_defaults(func=postprocess_portal)
