"""
BFS + min-cut reference layout for PGAP2 annotated pangenome graph.

Reads ``pgap2.partition.map_annotated.gml`` (produced by the annotation step
in :mod:`pgap2.postprocess.graph`) and severs long-range connections:

  1. Build a reference-gene mapping from ``genomeIDs`` / ``geneIDs``.
  2. Optionally add edges between consecutive reference genes.
  3. BFS from each reference node to discover long-range sinks.
  4. Min-cut on the induced subgraph to sever long-range connections.
  5. Write processed graph and cut-edge information.

Input
-----
- ``pgap2.partition.map_annotated.gml``
  Nodes carry ``genomeIDs`` (semicolon-separated strain names),
  ``geneIDs`` (semicolon-separated ``strain_idx:contig:gene_order``),
  ``gene_name``, ``product``, ``group``, ``size``, etc.
  Edges carry ``size`` (# shared strains).

Output
------
- ``pgap2.partition.map_annotated_layout.gml``   – long-range edges removed
- ``pgap2.partition.map_annotated_complete.gml``  – all edges, ``is_cut`` attr
"""

import collections
from typing import Dict, List, Set, Tuple

import networkx as nx
from loguru import logger


# ===========================================================================
# GML Reader
# ===========================================================================

def read_gml(gml_path: str) -> nx.Graph:
    """Read PGAP2 annotated GML.

    Attaches per-node:
    - ``strain_set``  : set of strain *name* strings (e.g. {"EC01","EC02"})
    - ``orig_label``  : the ``label`` attribute (e.g. ``"clust_0"``)
    """
    G = nx.read_gml(gml_path, label="id")

    for n in G.nodes():
        nd = G.nodes[n]
        nd["orig_label"] = nd.get("label", str(n))

        genome_ids_str = nd.get("genomeIDs", "")
        nd["strain_set"] = (
            set(genome_ids_str.split(";")) if genome_ids_str else set()
        )

    return G


# ===========================================================================
# Reference-gene mapping
# ===========================================================================

def create_mapping(G: nx.Graph, ref_name: str) -> Dict[int, int]:
    """Build mapping ``{node_id: gene_order}`` for the reference strain.

    ``ref_name`` is the strain name (e.g. ``"EC01"``).
    ``geneIDs`` / ``genomeIDs`` are parallel semicolon-separated lists;
    gene indices have format ``strain_idx:contig:gene_order``.
    """
    # Build strain_name → strain_idx lookup from genomeIDs/geneIDs
    name_to_idx: Dict[str, Set[str]] = {}
    for n in G.nodes():
        nd = G.nodes[n]
        genome_ids = nd.get("genomeIDs", "").split(";")
        gene_ids = nd.get("geneIDs", "").split(";")
        for gname, gi in zip(genome_ids, gene_ids):
            parts = gi.split(":")
            if len(parts) == 3:
                name_to_idx.setdefault(gname.strip(), set()).add(parts[0])

    ref_strain_idxs = name_to_idx.get(ref_name, set())
    if not ref_strain_idxs:
        raise ValueError(
            f"Strain name '{ref_name}' not found in genomeIDs. "
            f"Available: {sorted(name_to_idx.keys())}"
        )

    mapping: Dict[int, int] = {}
    for n in G.nodes():
        nd = G.nodes[n]
        if ref_name not in nd.get("strain_set", set()):
            continue
        gene_ids_str = nd.get("geneIDs", "")
        if not gene_ids_str:
            continue
        for gi in gene_ids_str.split(";"):
            parts = gi.split(":")
            if len(parts) == 3 and parts[0] in ref_strain_idxs:
                mapping[n] = int(parts[2])
                break
    return mapping


# ===========================================================================
# Resolve reference strain
# ===========================================================================

def resolve_ref(G: nx.Graph, ref_name: str = None) -> str:
    """Return the reference strain name.

    If *ref_name* is ``None``, auto-select the strain appearing in the
    most nodes.
    """
    if ref_name is not None:
        # Validate
        for n in G.nodes():
            if ref_name in G.nodes[n].get("strain_set", set()):
                return ref_name
        all_strains: Set[str] = set()
        for n in G.nodes():
            all_strains.update(G.nodes[n].get("strain_set", set()))
        raise ValueError(
            f"Strain '{ref_name}' not found. "
            f"Available: {sorted(all_strains)}"
        )

    # Auto-select
    genome_count: Dict[str, int] = collections.Counter()
    for n in G.nodes():
        gids = G.nodes[n].get("genomeIDs", "")
        if gids:
            for gid in gids.split(";"):
                g = gid.strip()
                if g:
                    genome_count[g] += 1

    best = genome_count.most_common(1)[0][0]
    logger.info(
        f"Auto-selected reference: strain='{best}' "
        f"({genome_count[best]} nodes)")
    return best


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _get_dist(pos_s: int, pos_t: int, max_dist: int) -> int:
    """Circular distance between two gene positions."""
    return min(abs(pos_s - pos_t), abs(abs(pos_s - pos_t) - max_dist))


def _add_to_queue(
    G: nx.Graph,
    source_node: int,
    neighbor_nodes: List[int],
    visited: Set[int],
    sink: Dict,
    mapping: Dict[int, int],
    ref_id: str,
    max_dist: int,
    distance_threshold: int = 100,
) -> List[int]:
    """Explore neighbours of *source_node* during BFS."""
    add: List[int] = []
    for i in neighbor_nodes:
        if i in visited:
            continue

        if ref_id not in G.nodes[i]["strain_set"]:
            add.append(i)
        else:
            if i not in mapping:
                sink["sink"] = i
                visited.add(i)
                continue
            dist = _get_dist(mapping[source_node], mapping[i], max_dist)
            if dist > distance_threshold:
                sink["sink"] = i
        visited.add(i)
    return add


# ===========================================================================
# GML writing helpers
# ===========================================================================

def _gml_attr(fh, key, val, indent="    "):
    """Write one GML key-value pair."""
    if isinstance(val, list):
        for item in val:
            _gml_attr(fh, key, item, indent)
    elif isinstance(val, str):
        val_esc = val.replace("\\", "\\\\").replace('"', '\\"')
        fh.write(f'{indent}{key} "{val_esc}"\n')
    elif isinstance(val, bool):
        fh.write(f"{indent}{key} {1 if val else 0}\n")
    elif isinstance(val, int):
        fh.write(f"{indent}{key} {val}\n")
    elif isinstance(val, float):
        fh.write(f"{indent}{key} {val}\n")


def _write_gml(fh, graph, extra_edge_fn=None):
    """Write *graph* as GML preserving all primitive attributes."""
    fh.write("graph [\n  directed 0\n")
    for n in graph.nodes():
        fh.write("  node [\n")
        fh.write(f"    id {n}\n")
        for k, v in graph.nodes[n].items():
            _gml_attr(fh, k, v)
        fh.write("  ]\n")
    for s, t, data in graph.edges(data=True):
        fh.write("  edge [\n")
        fh.write(f"    source {s}\n")
        fh.write(f"    target {t}\n")
        for k, v in data.items():
            if k == "capacity":
                continue
            _gml_attr(fh, k, v)
        if extra_edge_fn:
            extra_edge_fn(fh, s, t)
        fh.write("  ]\n")
    fh.write("]\n")


# ===========================================================================
# Main layout algorithm
# ===========================================================================

def layout(
    gml_path: str,
    ref_name: str,
    output_prefix: str,
    add_reference_edges: bool = False,
    distance_threshold: int = 100,
) -> Tuple[nx.Graph, List[Tuple]]:
    """Detect and sever long-range connections in the annotated pangenome graph.

    Parameters
    ----------
    gml_path : str
        Path to ``pgap2.partition.map_annotated.gml``.
    ref_name : str
        Strain name used as reference (e.g. ``"EC01"``).
    output_prefix : str
        Prefix for output files (without ``.gml``).
    add_reference_edges : bool
        Add edges between consecutive reference genes.
    distance_threshold : int
        Gene-order distance beyond which a ref node is a long-range sink.

    Returns
    -------
    G : nx.Graph
        The graph with long-range edges removed.
    cut_edges : list[tuple]
        List of ``(source, target)`` edges that were cut.
    """
    logger.info(f"Reading annotated graph: {gml_path}")
    G = read_gml(gml_path)
    logger.info(
        f"  Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

    # Build reference mapping
    mapping = create_mapping(G, ref_name)
    if not mapping:
        all_strains: Set[str] = set()
        for n in G.nodes():
            all_strains.update(G.nodes[n]["strain_set"])
        raise ValueError(
            f"No reference nodes found for strain='{ref_name}'. "
            f"Available: {sorted(all_strains)}"
        )
    logger.info(f"  Reference strain: {ref_name}, nodes: {len(mapping)}")

    gene_orders = list(mapping.values())
    max_dist = max(gene_orders)
    logger.info(f"  Max gene order (genome length): {max_dist}")

    # Optionally add reference edges
    if add_reference_edges:
        sorted_ref = sorted(mapping.items(), key=lambda x: x[1])
        added = 0
        for j in range(1, len(sorted_ref)):
            n1, n2 = sorted_ref[j - 1][0], sorted_ref[j][0]
            if not G.has_edge(n1, n2):
                G.add_edge(n1, n2)
                added += 1
        logger.info(f"  Added {added} reference edges")

    # Set edge capacities for min-cut (use edge 'size' if available)
    for e in G.edges:
        G.edges[e]["capacity"] = G.edges[e].get("size", 1)

    # ----------------------------------------------------------------
    # BFS + min-cut loop
    # ----------------------------------------------------------------
    cut_edges: List[Tuple] = []
    sorted_ref = sorted(mapping.items(), key=lambda x: x[1])
    ref_node_list = [n for n, _ in sorted_ref]

    i = 0
    while i < len(ref_node_list):
        nid = ref_node_list[i]

        if nid not in G or nid not in mapping:
            i += 1
            continue

        visited = set([nid])
        sink: Dict = {"sink": None}

        queue = _add_to_queue(
            G, nid, list(G.neighbors(nid)), visited, sink,
            mapping, ref_name, max_dist, distance_threshold,
        )

        while queue:
            target = queue.pop(0)
            visited.add(target)
            queue += _add_to_queue(
                G, nid, list(G.neighbors(target)), visited, sink,
                mapping, ref_name, max_dist, distance_threshold,
            )

        if sink["sink"] is not None:
            logger.debug(
                f"  Node {nid} (order={mapping[nid]}): "
                f"long-range sink {sink['sink']}")
            visited.add(sink["sink"])

            s_t_graph = G.subgraph(visited).copy()

            # Remove nearby ref-ref edges (should not be cut)
            remove = []
            for e in s_t_graph.edges:
                n0, n1 = e[0], e[1]
                if (ref_name in G.nodes[n0]["strain_set"]
                        and ref_name in G.nodes[n1]["strain_set"]):
                    if n0 in mapping and n1 in mapping:
                        if abs(mapping[n0] - mapping[n1]) < distance_threshold:
                            remove.append(e)
            s_t_graph.remove_edges_from(remove)

            try:
                _, partitions = nx.algorithms.flow.minimum_cut(
                    s_t_graph, nid, sink["sink"]
                )
            except nx.NetworkXError as exc:
                logger.warning(f"  Min-cut failed at node {nid}: {exc}")
                i += 1
                continue

            cut = [
                (p1, p2)
                for p1 in partitions[0]
                for p2 in partitions[1]
                if s_t_graph.has_edge(p1, p2)
            ]

            if not cut:
                logger.warning(f"  No min-cut edges found at node {nid}")
                i += 1
                continue

            for e in cut:
                lbl0 = G.nodes[e[0]].get("label", str(e[0]))
                lbl1 = G.nodes[e[1]].get("label", str(e[1]))
                logger.debug(f"    Cut: {lbl0} -- {lbl1}")
                cut_edges.append(e)

            G.remove_edges_from(cut)
            sink["sink"] = None
        else:
            i += 1
            sink["sink"] = None

    logger.info(f"Total cut edges: {len(cut_edges)}")

    # ----------------------------------------------------------------
    # Prepare output
    # ----------------------------------------------------------------
    cut_edge_set = set()
    for e in cut_edges:
        cut_edge_set.add((e[0], e[1]))
        cut_edge_set.add((e[1], e[0]))

    # Clean up internal-only attributes
    for n in G.nodes():
        G.nodes[n]["label"] = G.nodes[n].get("orig_label", str(n))
        for attr in ("strain_set", "orig_label"):
            if attr in G.nodes[n]:
                del G.nodes[n][attr]

    # ---- Output 1: _layout.gml (cut edges REMOVED) ----
    out_layout = output_prefix + "_layout.gml"
    with open(out_layout, "w") as fh:
        _write_gml(fh, G)
    logger.info(f"Wrote layout GML (cut edges removed): {out_layout}")

    # ---- Output 2: _complete.gml (ALL edges, is_cut attribute) ----
    G_orig = nx.read_gml(gml_path, label="id")
    out_complete = output_prefix + "_complete.gml"
    with open(out_complete, "w") as fh:
        def _add_is_cut(fh, s, t):
            is_cut = 1 if (s, t) in cut_edge_set else 0
            fh.write(f"    is_cut {is_cut}\n")
        _write_gml(fh, G_orig, extra_edge_fn=_add_is_cut)
    logger.info(f"Wrote complete GML (with is_cut attr): {out_complete}")

    return G, cut_edges
