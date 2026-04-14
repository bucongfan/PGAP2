"""
Orchestrator for the 'neigh' (neighborhood) subcommand.

This module reuses the entire partition pipeline from utils/partition.py,
then appends a conserved gene neighborhood discovery step.

The workflow is:
  1. Identical to `pgap2 main`: load data → cluster → align → MCL →
     build network → merge by similarity → merge by synteny → organize results.
  2. After partitioning, run the synteny block detection algorithm
     (lib/synteny.py) on the resulting pangenome graph.
  3. Output neighborhood results (TSV + GML).
"""

import os
import pickle
import networkx as nx

from loguru import logger

from pgap2.lib.basic import Basic
from pgap2.lib.pklcheck import PklCheck
from pgap2.lib.pangenome import Pangenome
from pgap2.lib.tree import Tree

from pgap2.lib.synteny import (
    discover_neighborhoods,
    write_neighborhood_tsv,
    write_neighborhood_detail_tsv,
    write_neighborhood_gml,
)

from pgap2.utils.supply import tqdm_


def main(
    indir: str, outdir: str, evalue: float, hconf_thre: float,
    aligner: str, clust_method: str, falen: int, fast_mode: bool,
    threads: int, orth_id: float, para_id: float, dup_id: float,
    id_attr_key: str, type_filter: str, max_targets: int,
    coverage: float, LD: float, AS: float, AL: float,
    context_similarity: float, accurate: bool, exhaust_orth: bool,
    flank: int, disable: bool, annot: bool, gcode: int,
    retrieve: bool, radius: int, sensitivity: int, ins: bool,
    # neighborhood-specific parameters
    max_gap: int = 3,
    min_genes: int = 2,
    pval_threshold: float = 0.01,
    min_genomes: int = 2,
    resolution: float = 1.0,
    use_order: bool = True,
):
    """
    Run the full partition pipeline, then discover conserved gene neighborhoods.
    """

    from pgap2.utils.partition import (
        main as partition_main,
    )

    # We call partition_main to produce all standard outputs.
    # After it returns, we reload the pickle files it saved to get
    # the G, pg, tree objects for neighborhood analysis.
    partition_main(
        indir=indir, outdir=outdir, evalue=evalue,
        hconf_thre=hconf_thre, aligner=aligner,
        clust_method=clust_method, falen=falen,
        fast_mode=fast_mode, threads=threads,
        orth_id=orth_id, para_id=para_id, dup_id=dup_id,
        id_attr_key=id_attr_key, type_filter=type_filter,
        max_targets=max_targets, coverage=coverage,
        LD=LD, AS=AS, AL=AL,
        context_similarity=context_similarity,
        accurate=accurate, exhaust_orth=exhaust_orth,
        flank=flank, disable=disable, annot=annot,
        gcode=gcode, retrieve=retrieve,
        radius=radius, sensitivity=sensitivity, ins=ins,
    )


    logger.info("=" * 60)
    logger.info("Phase 2: Conserved Gene Neighborhood Discovery")
    logger.info("=" * 60)

    # Load the graph object saved by partition
    graph_pkl = os.path.join(outdir, "graph.pkl")
    if not os.path.exists(graph_pkl):
        logger.error(f"Cannot find {graph_pkl}. Partition may have failed.")
        raise FileNotFoundError(graph_pkl)

    with open(graph_pkl, 'rb') as fh:
        graph_check: PklCheck = pickle.load(fh)
        G: nx.Graph = graph_check.data_dump('graph')

    # Load the basic object (contains pg reference)
    basic_pkl = os.path.join(outdir, "basic.pkl")
    if not os.path.exists(basic_pkl):
        logger.error(f"Cannot find {basic_pkl}. Partition may have failed.")
        raise FileNotFoundError(basic_pkl)

    with open(basic_pkl, 'rb') as fh:
        basic_check: PklCheck = pickle.load(fh)
        basic: Basic = basic_check.data_dump('basic')

    # Load the preprocess pickle (contains tree)
    prep_pkl = os.path.join(outdir, "preprocess.pkl")
    if not os.path.exists(prep_pkl):
        logger.error(f"Cannot find {prep_pkl}. Partition may have failed.")
        raise FileNotFoundError(prep_pkl)

    with open(prep_pkl, 'rb') as fh:
        prep_check: PklCheck = pickle.load(fh)
        pg: Pangenome = prep_check.data_dump('pangenome')
        tree: Tree = prep_check.data_dump('tree')

    # Reload annotation so pg.annot is available
    pg.reload_annot_file(retrieve=retrieve)

    neighborhoods = discover_neighborhoods(
        G=G, pg=pg, tree=tree,
        max_gap=max_gap,
        min_genes=min_genes,
        pval_threshold=pval_threshold,
        min_genomes=min_genomes,
        resolution=resolution,
        use_order=use_order,
        disable=disable,
        threads=threads,
    )

    prefix = "pgap2.neighborhood"
    write_neighborhood_tsv(neighborhoods, pg, G, outdir, prefix)
    write_neighborhood_detail_tsv(neighborhoods, pg, G, outdir, prefix)
    write_neighborhood_gml(neighborhoods, outdir, prefix)

    logger.success(
        f"Conserved gene neighborhood analysis complete. "
        f"Found {len(neighborhoods)} neighborhoods in {outdir}/")

    return 0
