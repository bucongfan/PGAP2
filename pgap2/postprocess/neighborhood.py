import os
import pickle
import argparse

import networkx as nx
from loguru import logger
from argparse import ArgumentParser, _SubParsersAction

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

from pgap2.utils.supply import set_verbosity_level


def main(
    indir: str,
    outdir: str,
    threads: int = 1,
    disable: bool = False,
    retrieve: bool = False,
    max_gap: int = 3,
    min_genes: int = 2,
    pval_threshold: float = 0.01,
    min_genomes: int = 2,
    resolution: float = 1.0,
    use_order: bool = False,
):
    """Discover conserved gene neighborhoods from a completed partition."""

    logger.info("=" * 60)
    logger.info("Conserved Gene Neighborhood Discovery (post-process)")
    logger.info("=" * 60)

    # ── Load pickle files produced by partition ──
    for name in ("graph.pkl", "basic.pkl", "preprocess.pkl"):
        path = os.path.join(indir, name)
        if not os.path.exists(path):
            logger.error(f"Cannot find {path}. "
                         "The input dir should be the output dir of the partition step.")
            raise FileNotFoundError(path)

    annot_file = os.path.join(indir, "total.involved_annot.tsv")
    if not os.path.exists(annot_file):
        logger.error(f"Cannot find {annot_file}.")
        raise FileNotFoundError(annot_file)

    with open(os.path.join(indir, "graph.pkl"), 'rb') as fh:
        graph_check: PklCheck = pickle.load(fh)
        G: nx.Graph = graph_check.data_dump('graph')

    with open(os.path.join(indir, "preprocess.pkl"), 'rb') as fh:
        prep_check: PklCheck = pickle.load(fh)
        pg: Pangenome = prep_check.data_dump('pangenome')
        tree: Tree = prep_check.data_dump('tree')

    # Reload annotation so pg.annot is available
    pg.reload_annot_file(retrieve=retrieve)

    # ── Discover neighborhoods ──
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

    # ── Write outputs ──
    os.makedirs(outdir, exist_ok=True)
    prefix = "pgap2.neighborhood"
    write_neighborhood_tsv(neighborhoods, pg, G, outdir, prefix)
    write_neighborhood_detail_tsv(neighborhoods, pg, G, outdir, prefix)
    write_neighborhood_gml(neighborhoods, outdir, prefix)

    logger.success(
        f"Conserved gene neighborhood analysis complete. "
        f"Found {len(neighborhoods)} neighborhoods in {outdir}/")

    return neighborhoods


def launch(args: argparse.Namespace):
    outdir = os.path.abspath(
        args.outdir) if args.outdir else os.path.abspath(args.indir)
    os.makedirs(outdir, exist_ok=True)
    indir = os.path.abspath(args.indir)
    main(
        indir=indir,
        outdir=outdir,
        threads=args.threads,
        disable=args.disable,
        retrieve=getattr(args, 'retrieve', False),
        max_gap=args.max_gap,
        min_genes=args.min_genes,
        pval_threshold=args.pval,
        min_genomes=args.min_genomes,
        resolution=args.resolution,
        use_order=args.use_order,
    )


def postprocess_portal(args):
    outdir = args.outdir if args.outdir else args.indir
    set_verbosity_level(outdir, args.verbose,
                        args.debug, 'postprocess_neighborhood')
    launch(args)


def neigh_cmd(subparser: _SubParsersAction):
    p: ArgumentParser = subparser.add_parser(
        'neigh',
        help='Conserved gene neighborhood (synteny block) analysis '
             'on an existing partition output.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        '--indir', '-i', required=True,
        help='Input directory (output of partition step, containing '
             'graph.pkl, preprocess.pkl, and total.involved_annot.tsv)',
    )
    p.add_argument(
        '--outdir', '-o', required=False, default=None,
        help='Output directory (default: same as input directory)',
    )
    p.add_argument(
        '--threads', '-t', required=False, default=1, type=int,
        help='Number of threads',
    )
    p.add_argument(
        '--retrieve', '-r', required=False, action='store_true',
        default=False,
        help='Retrieve genes lost during annotation (passed to reload_annot_file)',
    )
    p.add_argument(
        '--max_gap', required=False, type=int, default=3,
        help='Maximum gene gap allowed between hits for merging into a '
             'syntenic block.',
    )
    p.add_argument(
        '-mg', '--min_genes', required=False, type=int, default=2,
        help='Minimum number of genes to define a syntenic block.',
    )
    p.add_argument(
        '--pval', required=False, type=float, default=0.01,
        help='P-value threshold for gene-order filtering.',
    )
    p.add_argument(
        '--use_order', required=False, action='store_true', default=False,
        help='Enable gene-order conservation filtering.',
    )
    p.add_argument(
        '-mgs', '--min_genomes', required=False, type=int, default=2,
        help='Minimum number of genomes a neighborhood must span.',
    )
    p.add_argument(
        '-re', '--resolution', required=False, type=float, default=1.0,
        help='Louvain resolution parameter.',
    )
    p.add_argument(
        '--verbose', '-V', required=False, action='store_true',
        default=False, help='Verbose output',
    )
    p.add_argument(
        '--debug', required=False, action='store_true',
        default=False, help='Debug mode',
    )
    p.add_argument(
        '--disable', required=False, action='store_true',
        default=False, help='Disable progress bar',
    )
    p.set_defaults(func=postprocess_portal)
