"""
CLI definitions for the 'neigh' (conserved gene neighborhood) subcommand.

This subcommand runs the full partition pipeline (identical to 'main')
followed by a conserved gene neighborhood discovery step inspired by
spacedust's agglomerative clustering algorithm.
"""

import os
import argparse

from argparse import ArgumentParser, _SubParsersAction

from pgap2.utils.supply import sfw, tqdm_
from pgap2.utils.supply import set_verbosity_level
from pgap2.utils.graph_utils import check_min_falen, check_gcode


def launch(args: argparse.Namespace):
    from pgap2.utils.neighborhood import main
    main(
        indir=os.path.abspath(args.indir),
        outdir=os.path.abspath(args.outdir),
        falen=args.min_falen, threads=args.threads, evalue=args.evalue,
        aligner=args.aligner, clust_method=args.clust_method,
        orth_id=args.orth_id, para_id=args.para_id, dup_id=args.dup_id,
        id_attr_key=args.id_attr_key, type_filter=args.type_filter,
        coverage=0.98, fast_mode=args.fast_mode, hconf_thre=args.hconf_thre,
        LD=args.LD, AS=args.AS, AL=args.AL, max_targets=args.max_targets,
        context_similarity=0, flank=5,
        accurate=args.accurate,
        exhaust_orth=args.exhaust_orth,
        gcode=args.gcode,
        disable=args.disable, annot=args.annot, retrieve=args.retrieve,
        radius=args.radius, sensitivity=args.sensitivity, ins=args.ins,
        # neighborhood-specific parameters
        max_gap=args.max_gap,
        min_genes=args.min_genes,
        pval_threshold=args.pval,
        min_genomes=args.min_genomes,
        resolution=args.resolution,
        use_order=args.use_order,
    )


def neighborhood_portal(args):
    set_verbosity_level(args.outdir, args.verbose,
                        args.debug, 'neighborhood')

    if args.clust_method == 'mmseqs2':
        sfw.check_dependency("mmseqs2")
    elif args.clust_method == 'cdhit':
        sfw.check_dependency("cdhit")

    sfw.check_dependency("mcl")
    sfw.check_dependency("mcxdeblast")
    if args.aligner == 'diamond':
        sfw.check_dependency("diamond")
    elif args.aligner == 'blastp':
        sfw.check_dependency("blastp")
        sfw.check_dependency("makeblastdb")

    if args.retrieve:
        sfw.check_dependency("miniprot")
        sfw.check_dependency("seqtk")
        tqdm_.set_total_step(14)  # 12 partition + 2 neighborhood steps
    else:
        tqdm_.set_total_step(9)   # 7 partition + 2 neighborhood steps

    if args.annot:
        sfw.check_dependency("prodigal")

    launch(args)


def neighborhood_cmd(subparser: _SubParsersAction):

    subparser_neigh: ArgumentParser = subparser.add_parser(
        'neigh',
        help='Conserved gene neighborhood (synteny block) analysis. '
             'Runs the full partition pipeline then discovers conserved '
             'gene neighborhoods across genomes (spacedust-style).',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # ---- Partition parameters (same as 'main') ----
    subparser_neigh.add_argument(
        '--indir', '-i', required=True,
        help='Input directory containing genome files. Same prefix = same strain.')
    subparser_neigh.add_argument(
        '--outdir', '-o', required=False, default='.',
        help='Output directory')
    subparser_neigh.add_argument(
        '--dup_id', required=False, type=float, default=0.99,
        help='Maximum identity for most recent duplication events.')
    subparser_neigh.add_argument(
        '--orth_id', required=False, type=float, default=0.98,
        help='Maximum identity between most similar pan-clusters.')
    subparser_neigh.add_argument(
        '--para_id', required=False, type=float, default=0.7,
        help='Paralogous identity threshold.')
    subparser_neigh.add_argument(
        '--type-filter', required=False, type=str, default='CDS',
        help='GFF feature type (3rd column) to include.')
    subparser_neigh.add_argument(
        '--id-attr-key', required=False, type=str, default='ID',
        help='GFF attribute key for record ID.')
    subparser_neigh.add_argument(
        '--hconf_thre', required=False, type=float, default=1,
        help='High confidence cluster threshold.')
    subparser_neigh.add_argument(
        '--exhaust_orth', '-e', required=False, action='store_true',
        default=False, help='Exhaustive paralog splitting.')
    subparser_neigh.add_argument(
        '--sensitivity', '-s', required=False, type=str, default='strict',
        choices=('soft', 'moderate', 'strict'),
        help='Node connectivity sensitivity.')
    subparser_neigh.add_argument(
        '--ins', '-n', required=False, action='store_true', default=False,
        help='Ignore insertion sequence influence.')
    subparser_neigh.add_argument(
        '--fast', '-f', dest='fast_mode', required=False,
        action='store_true', default=False,
        help='Fast mode (skip fine feature analysis).')
    subparser_neigh.add_argument(
        '--accurate', '-a', required=False, action='store_true',
        default=False, help='Bidirection check for paralog partition.')
    subparser_neigh.add_argument(
        '--retrieve', '-r', required=False, action='store_true',
        default=False, help='Retrieve genes lost during annotation.')
    subparser_neigh.add_argument(
        '--threads', '-t', required=False, default=1, type=int,
        help='Number of threads.')
    subparser_neigh.add_argument(
        '--max_targets', '-k', required=False, type=int, default=2000,
        help='Maximum alignment targets per query.')
    subparser_neigh.add_argument(
        '--LD', required=False, type=float, default=0.6,
        help='Minimum gene length difference proportion.')
    subparser_neigh.add_argument(
        '--AS', required=False, type=float, default=0.6,
        help='Coverage for shorter sequence.')
    subparser_neigh.add_argument(
        '--AL', required=False, type=float, default=0.6,
        help='Coverage for longer sequence.')
    subparser_neigh.add_argument(
        '--evalue', required=False, type=float, default=1E-5,
        help='Aligner E-value threshold.')
    subparser_neigh.add_argument(
        '--aligner', required=False, type=str, default='diamond',
        choices=('diamond', 'blastp'), help='Pairwise aligner.')
    subparser_neigh.add_argument(
        '--clust_method', required=False, type=str, default='mmseqs2',
        choices=('cdhit', 'mmseqs2'), help='Gene clustering method.')
    subparser_neigh.add_argument(
        '--radius', required=False, type=int, default=3,
        help='Search region radius.')
    subparser_neigh.add_argument(
        '--min_falen', '-m', required=False, type=check_min_falen,
        default=20, help='Minimum protein length (at least 11).')
    subparser_neigh.add_argument(
        '--gcode', required=False, type=check_gcode, default=11,
        help='Genetic code.')
    subparser_neigh.add_argument(
        '--annot', required=False, action='store_true', default=False,
        help='Re-annotate genomes with Prodigal.')

    # ---- Neighborhood-specific parameters ----
    neigh_group = subparser_neigh.add_argument_group(
        'Neighborhood parameters',
        'Parameters specific to conserved gene neighborhood detection.')
    neigh_group.add_argument(
        '--max_gap', required=False, type=int, default=3,
        help='Maximum gene gap allowed between hits for merging into a '
             'syntenic block. Larger values allow more intervening genes.')
    neigh_group.add_argument(
        '-mg','--min_genes', required=False, type=int, default=2,
        help='Minimum number of genes to define a syntenic block.')
    neigh_group.add_argument(
        '--pval', required=False, type=float, default=0.01,
        help='P-value threshold for gene-order filtering (only used '
             'when --use_order is enabled).')
    neigh_group.add_argument(
        '--use_order', required=False, action='store_true', default=False,
        help='Enable gene-order conservation filtering. When set, blocks '
             'are filtered by P_order <= --pval, edge weights use '
             '-log10(P_order), and neighborhood P-values are computed '
             'via Fisher\'s combined test. Default: off (accept blocks '
             'purely based on --min_genes co-localization).')
    neigh_group.add_argument(
        '-mgs','--min_genomes', required=False, type=int, default=2,
        help='Minimum number of genomes a neighborhood must span to be reported.')
    neigh_group.add_argument(
        '-re','--resolution', required=False, type=float, default=1.0,
        help='Louvain resolution parameter. Higher values produce smaller, '
             'more focused neighborhoods; lower values produce larger, '
             'more inclusive ones. Try 2.0-10.0 if neighborhoods are too large.')

    # ---- Common flags ----
    subparser_neigh.add_argument(
        '--verbose', '-v', required=False, action='store_true',
        default=False, help='Verbose output')
    subparser_neigh.add_argument(
        '--debug', '-D', required=False, action='store_true',
        default=False, help='Debug mode (very verbose)')
    subparser_neigh.add_argument(
        '--disable', required=False, action='store_true',
        default=False, help='Disable progress bar')

    subparser_neigh.set_defaults(func=neighborhood_portal)
