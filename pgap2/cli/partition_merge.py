"""
CLI definitions for the 'merge' subcommand.
"""
import os
import argparse

from argparse import ArgumentParser, _SubParsersAction

from pgap2.utils.supply import sfw, tqdm_
from pgap2.utils.supply import set_verbosity_level
from pgap2.utils.graph_utils import check_min_falen, check_gcode


def launch(args: argparse.Namespace):
    from pgap2.utils.partition_merge import main
    main(dir_a=os.path.abspath(args.dir_a),
         dir_b=os.path.abspath(args.dir_b),
         outdir=os.path.abspath(args.outdir),
         aligner=args.aligner, clust_method=args.clust_method,
         falen=args.min_falen, threads=args.threads,
         id_attr_key=args.id_attr_key, type_filter=args.type_filter,
         annot=args.annot, gcode=args.gcode, retrieve=args.retrieve,
         disable=args.disable)


def merge_portal(args):
    set_verbosity_level(args.outdir, args.verbose,
                        args.debug, 'partition')

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

    if args.retrieve:
        tqdm_.set_total_step(12)
    else:
        tqdm_.set_total_step(7)
    launch(args)


def merge_cmd(subparser: _SubParsersAction):

    subparser_merge: ArgumentParser = subparser.add_parser(
        'merge', help='Merge two independent PGAP2 results into one pangenome',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparser_merge.add_argument(
        '--dir_a', '-a', required=True,
        help='First PGAP2 result directory (output from pgap2 main or pgap2 add).')
    subparser_merge.add_argument(
        '--dir_b', '-b', required=True,
        help='Second PGAP2 result directory (output from pgap2 main or pgap2 add).')
    subparser_merge.add_argument(
        '--outdir', '-o', required=False, default='.',
        help='Output directory for merged results.')
    subparser_merge.add_argument(
        '--aligner', required=False, type=str,
        default='diamond', choices=('diamond', 'blastp'),
        help='The aligner used for pairwise alignment.')
    subparser_merge.add_argument(
        '--clust_method', required=False, type=str,
        default='mmseqs2', choices=('cdhit', 'mmseqs2'),
        help='The clustering method used in previous runs.')
    subparser_merge.add_argument(
        '--type-filter', required=False, type=str,
        default='CDS',
        help='Only for gff file as input, feature type to filter.')
    subparser_merge.add_argument(
        '--id-attr-key', required=False, type=str,
        default='ID',
        help='Only for gff file as input, attribute key for record ID.')
    subparser_merge.add_argument(
        '--retrieve', '-r', required=False,
        action='store_true', default=False,
        help='Retrieve genes that may be lost with annotations.')
    subparser_merge.add_argument(
        '--threads', '-t', required=False, default=1, type=int,
        help='Threads used in parallel.')
    subparser_merge.add_argument(
        '--min_falen', '-m', required=False, type=check_min_falen,
        default=20,
        help='Protein length of throw_away_sequences, at least 11.')
    subparser_merge.add_argument(
        '--gcode', required=False, type=check_gcode,
        default=11,
        help='The genetic code of your species.')
    subparser_merge.add_argument(
        '--annot', required=False,
        action='store_true', default=False,
        help='Discard original annotation, re-annotate using prodigal.')
    subparser_merge.add_argument(
        '--verbose', '-v', required=False,
        action='store_true', default=False,
        help='Verbose output.')
    subparser_merge.add_argument(
        '--debug', '-D', required=False,
        action='store_true', default=False,
        help='Debug mode. Note: very verbose.')
    subparser_merge.add_argument(
        '--disable', required=False,
        action='store_true', default=False,
        help='Disable progress bar.')

    subparser_merge.set_defaults(func=merge_portal)
