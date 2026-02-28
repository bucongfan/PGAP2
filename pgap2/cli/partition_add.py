"""
CLI definitions for the 'add' subcommand.
"""
import os
import argparse

from argparse import ArgumentParser, _SubParsersAction

from pgap2.utils.supply import sfw, tqdm_
from pgap2.utils.supply import set_verbosity_level
from pgap2.utils.graph_utils import check_min_falen, check_gcode


def launch(args: argparse.Namespace):
    from pgap2.utils.partition_add import main
    main(indir=os.path.abspath(args.indir), outdir=os.path.abspath(args.outdir),
         previous_dir=os.path.abspath(args.previous),
         aligner=args.aligner, clust_method=args.clust_method,
         falen=args.min_falen, threads=args.threads,
         id_attr_key=args.id_attr_key, type_filter=args.type_filter,
         annot=args.annot, gcode=args.gcode, retrieve=args.retrieve,
         disable=args.disable)


def add_portal(args):
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
    if args.annot:
        sfw.check_dependency("prodigal")
    if args.retrieve:
        tqdm_.set_total_step(12)
    else:
        tqdm_.set_total_step(7)
    launch(args)


def add_cmd(subparser: _SubParsersAction):

    subparser_add: ArgumentParser = subparser.add_parser(
        'add', help='Add one or more genomes to an existing PGAP2 project', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparser_add.add_argument(
        '--indir', '-i', required=True, help='Input file contained, same prefix seems as the same strain.',)
    subparser_add.add_argument(
        '--outdir', '-o', required=False, help='Output directory', default='.',)
    subparser_add.add_argument(
        '--previous', '-p', required=False, help='The previous PGAP2 result directory, used to resume the partition step quickly.', default=None,)
    subparser_add.add_argument('--aligner', required=False, type=str,
                               default='diamond', choices=('diamond', 'blastp'), help='The aligner used to pairwise alignment.')
    subparser_add.add_argument('--clust_method', required=False, type=str,
                               default='mmseqs2', choices=('cdhit', 'mmseqs2'), help='The method used to cluster the genes.')

    subparser_add.add_argument("--type-filter", required=False, type=str,
                               default='CDS', help="Only for gff file as input, feature type (3rd column) to include, Only lines matching these types will be processed.")
    subparser_add.add_argument("--id-attr-key", required=False, type=str,
                               default='ID', help="Only for gff file as input, Attribute key to extract from the 9th column as the record ID (e.g., 'ID', 'gene', 'locus_tag').")
    subparser_add.add_argument('--retrieve', '-r', required=False,
                               action='store_true', default=False, help='Retrieve gene that may lost with annotations')
    subparser_add.add_argument(
        '--threads', '-t', required=False, default=1, help='threads used in parallel', type=int)
    subparser_add.add_argument('--min_falen', '-m', required=False, type=check_min_falen,
                               default=20, help='protein length of throw_away_sequences, at least 11')
    subparser_add.add_argument('--gcode', required=False, type=check_gcode,
                               default=11, help='The genetic code of your species.')
    subparser_add.add_argument('--annot', required=False,
                               action='store_true', default=False, help='Discard original annotation, and re-annote the genome privately using prodigal')

    subparser_add.add_argument(
        '--verbose', '-v', required=False, action='store_true', default=False, help='Verbose output')
    subparser_add.add_argument(
        '--debug', '-D', required=False, action='store_true', default=False, help='Debug mode. Note: very verbose')
    subparser_add.add_argument(
        '--disable', required=False, action='store_true', default=False, help='Disable progress bar')

    subparser_add.set_defaults(func=add_portal)
