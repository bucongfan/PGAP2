"""
CLI definitions for the 'main' (partition) subcommand.
"""
import os
import argparse

from argparse import ArgumentParser, _SubParsersAction

from pgap2.utils.supply import sfw, tqdm_
from pgap2.utils.supply import set_verbosity_level
from pgap2.utils.graph_utils import check_min_falen, check_gcode


def launch(args: argparse.Namespace):
    from pgap2.utils.partition import main
    main(indir=os.path.abspath(args.indir), outdir=os.path.abspath(args.outdir),
         falen=args.min_falen, threads=args.threads, evalue=args.evalue,
         aligner=args.aligner, clust_method=args.clust_method,
         orth_id=args.orth_id, para_id=args.para_id, dup_id=args.dup_id,
         id_attr_key=args.id_attr_key, type_filter=args.type_filter,
         coverage=0.98, fast_mode=args.fast_mode, hconf_thre=args.hconf_thre,
         #  coverage=args.coverage,
         LD=args.LD, AS=args.AS, AL=args.AL, max_targets=args.max_targets,
         # notused parameters, set a default value and will discard in the next release if dont use for sure
         context_similarity=0, flank=5,
         accurate=args.accurate,
         exhaust_orth=args.exhaust_orth,
         gcode=args.gcode,
         disable=args.disable, annot=args.annot, retrieve=args.retrieve,
         radius=args.radius, sensitivity=args.sensitivity, ins=args.ins)


def partition_portal(args):
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
        tqdm_.set_total_step(12)
    else:
        tqdm_.set_total_step(7)

    if args.annot:
        sfw.check_dependency("prodigal")

    launch(args)


def partition_cmd(subparser: _SubParsersAction):

    subparser_partition: ArgumentParser = subparser.add_parser(
        'main', help='Core functions of PGAP2', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparser_partition.add_argument(
        '--indir', '-i', required=True, help='Input file contained, same prefix seems as the same strain.',)
    subparser_partition.add_argument(
        '--outdir', '-o', required=False, help='Output directory', default='.',)
    subparser_partition.add_argument('--dup_id', required=False, type=float,
                                     default=0.99, help='The maximum identity between the most recent duplication envent.')
    subparser_partition.add_argument('--orth_id', required=False, type=float,
                                     default=0.98, help='The maximum identity between the most similar panclusters.')
    subparser_partition.add_argument('--para_id', required=False, type=float,
                                     default=0.7, help='Use this identity as the paralogous identity.')
    subparser_partition.add_argument("--type-filter", required=False, type=str,
                                     default='CDS', help="Only for gff file as input, feature type (3rd column) to include, Only lines matching these types will be processed.")
    subparser_partition.add_argument("--id-attr-key", required=False, type=str,
                                     default='ID', help="Only for gff file as input, Attribute key to extract from the 9th column as the record ID (e.g., 'ID', 'gene', 'locus_tag').")
    subparser_partition.add_argument('--hconf_thre', required=False, type=float,
                                     default=1, help='The threshold to define high confidence cluster which is used to evaluate the cluster diversity. Loose this value when your input is too large or too diverse, such as 0.95.')
    subparser_partition.add_argument('--exhaust_orth', '-e', required=False, action='store_true',
                                     default=False, help='Try to split every paralogs gene exhausted')
    subparser_partition.add_argument('--sensitivity', '-s', required=False, type=str,
                                     default='strict', choices=('soft', 'moderate', 'strict'), help='The degree of connectedness between each node of clust.')
    subparser_partition.add_argument('--ins', '-n', required=False,
                                     action='store_true', default=False, help='Ignore the influence of insertion sequence.')
    subparser_partition.add_argument('--fast', '-f', dest='fast_mode', required=False,
                                     action='store_true', default=False, help='Do not apply fine feature analysis just partition according to the gene identity and synteny.')
    subparser_partition.add_argument('--accurate', '-a', required=False,
                                     action='store_true', default=False, help='Apply bidirection check for paralogous gene partition (useless if exhaust_orth asigned).')
    subparser_partition.add_argument('--retrieve', '-r', required=False,
                                     action='store_true', default=False, help='Retrieve gene that may lost with annotations')
    subparser_partition.add_argument(
        '--threads', '-t', required=False, default=1, help='threads used in parallel', type=int)
    subparser_partition.add_argument('--max_targets', '-k', required=False, type=int,
                                     default=2000, help='The maximum targets for each query in alignment. Improves accuracy for large-scale analyses, but increases runtime and memory usage.')
    subparser_partition.add_argument('--LD', required=False, type=float,
                                     default=0.6, help='Minimum gene length difference proportion between two genes.')
    subparser_partition.add_argument('--AS', required=False, type=float,
                                     default=0.6, help='Coverage for the shorter sequence.')
    subparser_partition.add_argument('--AL', required=False, type=float,
                                     default=0.6, help='Coverage for the longer sequence.')
    subparser_partition.add_argument('--evalue', required=False, type=float,
                                     default=1E-5, help='The evalue of aligner.')
    subparser_partition.add_argument('--aligner', required=False, type=str,
                                     default='diamond', choices=('diamond', 'blastp'), help='The aligner used to pairwise alignment.')
    subparser_partition.add_argument('--clust_method', required=False, type=str,
                                     default='mmseqs2', choices=('cdhit', 'mmseqs2'), help='The method used to cluster the genes.')
    subparser_partition.add_argument('--radius', required=False, type=int,
                                     default=3, help='The radius of search region.')
    subparser_partition.add_argument('--min_falen', '-m', required=False, type=check_min_falen,
                                     default=20, help='protein length of throw_away_sequences, at least 11')
    subparser_partition.add_argument('--gcode', required=False, type=check_gcode,
                                     default=11, help='The genetic code of your species.')
    subparser_partition.add_argument('--annot', required=False,
                                     action='store_true', default=False, help='Discard original annotation, and re-annote the genome privately using prodigal')

    subparser_partition.add_argument(
        '--verbose', '-v', required=False, action='store_true', default=False, help='Verbose output')
    subparser_partition.add_argument(
        '--debug', '-D', required=False, action='store_true', default=False, help='Debug mode. Note: very verbose')
    subparser_partition.add_argument(
        '--disable', required=False, action='store_true', default=False, help='Disable progress bar')

    subparser_partition.set_defaults(func=partition_portal)
