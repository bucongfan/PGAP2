"""
CLI definitions for the preprocess subcommand.
"""
import os
import argparse

from argparse import ArgumentParser, _SubParsersAction

from pgap2.utils.supply import sfw, tqdm_
from pgap2.utils.supply import set_verbosity_level
from pgap2.utils.graph_utils import check_min_falen, check_gcode


def launch(args: argparse.Namespace):
    from pgap2.utils.preprocess import main
    indir = os.path.abspath(args.indir)
    outdir = os.path.abspath(args.outdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    main(indir=indir, outdir=outdir,
         orth_id=args.orth_id, para_id=args.para_id, dup_id=args.dup_id, accurate=args.accurate,
         id_attr_key=args.id_attr_key, type_filter=args.type_filter, max_targets=args.max_targets,
         LD=args.LD, AS=args.AS, AL=args.AL,
         evalue=args.evalue,
         aligner=args.aligner, clust_method=args.clust_method,
         coverage=0.98,
         #  coverage=args.coverage,
         nodraw=args.nodraw, single_file=args.single_file,
         marker_file=args.marker, ani_thre=args.ani_thre,
         annot=args.annot, threads=args.threads, gcode=args.gcode,
         disable=args.disable, retrieve=args.retrieve, falen=args.min_falen,
         exclude_outlier=args.exclude_outlier,)
    return 0


def preprocess_portal(args):
    set_verbosity_level(args.outdir, args.verbose, args.debug, 'preprocess')
    if not args.nodraw:
        sfw.check_dependency('draw_prep')

    if args.clust_method == 'mmseqs2':
        sfw.check_dependency("mmseqs2")
    elif args.clust_method == 'cdhit':
        sfw.check_dependency("cdhit")

    if args.aligner == 'diamond':
        sfw.check_dependency("diamond")
    elif args.aligner == 'blast':
        sfw.check_dependency("blastp")
        sfw.check_dependency("makeblastdb")

    if args.retrieve:
        sfw.check_dependency("miniprot")
        sfw.check_dependency("seqtk")
    if args.annot:
        sfw.check_dependency("prodigal")

    tqdm_.set_total_step(5)
    launch(args)


def preprocess_cmd(subparser: _SubParsersAction):

    subparser_preprocess: ArgumentParser = subparser.add_parser(
        'prep', help='Preprocess the input files', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparser_preprocess.add_argument(
        '--indir', '-i', required=True, help='Input file contained, same prefix seems as the same strain.',)
    subparser_preprocess.add_argument(
        '--outdir', '-o', required=False, help='Output directory', default='.',)
    subparser_preprocess.add_argument('--dup_id', required=False, type=float,
                                      default=0.99, help='The maximum identity between the most recent duplication envent.')
    subparser_preprocess.add_argument('--orth_id', required=False, type=float,
                                      default=0.98, help='The maximum identity between the most similar panclusters, 0 means automatic selection.')
    subparser_preprocess.add_argument('--para_id', required=False, type=float,
                                      default=0.7, help='Use this identity as the paralogous identity, 0 means automatic selection.')
    subparser_preprocess.add_argument("--type-filter", required=False, type=str,
                                      default='CDS', help="Only for gff file as input, feature type (3rd column) to include, Only lines matching these types will be processed.")
    subparser_preprocess.add_argument("--id-attr-key", required=False, type=str,
                                      default='ID', help="Only for gff file as input, Attribute key to extract from the 9th column as the record ID (e.g., 'ID', 'gene', 'locus_tag').")
    # subparser_preprocess.add_argument('--coverage', required=False, type=float,
    #                                   default=0.98, help='The least coverage of each gene.')
    subparser_preprocess.add_argument('--min_falen', '-l', required=False, type=check_min_falen,
                                      default=20, help='protein length of throw_away_sequences, at least 11')
    subparser_preprocess.add_argument('--accurate', '-a', required=False,
                                      action='store_true', default=False, help='Apply bidirection check for paralogous gene partition.')
    subparser_preprocess.add_argument('--max_targets', '-k', required=False, type=int,
                                      default=2000, help='The maximum targets for each query in alignment. Improves accuracy for large-scale analyses, but increases runtime and memory usage.')
    subparser_preprocess.add_argument('--LD', required=False, type=float,
                                      default=0.6, help='Minimum gene length difference proportion between two genes.')
    subparser_preprocess.add_argument('--AS', required=False, type=float,
                                      default=0.6, help='Coverage for the shorter sequence.')
    subparser_preprocess.add_argument('--AL', required=False, type=float,
                                      default=0.6, help='Coverage for the longer sequence.')
    subparser_preprocess.add_argument('--evalue', required=False, type=float,
                                      default=1e-5, help='The evalue of aligner.')
    subparser_preprocess.add_argument('--aligner', required=False, type=str,
                                      default='diamond', choices=('diamond', 'blast'), help='The aligner used to pairwise alignment.')
    subparser_preprocess.add_argument('--clust_method', required=False, type=str,
                                      default='mmseqs2', choices=('cdhit', 'mmseqs2'), help='The method used to cluster the genes.')
    subparser_preprocess.add_argument('--marker', required=False,
                                      help='Assigned darb or outlier strain used to filter the input. See detail in marker.cfg in the main path', default=None)
    subparser_preprocess.add_argument('--ani_thre', required=False, type=float,
                                      help='Expect ani threshold', default=95)
    subparser_preprocess.add_argument('--annot', required=False,
                                      action='store_true', default=False, help='Discard original annotation, and re-annote the genome privately using prodigal')
    subparser_preprocess.add_argument('--retrieve', required=False,
                                      action='store_true', default=False, help='Retrieving gene that may lost with annotations')
    subparser_preprocess.add_argument('--gcode', required=False, type=check_gcode,
                                      default=11, help='The genetic code of your species. Default is [11] (bacteria).')
    subparser_preprocess.add_argument('--exclude_outlier', required=False,
                                      action='store_true', default=False, help='Exclude outlier strains from the output pkl. When set, low-quality or dissimilar strains will be removed from the pangenome data. By default, all strains are kept in the pkl (outliers are flagged but not removed).')
    subparser_preprocess.add_argument('--nodraw', required=False,
                                      action='store_true', default=False, help='Only output flat file, but no drawing plot')
    subparser_preprocess.add_argument(
        '--single_file', '-s', action='store_true', default=False, help='Output each vector plot as a single file')
    subparser_preprocess.add_argument(
        '--verbose', '-V', required=False, action='store_true', default=False, help='Verbose output')
    subparser_preprocess.add_argument(
        '--debug', '-D', required=False, action='store_true', default=False, help='Debug mode. Note: very verbose')
    subparser_preprocess.add_argument(
        '--disable', required=False, action='store_true', default=False, help='Disable progress bar')
    subparser_preprocess.add_argument(
        '--threads', '-t', required=False, default=1, help='threads used in parallel', type=int)
    subparser_preprocess.set_defaults(func=preprocess_portal)
