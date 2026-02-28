"""
CLI definitions for the 'minus' subcommand.
"""
import os
import pickle
import argparse

from argparse import ArgumentParser, _SubParsersAction

from pgap2.utils.supply import sfw, tqdm_
from pgap2.utils.supply import set_verbosity_level
from pgap2.lib.pklcheck import PklCheck


def launch(args: argparse.Namespace):
    from pgap2.utils.partition_minus import main
    main(prev_dir=os.path.abspath(args.prev_dir),
         removal_file=os.path.abspath(args.removal_file),
         outdir=os.path.abspath(args.outdir),
         threads=args.threads,
         disable=args.disable)


def minus_portal(args):
    set_verbosity_level(args.outdir, args.verbose,
                        args.debug, 'partition')

    sfw.check_dependency("mcl")
    sfw.check_dependency("mcxdeblast")

    # Load params to check if retrieve is needed
    prev_dir = os.path.abspath(args.prev_dir)
    with open(f'{prev_dir}/basic.pkl', 'rb') as fh:
        pkl: PklCheck = pickle.load(fh)
        basic = pkl.data_dump('basic')
        params = basic.params

    retrieve = params.get('retrieve', False)

    if retrieve:
        sfw.check_dependency("miniprot")
        sfw.check_dependency("seqtk")
        tqdm_.set_total_step(12)
    else:
        tqdm_.set_total_step(7)
    launch(args)


def minus_cmd(subparser: _SubParsersAction):

    subparser_minus: ArgumentParser = subparser.add_parser(
        'minus', help='Remove strains from an existing PGAP2 result',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparser_minus.add_argument(
        '--prev_dir', '-p', required=True,
        help='Previous PGAP2 result directory (output from pgap2 main, add, or merge).')
    subparser_minus.add_argument(
        '--removal_file', '-r', required=True,
        help='Text file listing strain names to remove (one per line).')
    subparser_minus.add_argument(
        '--outdir', '-o', required=False, default='.',
        help='Output directory for results after removal.')
    subparser_minus.add_argument(
        '--threads', '-t', required=False, default=1, type=int,
        help='Threads used in parallel.')
    subparser_minus.add_argument(
        '--verbose', '-v', required=False,
        action='store_true', default=False,
        help='Verbose output.')
    subparser_minus.add_argument(
        '--debug', '-D', required=False,
        action='store_true', default=False,
        help='Debug mode. Note: very verbose.')
    subparser_minus.add_argument(
        '--disable', required=False,
        action='store_true', default=False,
        help='Disable progress bar.')

    subparser_minus.set_defaults(func=minus_portal)
