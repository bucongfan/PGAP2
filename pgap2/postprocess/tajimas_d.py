import os
import argparse
import pickle

from loguru import logger
from argparse import ArgumentParser, _SubParsersAction

from pgap2.lib.basic import Basic
from pgap2.lib.pklcheck import PklCheck
from pgap2.lib.phylogeny import Phylogeny
from pgap2.utils.supply import sfw, tqdm_
from pgap2.utils.supply import set_verbosity_level

"""
Statistical analysis of Tajima's D based on gene content.
Tajima's D is a measure of genetic diversity within a population.
This module provides a command-line interface for running the Tajima's D test on gene clusters.

input:
- indir: Input directory containing the necessary files for the Tajima's D test.
- outdir: Output directory where the results will be saved.
params:
- threads: Number of threads to use for parallel processing.
- disable: If set, disables the progress bar.
- msa_method: Method for multiple sequence alignment.
- para_strategy: Strategy for handling paralogs.
- clusts: A file containing one line per cluster for which Tajima's D test needs to be performed. By default, all clusters will be used.
- step: Specifies the step at which to terminate the workflow.
- add_paras: Additional parameters for specific steps in the workflow
output:
- a TSV file containing the results of the Tajima's D test for each cluster.
TODO:
- Add visualization capabilities for the Tajima's D results.
"""


def main(indir: str, outdir: str, nodraw: bool, threads: int, disable: bool = False, clusts: str = '', step: int = 5, para_strategy: str = 'best', msa_method: str = 'mafft', add_paras: list = []):
    detail_file = f'{indir}/pgap2.partition.gene_content.detail.tsv'
    fa_file = f'{indir}/total.involved_annot.tsv'

    with open(f'{indir}/basic.pkl', 'rb') as fh:
        previous: PklCheck = pickle.load(fh)
        decode_status = previous.decode()
        if decode_status:
            basic: Basic = previous.data_dump('basic')
            for rec in basic.dumper():
                logger.info(rec)
        else:
            logger.error('Failed to decode the basic.pkl')
            raise ValueError('Failed to decode the basic.pkl')
    logger.info(f'Reading the pangenome information...')
    basic.load_used_clusters(clusts=clusts, file=detail_file)
    basic.phylogeny_from_detail_pav(
        file=detail_file, core_thre=0, also_pan=False, para_strategy=para_strategy)
    basic.phylogeny_from_id(file=fa_file)
    phy = Phylogeny(basic=basic, outdir=outdir, threads=threads, disable=disable,
                    msa_method=msa_method, tree_method=None,
                    fastbaps_levels=None, fastbaps_prior=None, add_paras=add_paras)
    for step_i in range(step+1):
        if step_i <= 4:
            phy.start_at(step_i)
        elif step_i == 5:
            phy.start_at(11)
        else:
            logger.error(f'Invalid step number: {step_i}')
            raise ValueError(f'Invalid step number: {step_i}')
    logger.info('Dumping results...')
    phy.dump_results()
    logger.success('All steps done')

    # logger.info('Drawing...')
    # html_report = postprocess_draw(target='phylogeny', outdir=outdir)
    # logger.info(f'Report at {html_report}')
    # logger.success('Done')


def launch(args: argparse.Namespace):
    outdir = os.path.abspath(args.outdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    assert os.path.exists(
        f'{args.indir}/basic.pkl'), logger.error(f'basic.pkl not found in {args.indir}, the input dir should be the output dir of the main step')
    assert os.path.exists(
        f'{args.indir}/pgap2.partition.gene_content.detail.tsv'), logger.error(f'pgap2.partition.gene_content.detail.tsv not found in {args.indir}, the input dir should be the output dir of the main step')
    assert os.path.exists(
        f'{args.indir}/pgap2.partition.gene_content.pav'), logger.error(f'pgap2.partition.gene_content.pav not found in {args.indir}, the input dir should be the output dir of the main step')
    main(indir=args.indir, outdir=outdir, nodraw=args.nodraw,
         threads=args.threads, disable=args.disable,
         step=args.step, para_strategy=args.para_strategy,
         clusts=args.clusts,
         msa_method=args.msa_method,
         add_paras=args.add_paras)


def postprocess_portal(args):
    set_verbosity_level(args.outdir, args.verbose,
                        args.debug, 'postprocess_tree')
    tqdm_.set_total_step(args.step)

    if args.step >= 2:
        if args.msa_method == 'muscle':
            sfw.check_dependency('muscle')
        elif args.msa_method == 'mafft':
            sfw.check_dependency('mafft')
        elif args.msa_method == 'tcoffee':
            sfw.check_dependency('tcoffee')
    if args.step >= 4:
        sfw.check_dependency('clipkit')

    launch(args)


def tajimas_d_cmd(subparser: _SubParsersAction):
    subparser_postprocess: ArgumentParser = subparser.add_parser(
        'tajimas_d', help="Workflow for Tajima's D test", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparser_postprocess.add_argument(
        '--indir', '-i', required=True, help='Input directory generated by partition',)
    subparser_postprocess.add_argument(
        '--outdir', '-o', required=False, help='Output directory', default='.',)
    subparser_postprocess.add_argument('--nodraw', required=False,
                                       action='store_true', default=False, help='Only output flat file, but no drawing plot')
    subparser_postprocess.add_argument(
        '--verbose', '-V', required=False, action='store_true', default=False, help='Verbose output')
    subparser_postprocess.add_argument(
        '--debug', required=False, action='store_true', default=False, help='Debug mode. Note: very verbose')
    subparser_postprocess.add_argument(
        '--disable', required=False, action='store_true', default=False, help='Disable progress bar')
    subparser_postprocess.add_argument(
        '--threads', '-t', required=False, default=1, help='threads used in parallel', type=int)
    subparser_postprocess.add_argument('--msa_method', required=False, default='mafft', choices=(
        'mafft', 'muscle', 'tcoffee'), help='The method of multiple sequence alignment.')
    subparser_postprocess.add_argument('--para_strategy', required=False, default='best', choices=('drop', 'best'),
                                       help='The strategy of paralog including cluster. best: keep the best one; drop: drop all paralogs contained clusters.')
    subparser_postprocess.add_argument(
        '--clusts', required=False, help="A file containing one line per cluster for which Tajima's D test needs to be performed. By default, all clusters will be used.")
    subparser_postprocess.add_argument('--step', required=False, default=5, help='''Terminate at this step. 1. Extract core cds and prot.
                                       2. Multiple sequence alignment using --msa_method assigned.
                                       3. Codon alignment: Trans multiple protein alignment to corresponding nucleotide alignment.
                                       4. Trim alignment using ClipKit.
                                       5. Tajima's D Test.
                                       ''', type=int, choices=[1, 2, 3, 4, 5])
    subparser_postprocess.add_argument(
        '--add_paras', action='append', help='Add additional parameters in the step. Format like "Step number:Parameters.", such as: 4:-g 0.8')

    subparser_postprocess.set_defaults(func=postprocess_portal)
