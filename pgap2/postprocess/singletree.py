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
Post-process the single copy tree results.
This module provides a command-line interface for running the single copy tree workflow, which includes
multiple steps such as multiple sequence alignment, tree construction.
It allows users to specify various parameters for the analysis, such as the method of multiple sequence alignment
and the method of tree construction.

input:
- indir: Input directory containing the necessary files for the single copy tree workflow.
- outdir: Output directory where the results will be saved.
params:
- nodraw: If set, the workflow will not generate any graphical output.
- threads: Number of threads to use for parallel processing.
- disable: If set, disables the progress bar.
- also_pan: If set, includes pan-genome alignment in the output.
- core_thre: Core threshold for tree construction.
- step: Specifies the step at which to terminate the workflow.
- para_strategy: Strategy for handling paralogs.
- msa_method: Method for multiple sequence alignment.
- tree_method: Method for tree construction.
- add_paras: Additional parameters for specific steps in the workflow.
output:
- Phylogenetic tree and population structure analysis results.

TODO:
- Add visualization capabilities for the phylogenetic tree.
"""


def main(indir: str, outdir: str, nodraw: bool, threads: int, disable: bool = False, also_pan: bool = False, core_thre: float = 0.95, step: int = 9, para_strategy: str = 'best', msa_method: str = 'mafft', tree_method: str = 'fasttree', add_paras: list = []):
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
    basic.load_used_clusters(clusts=None, file=detail_file)
    basic.phylogeny_from_detail_pav(
        file=detail_file, core_thre=core_thre, also_pan=also_pan, para_strategy=para_strategy)
    basic.phylogeny_from_id(file=fa_file)
    phy = Phylogeny(basic=basic, outdir=outdir, threads=threads, disable=disable,
                    msa_method=msa_method, tree_method=tree_method,
                    fastbaps_levels=None, fastbaps_prior=None, add_paras=add_paras)
    for step_i in range(step):
        step_i += 1
        phy.start_at(step_i)
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
         threads=args.threads, disable=args.disable, also_pan=args.also_pan,
         core_thre=args.core_threshold, step=args.step, para_strategy=args.para_strategy,
         msa_method=args.msa_method,
         tree_method=args.tree_method, add_paras=args.add_paras)


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
    if args.step >= 6:
        if args.tree_method == 'fasttree':
            sfw.check_dependency('fasttree')
        elif args.tree_method == 'raxml':
            sfw.check_dependency('raxml')
        elif args.tree_method == 'iqtree':
            sfw.check_dependency('iqtree')

    launch(args)


def singletree_cmd(subparser: _SubParsersAction):
    subparser_postprocess: ArgumentParser = subparser.add_parser(
        'singletree', help='Workflow for bayesian analysis of population structure', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
        '--threads', '-t', required=False, default=8, help='threads used in parallel', type=int)
    subparser_postprocess.add_argument('--msa_method', required=False, default='mafft', choices=(
        'mafft', 'muscle', 'tcoffee'), help='The method of multiple sequence alignment.')
    subparser_postprocess.add_argument('--tree_method', required=False, default='fasttree', choices=(
        'fasttree', 'raxml', 'iqtree'), help='The method of tree building.')
    subparser_postprocess.add_argument('--para_strategy', required=False, default='best', choices=('drop', 'best'),
                                       help='The strategy of paralog including cluster. best: keep the best one; drop: drop all paralogs contained clusters.')
    subparser_postprocess.add_argument(
        '--also_pan', required=False, default=False, action='store_true', help='Also output pan genome alignment')

    subparser_postprocess.add_argument('--core_threshold', required=False, default=0.95,
                                       help='Core threshold to extract for tree construction', type=float)
    subparser_postprocess.add_argument('--step', required=False, default=6, help='''Terminate at this step. 1. Extract core cds and prot.
                                       2. Multiple sequence alignment using --msa_method assigned.
                                       3. Codon alignment: Trans multiple protein alignment to corresponding nucleotide alignment.
                                       4. Trim alignment using ClipKit.
                                       5. Core alignment concatenation
                                       6. Phylogeny tree construction using tree_method assigned software.
                                       ''', type=int, choices=[1, 2, 3, 4, 5, 6])
    subparser_postprocess.add_argument(
        '--add_paras', action='append', help='Add additional parameters in the step. Format like "Step number:Parameters.", such as: 6:-pseudo')

    subparser_postprocess.set_defaults(func=postprocess_portal)
