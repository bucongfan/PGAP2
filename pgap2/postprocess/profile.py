import os
import argparse
import pandas as pd

from loguru import logger
from argparse import ArgumentParser, _SubParsersAction

from pgap2.utils.supply import sfw, tqdm_
from pgap2.utils.tools import detect_separator, is_numeric_pd
from pgap2.utils.supply import set_verbosity_level
from pgap2.utils.draw import postprocess_draw, postprocess_draw_vector
from pgap2.postprocess.stat import get_pan_group, get_rarefaction, fit_rerefaction


def main(pav_file: str, outdir: str, nodraw: bool, threads: int = 1, disable: bool = False, single_file: bool = False, S: int = 100, N: int = 500, R: int = 3, K: int = 3, method: str = 'DG'):

    logger.info(f'Reading the pangenome information...')
    sep = detect_separator(pav_file)
    pav = pd.read_csv(pav_file, sep=sep, index_col=0)
    is_number = is_numeric_pd(pav)
    if not is_number:
        logger.warning(
            'The pav matrix contains non-numeric values, counting the number of ";" in each cell')
        pav = pav.map(lambda x: 0 if pd.isna(x) else len(str(x).split(';')))

    pan_profile, new_clusters = get_rarefaction(
        pav, outdir, disable=disable, threads=threads, S=S, N=N, R=R, K=K, method=method)
    group_freq, pan_para_stat = get_pan_group(pav, outdir)

    logger.info("Starting rarefaction curve fitting.")
    fit_rerefaction(pan_profile, outdir)

    if nodraw:
        logger.info('Drawing is disabled')
    else:
        logger.info('Drawing...')
        html_report = postprocess_draw(
            target='profile', pan_profile=pan_profile, new_clusters=new_clusters, group_freq=group_freq, pan_para_stat=pan_para_stat, outdir=outdir)
        logger.info(f'Report at {html_report}')

        logger.info('Draw the vector report...')
        postprocess_draw_vector(
            target='profile', sfw=sfw.draw_post_profile, outdir=outdir, single_file=single_file)
        logger.info(f'Vector report at {outdir}:')
        if single_file:
            logger.info(
                f'[1/6] postprocess.pan_group_stat.pdf')
            logger.info(
                f'[2/6] postprocess.clust_strain_freq.pdf')
            logger.info(
                f'[3/6] postprocess.rarefaction.pdf')
            logger.info(
                f'[4/6] postprocess.new_clusters.pdf')
            logger.info(
                f'[5/6] postprocess.para_stat.pdf')
            logger.info(
                f'[6/6] postprocess.para_stat_facet.pdf')
        else:
            logger.info(
                f'[1/2] pgap2.postprocess_stat_para.pdf')
            logger.info(
                f'[2/2] pgap2.postprocess_profile.pdf')

        logger.info(f'Vector report at {outdir}:')

    logger.success('Done')


def launch(args: argparse.Namespace):
    outdir = os.path.abspath(args.outdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    tqdm_.set_total_step(1)
    assert os.path.exists(
        f'{args.pav}'), logger.error(f'pav file not found: {args.pav}')
    tqdm_.set_total_step(1)
    main(pav_file=args.pav, outdir=outdir, nodraw=args.nodraw,
         disable=args.disable, threads=args.threads, single_file=args.single_file,
         S=args.S, N=args.N, R=args.R, K=args.K, method=args.method)


def postprocess_portal(args):
    set_verbosity_level(args.outdir, args.verbose,
                        args.debug, 'postprocess_profile')
    sfw.check_dependency('draw_post_profile')
    launch(args)


def profile_cmd(subparser: _SubParsersAction):
    subparser_postprocess: ArgumentParser = subparser.add_parser(
        'profile', help='To generate the pangenome profile using PAV matrix, it is the subset of [stat] module', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparser_postprocess.add_argument(
        '--pav', '-p', help='The pav matrix', required=True)
    subparser_postprocess.add_argument(
        '--outdir', '-o', required=False, help='Output directory', default='.',)
    subparser_postprocess.add_argument('--nodraw', required=False,
                                       action='store_true', default=False, help='Only output flat file, but no drawing plot')
    subparser_postprocess.add_argument(
        '--single_file', '-s', action='store_true', default=False, help='Output each vector plot as a single file')
    subparser_postprocess.add_argument(
        '--S', required=False, type=int, default=100, help='Number of strains to sample in each bin')
    subparser_postprocess.add_argument(
        '--N', required=False, type=int, default=100, help='Number of bins to sample')
    subparser_postprocess.add_argument(
        '--R', required=False, type=int, default=1, help='Number of repeats for every sample')
    subparser_postprocess.add_argument(
        '--K', required=False, type=int, default=3, help='Inflated cof for DG method')
    subparser_postprocess.add_argument(
        '--method', required=False, choices=['TR', 'DG'], default='TR', help='Sampling method. TR: total random; DG: diversity guided')
    subparser_postprocess.add_argument(
        '--threads', '-t', required=False, type=int, default=4, help='Number of threads')
    subparser_postprocess.add_argument(
        '--disable', required=False, action='store_true', default=False, help='Disable progress bar')
    subparser_postprocess.add_argument(
        '--verbose', required=False, action='store_true', default=False, help='Verbose output')
    subparser_postprocess.add_argument(
        '--debug', required=False, action='store_true', default=False, help='Debug mode. Note: very verbose')
    subparser_postprocess.set_defaults(func=postprocess_portal)
