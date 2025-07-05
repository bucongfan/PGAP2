import os
import pickle
import random
import argparse
import numpy as np
import pandas as pd

from tqdm import tqdm
from loguru import logger
from functools import partial
from multiprocessing import Pool
from scipy.optimize import curve_fit
from argparse import ArgumentParser, _SubParsersAction

from pgap2.lib.basic import Basic
from pgap2.lib.pangp import PanGP
from pgap2.lib.pklcheck import PklCheck
from pgap2.lib.pangenome import pan_judger
from pgap2.utils.supply import sfw, tqdm_
from pgap2.utils.supply import set_verbosity_level
from pgap2.utils.draw import postprocess_draw, postprocess_draw_vector


def generate_tasks(strain_num, N):
    if strain_num <= N:
        tasks = [i
                 for i in range(1, strain_num + 1)]
    else:
        indices = np.linspace(1, strain_num, N, dtype=int)
        tasks = [i for i in indices]

    return tasks


def core_gene_model(x, C0, a, C_inf):
    return C0 * np.exp(-a * x) + C_inf


def pan_gene_model(x, k, gamma):
    return k * x**gamma


def stat_pan_core(sampling_list, pav, i_strain, strain_num):
    core_list = []
    pan_list = []
    for this_index in sampling_list:
        this_pav = pav[:, sorted(this_index)]
        pan_num = np.sum(np.any(this_pav != 0, axis=1))
        core_num = np.sum(np.all(this_pav != 0, axis=1))
        core_list.append(core_num)
        pan_list.append(pan_num)
        if i_strain == strain_num:
            for temp_i in range(1, strain_num):
                core_list.append(core_num)
                pan_list.append(pan_num)
    return core_list, pan_list


def calculate_rarefaction(i_strain, pangp: PanGP):

    strain_num = pangp.strain_num
    sampling_record = {'core': [], 'pan': []}
    for _ in range(pangp.repeat):
        sampling_list = pangp.sampling(i_strain)
        core_list, pan_list = stat_pan_core(
            sampling_list, pangp.pav, i_strain, strain_num)
        sampling_record['core'].append(core_list)
        sampling_record['pan'].append(pan_list)

    record = {'strain_num': i_strain}
    for key in ('core', 'pan'):
        mean_list = np.round(
            np.mean(np.array(sampling_record[key]), axis=0)).astype(int).tolist()
        record.update({key: mean_list})

    random.seed(i_strain)
    strain_rank = np.arange(strain_num)
    random.shuffle(strain_rank)
    new_clusters = calculate_new_clusters(pangp.pav, strain_rank)

    return {'strain_num': i_strain, 'profile': record, 'new_clusters': new_clusters, }


def log_pan_genome_state(popt_pan, ofh):
    k, gamma = popt_pan
    if gamma < 1:
        genome_status = "open"
        status_description = (
            "The pan-genome is open, meaning that the discovery of new genes will continue as more genomes are added."
        )
    else:
        genome_status = "closed"
        status_description = (
            "The pan-genome is closed, meaning that the discovery of new genes will saturate as more genomes are added."
        )
    logger.info(
        f"The pan-genome is {genome_status}. {status_description}")
    ofh.write(
        f"The pan-genome is {genome_status}. {status_description}\n")

    logger.info(
        f"The coefficient k is approximately {k:.2f}, indicating the initial addition rate of new gene clusters per genome.\nLarger values of k suggest a higher initial rate of gene discovery. ")
    ofh.write(
        f"The coefficient k is approximately {k:.2f}, indicating the initial addition rate of new gene clusters per genome.\nLarger values of k suggest a higher initial rate of gene discovery. ")


def fit_rerefaction(pan_rarefaction, outdir):

    x_data = pan_rarefaction[0]
    core_gene_samples = pan_rarefaction[1]
    pan_gene_samples = pan_rarefaction[2]

    for i in range(len(core_gene_samples)):
        if 'NA' in core_gene_samples[i]:
            core_gene_samples[i] = [
                int(_) for _ in core_gene_samples[i] if _ != 'NA']
        if 'NA' in pan_gene_samples[i]:
            pan_gene_samples[i] = [
                int(_) for _ in pan_gene_samples[i] if _ != 'NA']

    core_gene_means = [np.mean(samples) for samples in core_gene_samples]
    pan_gene_means = [np.mean(samples) for samples in pan_gene_samples]

    ofh = open(file=f'{outdir}/postprocess.curve_fit.txt', mode='w')

    logger.info("Fitting core gene data.")
    popt_core, _ = curve_fit(core_gene_model, x_data,
                             core_gene_means, p0=[2000, 0.1, 500])
    logger.info(f"Optimal parameters for core genes: {popt_core}")
    ofh.write(f'Optimal parameters for core genes: {popt_core}\n')

    logger.info("Fitting pan gene data.")
    popt_pan, _ = curve_fit(pan_gene_model, x_data,
                            pan_gene_means, p0=[1000, 0.8])
    logger.info(f"Optimal parameters for pan genes: {popt_pan}")
    ofh.write(f'Optimal parameters for pan genes: {popt_pan}\n')
    ofh.write(
        f'Core Gene Fit: y={popt_core[0]:.3f}*exp({popt_core[1]:.3f}*x)+ {popt_core[2]:.3f}\n')
    ofh.write(
        f'Pan Gene Fit: y={popt_pan[0]:.3f}*x^{popt_pan[1]:.3f} \n')
    log_pan_genome_state(popt_pan, ofh)

    # plt.figure(figsize=(8, 6))

    # logger.info("Plotting data and fitted curves.")
    # plt.errorbar(x_data, core_gene_means, yerr=[np.std(samples) for samples in core_gene_samples],
    #              fmt='s', markerfacecolor='none', markeredgecolor='b', markersize=5, ecolor='b', elinewidth=1,
    #              label='Core Genes', color='blue', capsize=5, linestyle='None')
    # x_fine = np.linspace(min(x_data), max(x_data), 100)
    # plt.plot(x_fine, core_gene_model(x_fine, *popt_core), 'b--',
    #          label=f'Core Gene Fit: y={popt_core[0]:.3f}*exp({popt_core[1]:.3f}*x) + {popt_core[2]:.3f}')

    # plt.errorbar(x_data, pan_gene_means, yerr=[np.std(samples) for samples in pan_gene_samples],
    #              fmt='s', markerfacecolor='none', markeredgecolor='r', markersize=5, ecolor='r', elinewidth=1,
    #              label='Pan Genes', color='red', capsize=5, linestyle='None')
    # plt.plot(x_fine, pan_gene_model(x_fine, *popt_pan), 'r--',
    #          label=f'Pan Gene Fit: y={popt_pan[0]:.3f}*x^{popt_pan[1]:.3f}')

    # plt.xlabel('Genome number')
    # plt.ylabel('Gene cluster number')
    # plt.legend()
    # plt.title('Core and Pan Genes Curve Fitting')
    # plt.grid(True)
    # plot_path = f'{outdir}/postprocess.rarefaction.pdf'
    # plt.savefig(plot_path)
    # logger.info(f"Saved the rarefaction plot to {plot_path}")


def calculate_new_clusters(pav, genome_order):

    selected_pav = np.zeros(pav.shape[0], dtype=bool)
    new_clusters = []

    for i in range(len(genome_order)):
        genome_idx = genome_order[i]

        current_genes = pav[:, genome_idx].astype(bool)

        new_genes = current_genes & ~selected_pav

        selected_pav |= current_genes

        new_clusters.append(int(np.sum(new_genes)))

    return new_clusters


def get_rarefaction(pav: pd.DataFrame, outdir: str, disable: bool = False, threads: int = 1, S: int = 100, N: int = 100, R: int = 3, K: int = 3, method: str = 'DG'):

    logger.info(f"Start rarefaction analysis with {threads} threads")

    pav = pav.to_numpy()
    pangp = PanGP(pav, S=S, R=R, K=K, method=method)
    logger.info(f"{pav.shape[1]} strains were found in file")
    logger.info(f"Loaded {pav.shape[0]} gene clusters")
    strain_num = pangp.strain_num
    if strain_num <= 2:
        logger.error(
            "The number of strains is less than 2, so rarefaction analysis cannot be performed.")
        raise ValueError(
            "The number of strains is less than 2, so rarefaction analysis cannot be performed.")
    if strain_num < N:
        logger.warning(
            f"You assigned {N} strains, but only {strain_num} strains were found.")
        logger.warning(f"All {strain_num} strains will be used.")
    elif strain_num >= N:
        logger.info(
            f"Loaded {strain_num} strains, But only {N} strains will be used.")

    if method == 'DG':
        logger.info(
            f"Sampling {S} strain from {S*K} combinations for each strain {R} times.")
    elif method == 'TR':
        logger.info(f"Sampling {S} strains {R} times.")

    tasks = generate_tasks(strain_num, N)
    # Use multiprocessing to handle the computation
    pool = Pool(threads)
    results = list(tqdm(pool.imap(partial(calculate_rarefaction, pangp=pangp), tasks), total=len(
        tasks), unit=' strain', desc=tqdm_.step(1), disable=disable))

    pool.close()
    pool.join()

    rarefaction_record = {}
    new_clusters_record = {i: [] for i in tasks if i != 1}
    with open(file=f'{outdir}/postprocess.rarefaction.tsv', mode='w') as ofh_1, open(file=f'{outdir}/postprocess.new_clusters.tsv', mode='w') as ofh_2:
        ofh_1.write('Type\tStrain\tSampling\n')
        ofh_2.write('Strain\tSampling\n')
        for result in results:
            i_strain = result['strain_num']
            record = result['profile']
            new_clusters = result['new_clusters']
            for i, new_cluster in enumerate(new_clusters):
                if i+1 in new_clusters_record:
                    new_clusters_record[i+1].append(new_cluster)

            rarefaction_record[i_strain] = {
                'core': record['core'], 'pan': record['pan']}

            for key in ('core', 'pan'):
                intact_key = 'Core genome' if key == 'core' else 'Pan genome'
                for sample in record[key]:
                    ofh_1.write(f'{intact_key}\t{i_strain}\t{sample}\n')
        for i, new_cluster in new_clusters_record.items():
            for value in new_cluster:
                ofh_2.write(f'{i}\t{value}\n')

    pan_genome_profile = [[], [], []]

    for i_strain in sorted(rarefaction_record.keys()):
        core = rarefaction_record[i_strain]['core']
        pan = rarefaction_record[i_strain]['pan']
        pan_genome_profile[0].append(i_strain)
        pan_genome_profile[1].append(core)
        pan_genome_profile[2].append(pan)

    new_clusters_profile = [[], []]
    for i_strain, new_clusters in new_clusters_record.items():
        new_clusters_profile[0].append(i_strain)
        new_clusters_profile[1].append(new_clusters)
    return pan_genome_profile, new_clusters_profile


def get_pan_group(pav: pd.DataFrame, outdir: str, step=0.05):
    pan_group = {'Strict_core': 0, 'Core': 0,
                 'Soft_core': 0, 'Shell': 0, 'Cloud': 0}
    para_dict = {'Strict_core': [[], []], 'Core': [[], []],
                 'Soft_core': [[], []], 'Shell': [[], []], 'Cloud': [[], []]}
    total_strain_num = len(pav.columns)
    total_clust_num = len(pav)
    x = [str(round(f, 2)) for f in np.arange(0, 1.01, step)]
    x_index = {round(float(f), 2): i for i, f in enumerate(x)}
    y = [0]*len(x)
    strain_num = 0
    ofh = open(file=f'{outdir}/postprocess.clust_strain_freq.tsv', mode='w')
    ofh.write('Clust\tStrain_num\tFreq\n')
    for index, row in pav.iterrows():
        strain_num = (row != 0).sum()
        if strain_num > total_strain_num:
            logger.error(
                f'strain_num: {strain_num} total_num: {total_strain_num} row: {row}')
        freq = strain_num/total_strain_num
        ofh.write(f'{index}\t{strain_num}\t{freq}\n')
        group_freq = round(int(freq/step)*step, 2)
        assert group_freq in x_index, logger.error(
            f'group_freq: {group_freq} not in x: {x}')
        y[x_index[group_freq]] += 1
        pangroup = pan_judger(strain_num, total_strain_num)
        pan_group[pangroup] += 1

        para_strain = (row > 1).sum()
        if para_strain > 0:
            para_gene = row[row > 1].sum()-para_strain
            para_dict[pangroup][0].append(int(para_strain))
            para_dict[pangroup][1].append(int(para_gene))
    ofh.close()
    with open(f'{outdir}/postprocess.pan_group_stat.tsv', mode='w') as ofh:
        ofh.write('Group\tCount\tProportion\n')
        total_clusters = sum(pan_group.values())
        for k, v in pan_group.items():
            if k == 'Strict_core':
                k = 'Strict core'
            elif k == 'Core':
                k = 'Core'
            elif k == 'Soft_core':
                k = 'Soft core'
            elif k == 'Shell':
                k = 'Shell'
            elif k == 'Cloud':
                k = 'Cloud'
            proportion = round(100*v/total_clusters, 2)
            ofh.write(f'{k}\t{v}\t{proportion}\n')

    with open(f'{outdir}/postprocess.para_stat.tsv', mode='w') as ofh:
        ofh.write('Group\tPara_strain\tPara_gene\n')
        for group, (para_strains, para_genes) in para_dict.items():
            for para_strain, para_gene in zip(para_strains, para_genes):
                ofh.write(f'{group}\t{para_strain}\t{para_gene}\n')

    return (x, y, pan_group), para_dict


def main(indir: str, outdir: str, nodraw: bool, threads: int, disable: bool = False, single_file: bool = False, S: int = 100, N: int = 500, R: int = 3, K: int = 3, method: str = 'DG'):
    logger.info(f'Loading the basic.pkl...')
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

    basic.load_pav(file=f'{indir}/pgap2.partition.gene_content.pav')
    basic.stat_from_detail_pav(
        file=f'{indir}/pgap2.partition.gene_content.detail.tsv')
    pan_group_freq, pan_para_stat = get_pan_group(basic.pav, outdir)
    pan_profile, new_clusters = get_rarefaction(
        basic.pav, outdir, disable=disable, threads=threads, S=S, N=N, R=R, K=K, method=method)
    logger.info("Starting rarefaction curve fitting.")
    fit_rerefaction(pan_profile, outdir)

    if nodraw:
        logger.info('Drawing is disabled')
    else:
        logger.info('Generate the interactive report...')
        html_report = postprocess_draw(target='stat', pan_profile=pan_profile, new_clusters=new_clusters,
                                       group_freq=pan_group_freq, pan_para_stat=pan_para_stat, basic=basic, outdir=outdir)
        logger.info(f'Interactive report at {html_report}')

        logger.info('Draw the vector report...')
        postprocess_draw_vector(
            target='stat', sfw=sfw.draw_post_stat, outdir=outdir, single_file=single_file)
        postprocess_draw_vector(
            target='profile', sfw=sfw.draw_post_profile, outdir=outdir, single_file=single_file)

        logger.info(f'Vector report at {outdir}:')
        if single_file:
            if os.path.exists(f'{outdir}/postprocess.para_stat.pdf'):
                logger.info(
                    f'[1/10] postprocess.stat_attrs_mean.pdf')
                logger.info(
                    f'[2/10] postprocess.stat_attrs_min.pdf')
                logger.info(
                    f'[3/10] postprocess.stat_attrs_var.pdf')
                logger.info(
                    f'[4/10] postprocess.stat_attrs_uni.pdf')
                logger.info(
                    f'[5/10] postprocess.pan_group_stat.pdf')
                logger.info(
                    f'[6/10] postprocess.clust_strain_freq.pdf')
                logger.info(
                    f'[7/10] postprocess.rarefaction.pdf')
                logger.info(
                    f'[8/10] postprocess.new_clusters.pdf')
                logger.info(
                    f'[9/10] postprocess.para_stat.pdf')
                logger.info(
                    f'[10/10] postprocess.para_stat_facet.pdf')
            else:
                logger.info(
                    f'[1/8] postprocess.stat_attrs_mean.pdf')
                logger.info(
                    f'[2/8] postprocess.stat_attrs_min.pdf')
                logger.info(
                    f'[3/8] postprocess.stat_attrs_var.pdf')
                logger.info(
                    f'[4/8] postprocess.stat_attrs_uni.pdf')
                logger.info(
                    f'[5/8] postprocess.pan_group_stat.pdf')
                logger.info(
                    f'[6/8] postprocess.clust_strain_freq.pdf')
                logger.info(
                    f'[7/8] postprocess.rarefaction.pdf')
                logger.info(
                    f'[8/8] postprocess.new_clusters.pdf')
        else:
            if os.path.exists(f'{outdir}/pgap2.postprocess_stat_para.pdf'):
                logger.info(
                    f'[1/3] pgap2.postprocess_stat.pdf')
                logger.info(
                    f'[2/3] pgap2.postprocess_stat_para.pdf')
                logger.info(
                    f'[3/3] pgap2.postprocess_profile.pdf')
            else:
                logger.info(
                    f'[1/2] pgap2.postprocess_stat.pdf')
                logger.info(
                    f'[2/2] pgap2.postprocess_profile.pdf')

    logger.success('Done')


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
    tqdm_.set_total_step(1)
    main(indir=args.indir, outdir=outdir, nodraw=args.nodraw,
         threads=args.threads, disable=args.disable, single_file=args.single_file,
         S=args.S, N=args.N, R=args.R, K=args.K, method=args.method)


def postprocess_portal(args):
    set_verbosity_level(args.outdir, args.verbose,
                        args.debug, 'postprocess_stat')
    if not args.nodraw:
        sfw.check_dependency('draw_post_profile')
        sfw.check_dependency('draw_post_stat')
    launch(args)


def stat_cmd(subparser: _SubParsersAction):
    subparser_postprocess: ArgumentParser = subparser.add_parser(
        'stat', help='Statistical analysis')
    subparser_postprocess.add_argument(
        '--indir', '-i', required=True, help='Input directory generated by partition step',)
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
        '--threads', '-t', required=False, default=1, help='threads used in parallel', type=int)
    subparser_postprocess.add_argument(
        '--verbose', '-V', required=False, action='store_true', default=False, help='Verbose output')
    subparser_postprocess.add_argument(
        '--debug', required=False, action='store_true', default=False, help='Debug mode. Note: very verbose')
    subparser_postprocess.add_argument(
        '--disable', required=False, action='store_true', default=False, help='Disable progress bar')
    subparser_postprocess.set_defaults(func=postprocess_portal)
