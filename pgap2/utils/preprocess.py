
import os
import re
import pickle
import shutil
import argparse
import numpy as np
import networkx as nx

from Bio import SeqIO
from tqdm import tqdm
from loguru import logger
from collections import defaultdict
from argparse import ArgumentParser, _SubParsersAction

from pgap2.lib.species import Species
from pgap2.lib.pangenome import Pangenome
from pgap2.lib.tree import Tree
from pgap2.lib.pklcheck import PklCheck

from pgap2.utils.supply import sfw, tqdm_
from pgap2.utils.supply import set_verbosity_level
from pgap2.utils.generate_tree import generate_tree
from pgap2.utils.tools import check_min_falen, check_gcode
from pgap2.utils.data_loader import file_parser, get_file_dict
from pgap2.utils.draw import preprocess_draw, preprocess_draw_vector

"""
preprocsess.py
This module provides functionality for preprocessing genomic data, including loading marker files, indexing strains,
computing average nucleotide identity (ANI), and identifying outlier genomes.
It includes functions for genome feature statistics, pooling genomes, creating symlinks, and preprocessing output
to generate reports and visualizations.

input:
- marker_file: Path to the marker file containing strain information.
- ani: Expected average nucleotide identity threshold.
- outdir: Output directory for storing results.

output:
- Preprocessed data including ANI calculations, outlier identification, and genome statistics.
- postprocess.ANI.pdf: Vector plot of ANI results.
- postprocess.gene_number.pdf: Vector plot of gene number statistics.
- postprocess.half_core.pdf: Vector plot of half-core gene statistics.
- postprocess.proportion.pdf: Vector plot of base composition proportions.
- postprocess.gene_code.pdf: Vector plot of gene code statistics.
- preprocess.stat.tsv: Tab-separated values file containing strain statistics.
- preprocess.gene_code.csv: CSV file containing gene code statistics.
- html_report: HTML report summarizing preprocessing results.
- refined_input: Directory containing symlinks to passed high-quality input files.
- preprocess.pkl: Pickle file containing preprocessed data for future re-use.
"""


def trim_path(ani_result, prefix, suffix):
    trim_result = {}
    for query in ani_result:
        trim_query = re.match(f'^{prefix}(.*)\.{suffix}$', query).group(1)
        trim_query = trim_query.replace('/', '')
        trim_result[trim_query] = {}
        for reference in ani_result[query]:
            trim_reference = re.match(
                f'^{prefix}(.*).{suffix}$', reference).group(1)
            trim_reference = trim_reference.replace('/', '')
            trim_result[trim_query].update(
                {trim_reference: ani_result[query][reference]})
    return trim_result


def genome_feature_stat(genome):
    BaseSum, chr_count = 0, 0
    no_c, no_g, no_a, no_t, no_n = 0, 0, 0, 0, 0
    for id in genome:
        record = genome[id]
        this_len = len(record.seq)
        BaseSum += this_len
        chr_count += 1
        seq = record.seq.lower()
        no_c += seq.count('c')
        no_g += seq.count('g')
        no_a += seq.count('a')
        no_t += seq.count('t')
        no_n += seq.count('n')
    return (chr_count, BaseSum, {'a': no_a, 't': no_t, 'c': no_c, 'g': no_g, 'n': no_n})


def pool_genome(strain_dict, outdir):
    query = []
    for strain in strain_dict:
        genome_record = strain_dict[strain].genome
        genome_file = f"{outdir}/{strain}.fasta"
        SeqIO.write([genome_record[id]
                    for id in genome_record], genome_file, "fasta")
        query.append(genome_file)
        logger.debug(f"{strain} -- {genome_file}")
    return {'fafile': query}


def create_symlink(source, target):
    target_dir = target
    base_name = os.path.basename(source)
    target = f'{target_dir}/{base_name}'
    if not os.path.exists(target):
        os.symlink(source, target)
    else:
        logger.warning(f"Soft link {target} already exists.")


def preprocess_output(darb_dict, outdir, sp: Species, file_dict: dict):
    refined_input = f'{outdir}/refined_input'
    if not os.path.exists(refined_input):
        os.mkdir(refined_input)
    else:
        logger.warning(f"{refined_input} already exists, will overwrite it.")
        shutil.rmtree(refined_input)
        os.mkdir(refined_input)

    with open(f'{outdir}/preprocess.gene_code.csv', 'w') as fh:
        headers = []
        for i in range(len(sp.strain_dict)):
            strain_name = sp.strain_dict[i].strain_name
            headers.append(strain_name)
        header = ','.join(headers)
        fh.write(f"Gene_code,{header}\n")
        for gene_code, value_dict in sp.gene_code.items():
            counts = []
            for i in range(len(sp.strain_dict)):
                if i in value_dict:
                    counts.append(value_dict[i])
                else:
                    counts.append(0)
            count = ','.join([str(count) for count in counts])
            fh.write(f"{gene_code},{count}\n")

    with open(f'{outdir}/preprocess.stat.tsv', 'w') as fh:
        fh.write("strain\tcontig_num\ttotal_gene_num\tgene_incomplete\thalf_core\tsingle_cloud\tfalen\tA|T|C|G\tani\tis_darb\tis_outlier_ani\tis_outlier_gene\n")
        darb = sp.get_darb()
        genome_attr = sp.genome_attr
        for strain in range(len(sp.strain_dict)):
            genome_attr = sp.genome_attr[strain]
            total_gene_num = genome_attr['total_gene_num']
            value_incomplete = genome_attr['value_incomplete']
            half_core = genome_attr['half_core']
            single_cloud = genome_attr['single_cloud']
            genome_len = genome_attr['genome_len']
            chr_count = genome_attr['chr_count']
            base_a = genome_attr['content']['a']
            base_t = genome_attr['content']['t']
            base_c = genome_attr['content']['c']
            base_g = genome_attr['content']['g']

            is_darb = 1 if strain == darb else 0
            ani = sp.get_ani(strain) if not is_darb else 100
            is_outlier_ani = 1 if strain in darb_dict['ani'] else 0
            is_outlier_gene = 1 if strain in darb_dict['single_gene'] else 0
            is_outlier = True if is_outlier_ani or is_outlier_gene else False

            strain_name = sp.strain_dict[strain].strain_name
            fh.write(f"{strain_name}\t{chr_count}\t{total_gene_num}\t{value_incomplete}\t{half_core}\t{single_cloud}\t{genome_len}\t{base_a}|{base_t}|{base_c}|{base_g}\t{ani}\t{is_darb}\t{is_outlier_ani}\t{is_outlier_gene}\n")
            if not is_outlier:
                for _, file in file_dict[strain_name].items():
                    try:
                        create_symlink(source=file,
                                       target=refined_input)
                    except FileExistsError:
                        pass


def stat_core(tree: Tree, pg: Pangenome):
    STRAIN_NUM = pg.strain_num
    half_num = int(STRAIN_NUM/2)
    symbol_strain_core_stat = defaultdict(int)
    symbol_strain_single_stat = defaultdict(int)
    for nodes in tqdm(nx.connected_components(tree.raw_distance_graph), desc=tqdm_.step(3), unit=' cluster', disable=pg.disable_tqdm):
        strains = set()
        has_para = False
        for node in nodes:
            this_strains = tree.orth_identity_tree.nodes[node]['strains']
            if strains & this_strains:
                has_para = True
                break
            else:
                strains |= this_strains
        if has_para:
            continue
        if len(strains) > half_num:
            for strain in strains:
                symbol_strain_core_stat[int(strain)] += 1
        elif len(strains) == 1:
            for strain in strains:
                symbol_strain_single_stat[int(strain)] += 1

    return symbol_strain_core_stat, symbol_strain_single_stat


def find_outlier_from_cloud(genome_attr):
    data = [v['single_cloud'] for k, v in genome_attr.items()]
    Q1 = np.percentile(data, 1)
    Q3 = np.percentile(data, 99)
    IQR = Q3 - Q1
    # not use lower bound because the lower bound is seems to be a good feature
    upper_bound = Q3 + 1.5 * IQR
    outlier = [k for k, v in genome_attr.items() if v['single_cloud']
               > upper_bound]
    return outlier


def main(indir, outdir, orth_id, para_id, dup_id, id_attr_key, type_filter, evalue, aligner, clust_method, accurate, coverage, nodraw, single_file, LD, AS, AL, marker_file, ani_thre, annot, threads, disable, retrieve, falen, gcode):
    logger.debug(f'----------------')
    file_dict = get_file_dict(indir=indir)
    decode_status = False
    if os.path.exists(f'{outdir}/preprocess.pkl'):
        logger.info(f'Found {outdir}/preprocess.pkl')
        logger.info(f'Loding...')
        with open(f'{outdir}/preprocess.pkl', 'rb') as fh:
            previous: PklCheck = pickle.load(fh)
            logger.info(f'Check the previous file parameters...')
            decode_status = previous.decode(
                orth_id=orth_id, para_id=para_id, dup_id=dup_id, id_attr_key=id_attr_key, type_filter=type_filter, accurate=accurate, coverage=coverage, falen=falen, annot=annot, retrieve=retrieve, evalue=evalue, aligner=aligner, clust_method=clust_method, LD=LD, AS=AS, AL=AL,)

            if decode_status:
                # success
                pg: Pangenome = previous.data_dump('pangenome')
                tree = previous.data_dump('tree')
                previous_file_dict = previous.data_dump('file_dict')

                if previous_file_dict != file_dict:
                    logger.warning(
                        f'File structure has changed')
                    total_name = set(list(file_dict.keys()) +
                                     list(previous_file_dict.keys()))
                    max_width = max(
                        [len(name) for name in total_name.union(set(['previous', 'current']))])+2
                    logger.warning(
                        f'{"previous":<{max_width}}\t{"current":<{max_width}}')
                    new_add = 0
                    loaded_count = 0
                    for strain in total_name:
                        cur_name = strain if strain in file_dict else ''
                        pre_name = strain if strain in previous_file_dict else ''
                        if cur_name != pre_name:  # empty
                            logger.warning(
                                f'{pre_name:<{max_width}}\t{cur_name:<{max_width}}')
                        if cur_name and not pre_name:
                            new_add += 1
                        else:
                            loaded_count += 1
                    if new_add:
                        logger.warning(
                            f'Total {new_add} new strain added. Make sure the preprocess.pkl I loaded is the right one!!!')
                        logger.warning(
                            f'I will reload the file structure from the current input: {indir}')
                        decode_status = False
                    if loaded_count < 2:
                        logger.error(
                            f'Loaded file has less than 2 strains, it is not a valid file that may cause the --exclude_outlier parameter in preprocess step filtered much strains')
                        logger.error(
                            'Please check the input file quality and rerun the preprocess step or just begin from the partition step')
                        raise ValueError('Invalid preprocess.pkl file')

                    if decode_status:
                        file_dict = previous_file_dict
                else:
                    logger.info(
                        f'Load previous file structure from {outdir}/pgap2.pkl')
                    total_bad_gene_num = 0
                    for strain in tqdm(pg.strain_dict, unit=' strain', disable=disable, desc=tqdm_.step(1)):
                        bad_num = pg.strain_dict[strain].bad_gene_num
                        if bad_num > 0:
                            total_bad_gene_num += bad_num
                            logger.warning(
                                f'{strain} invalid gene count: {bad_num}')
                    if total_bad_gene_num > 0:
                        logger.info(
                            f'Total invalid gene count: {total_bad_gene_num}. Check it in log file: {outdir}/preprocess.log')

                    for _ in tqdm([dup_id, orth_id], unit=f" clust iteration", disable=disable, desc=tqdm_.step(2)):
                        ...
            else:
                logger.warning(
                    f'Previous file parameters is not match, start partition from the begining')

    if decode_status is False:
        logger.info(f'Load strain from {indir}')
        pg = file_parser(
            indir=indir, outdir=outdir, annot=annot, threads=threads, disable=disable, retrieve=retrieve, falen=falen, gcode=gcode, id_attr_key=id_attr_key, type_filter=type_filter, prefix='preprocess')
        file_prot = f'{outdir}/total.involved_prot.fa'
        file_annot = f'{outdir}/total.involved_annot.tsv'
        pg.load_annot_file(file_annot)
        pg.load_prot_file(file_prot)
        tree = generate_tree(
            input_file=file_prot, orth_list=[dup_id, orth_id], outdir=pg.outdir, coverage=coverage, LD=LD, AS=AS, AL=AL, falen=falen, disable=disable, threads=threads, evalue=evalue, aligner=aligner, clust_method=clust_method)

    sp = Species(marker_file=marker_file,
                 strain_dict=pg.strain_dict, ani=ani_thre, outdir=outdir)
    logger.info(f'Extract the feature of each strain ...')
    strain_core_stat, strain_single_stat = stat_core(tree, pg)
    genome_attr = {}
    logger.info(f'Extract the genome attributions of each strain ...')
    for strain_index in tqdm(pg.strain_dict, desc=tqdm_.step(4), unit=' strain', disable=disable):
        genome_record = pg.strain_dict[strain_index].genome
        chr_count, genome_len, atcg_dict = genome_feature_stat(
            genome_record)
        half_core = strain_core_stat[strain_index]
        single_cloud = strain_single_stat[strain_index]
        stat_dict = {'chr_count': chr_count, 'genome_len': genome_len,
                     'total_gene_num': sum(pg.strain_dict[strain_index].gene_num),
                     'value_incomplete': pg.strain_dict[strain_index].bad_gene_num,
                     'content': atcg_dict, 'half_core': half_core, 'single_cloud': single_cloud}
        genome_attr.update({strain_index: stat_dict})

    if sp.has_darb():
        darb_strain = sp.get_darb()
        logger.info(f'Representative strain {darb_strain} has been assigned')
    else:
        logger.info(f'Representative strain selecting ...')
        # darb is strain have the highest half core number
        darb = max(genome_attr, key=lambda x: genome_attr[x]['half_core'])
        sp.load_darb(darb)
        darb_strain_name = sp.strain_dict[darb].strain_name
        logger.info(f'Representative strain {darb_strain_name} selected')

    logger.info(f'Gene code and gene length statistics ...')
    sp.stat_gene_code(pg)

    if sp.has_outgroup():
        logger.info(
            f'Coreect the ani according to the outgroup strain: {sp.get_outgroup()}')
        outgroup = sp.get_outgroup()
        this_ani = float(ani_thre)
        this_outgroup = 'Default'
        for each_outgroup in outgroup:
            query = list(SeqIO.parse(each_outgroup, "fasta"))
            query = (bytes(record.seq) for record in query)
            hits = sp.mapper.query_draft(query)
            if hits:
                if hits[0].identity > this_ani:
                    this_ani = hits[0].ani_thre
        if this_ani == ani_thre:
            logger.warning(
                f'Still use default ANI threshold: {ani_thre} because the outgroup strain have no effect on it')
        else:
            logger.info(
                f'ANI threshold is corrected to {this_ani} according to the outgroup strain {this_outgroup}')
            sp.expect_ani = this_ani

    logger.info(f'ANI calculating that will take some time...')
    sp.find_outlier(threads=threads)

    gene_outlier = find_outlier_from_cloud(genome_attr)
    sp.load_genmoe_attr(genome_attr)
    sp.load_gene_outlier(gene_outlier, 'single_gene')

    outlier_dict = sp.get_outlier()

    logger.info(f'preprocess results output...')
    preprocess_output(outlier_dict, outdir, sp, file_dict)

    new_file_dict = get_file_dict(indir=f'{outdir}/refined_input')
    # reload the input file, due to unconsist number of input
    if len(new_file_dict) != len(file_dict):
        discard_strains = list(set(list(file_dict.keys())).difference(
            set(list(new_file_dict.keys()))))
        logger.warning('Total {} strain were discarded due to low quality or dissimilirity with others.'.format(
            len(discard_strains)))
        for discard_strain in discard_strains:
            logger.warning(f'discarded strain: {discard_strain}')
        indir = f'{outdir}/refined_input'
        file_dict = new_file_dict
    else:  # If it is same number after preprocess, it will not rerun the file_parser
        file_dict = file_dict

    pickle_preprocess = PklCheck(outdir=outdir, name='preprocess')
    pickle_preprocess.load('file_dict', main_data=file_dict)
    pickle_preprocess.load('pangenome', main_data=pg, parameter={'orth_id': orth_id, 'para_id': para_id, 'dup_id': dup_id, 'accurate': accurate,
                                                                 'id_attr_key': id_attr_key, 'type_filter': type_filter,
                                                                 'coverage': coverage, 'AS': AS, 'AL': AL, 'LD': LD, 'retrieve': retrieve,
                                                                 'evalue': evalue, 'aligner': aligner, 'clust_method': clust_method,
                                                                 'annot': annot, 'falen': falen})
    pickle_preprocess.load('tree', main_data=tree)
    pickle_preprocess.pickle_()

    if nodraw:
        logger.info('Drawing is disabled')
    else:
        logger.info(f'preprocess results drawing...')
        html_report = preprocess_draw(outlier_dict, outdir,
                                      sp, genome_attr)

        logger.info(f'Report at {html_report}')
        preprocess_draw_vector(sfw=sfw.draw_prep,
                               outdir=outdir, single_file=single_file,
                               ani_threshold=sp.expect_ani
                               )
        logger.info(f'Vector report at {outdir}:')
        if single_file:
            logger.info(
                f'[1/5] postprocess.ANI.pdf')
            logger.info(
                f'[2/5] postprocess.gene_number.pdf')
            logger.info(
                f'[3/5] postprocess.half_core.pdf')
            logger.info(
                f'[4/5] postprocess.proportion.pdf')
            logger.info(
                f'[5/5] postprocess.gene_code.pdf')
        else:
            logger.info(
                f'[1/1] pgap2.preprocess.pdf')
        logger.info('Done')

    return 0


def launch(args: argparse.Namespace):
    indir = os.path.abspath(args.indir)
    outdir = os.path.abspath(args.outdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    main(indir=indir, outdir=outdir,
         orth_id=args.orth_id, para_id=args.para_id, dup_id=args.dup_id, accurate=args.accurate,
         id_attr_key=args.id_attr_key, type_filter=args.type_filter,
         LD=args.LD, AS=args.AS, AL=args.AL,
         evalue=args.evalue,
         aligner=args.aligner, clust_method=args.clust_method,
         coverage=0.98,
         #  coverage=args.coverage,
         nodraw=args.nodraw, single_file=args.single_file,
         marker_file=args.marker, ani_thre=args.ani_thre,
         annot=args.annot, threads=args.threads, gcode=args.gcode,
         disable=args.disable, retrieve=args.retrieve, falen=args.min_falen,)
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
                                      default='cdhit', choices=('cdhit', 'mmseqs2'), help='The method used to cluster the genes.')
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
        '--threads', '-t', required=False, default=8, help='threads used in parallel', type=int)
    subparser_preprocess.set_defaults(func=preprocess_portal)
