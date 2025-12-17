import os
import re
import gzip
import shutil
import zipfile
import tempfile

import numpy as np

from tqdm import tqdm
from typing import IO
from loguru import logger
from functools import partial
from multiprocessing import get_context
from collections import OrderedDict, defaultdict

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, SimpleLocation

from pgap2.lib.pangenome import Pangenome
from pgap2.lib.strain import Strain
from pgap2.utils.supply import run_command
from pgap2.utils.supply import sfw, tqdm_

"""

Data loader module for PGAP2.
This module provides functions to load genomic data from various file formats,
including FASTA, GenBank, and GFF3 output.
It includes functions for parsing these files, extracting gene sequences,
and loading annotations into the PGAP2 data structures.

input: 
- indir: Input directory containing genomic files.
- outdir: Output directory for processed data.
- annot: Boolean indicating whether to generate annotation files.
- threads: Number of threads to use for parallel processing.
- disable: Boolean to disable progress bars.
- id_attr_key: Attribute key to use for gene IDs.
- type_filter: Type of features to filter (e.g., 'CDS').
- retrieve: Boolean indicating whether to retrieve sequences.
- falen: Minimum length of gene sequences to consider valid.
- gcode: Genetic code to use for translation.
- prefix: Prefix for output files.

output:
- processed: Boolean indicating whether the data was processed successfully.
- stats: Dictionary containing statistics about the processed data.
- pangenome: Pangenome object containing the loaded strains and genes.

"""


def set_logger(logger_):
    global logger
    logger = logger_


def extract_attr_from_prodigal(file):
    total_gene = 0
    complete_gene = 0
    total_score = []
    gene_dict = {}
    # stat_tmp = {'total_gene_num': 0, 'score': 0.0,
    #             'value_complete': 0, 'value_incomplete': 0, 'score_list': []}
    this_seq_id = 'ERROR_WHEN_SAW_ME'
    this_seq_num = 'ERROR_WHEN_SAW_ME'
    gene_num = 0
    with open(file) as fh:
        for line in fh:
            line = line.strip('\n')
            if re.match(r'DEFINITION\s+.*seqhdr=', line):
                try:
                    this_seq_id = re.findall(
                        r'DEFINITION.*seqhdr="(.*)";version', line)[0]
                    this_seq_num = re.findall(
                        r'DEFINITION.*seqnum=(\d+);seqlen', line)[0]
                    gene_num = 0
                except:
                    logger.error(
                        f'Cannot find the correct header from prodigal result {file}')
            elif re.match(r'\s+CDS\s+(\S+)\s*$', line):
                gene_num += 1
                strand, start, end = '', 0, 0
                strand = -1 if re.match(r'complement', line) else 1
                start = re.findall(r'CDS\s+.*?(\d+).*?(\d+)', line)[0][0]
                end = re.findall(r'CDS\s+.*?(\d+).*?(\d+)', line)[0][1]
                geneid = f'{this_seq_num}_{gene_num}'
                # biopython FeatureLocation is 0-based coordination.
                start = int(start)-1
                per = SeqFeature(SimpleLocation(int(start), int(end), strand=strand),
                                 type="CDS", id=geneid)
                if this_seq_id not in gene_dict:
                    gene_dict[this_seq_id] = []
                gene_dict[this_seq_id].append(per)

            elif re.match(r'\s+\/note=', line):
                total_gene += 1
                partial = re.match(
                    r'.*partial=(\d+).*score=(.*?);cscore', line.strip()).group(1)
                score = re.match(
                    r'.*partial=(\d+).*score=(.*?);cscore', line.strip()).group(2)
                total_score.append(float(score))
                if partial == '00':
                    complete_gene += 1

    # most_score = np.median(total_score)
    # portion_complete = complete_gene/total_gene
    # co = round(portion_complete*most_score, 3)
    # logger.debug(
    #     f"{file} co:{co} portion_complete:{portion_complete} median:{most_score} ")

    logger.debug(f'Reading the annotation result of {file}...')
    logger.debug(f'#Total gene: {total_gene}')
    logger.debug(f'#Complete gene: {complete_gene}')
    # logger.debug(f'coefficient: {co}')
    # stat_tmp['score'] = co
    # stat_tmp['total_gene_num'] = total_gene
    # stat_tmp['value_complete'] = complete_gene
    # stat_tmp['value_incomplete'] = total_gene-complete_gene
    # stat_tmp['score_list'] = total_score
    return gene_dict


def fa_parser(genome_file, strain_name, temp_out, strain_index: int, annot: bool, partial: bool = False, retrieve: bool = False, falen: int = 10, gcode: int = 11):
    logger.debug(f'Reding genome file: {genome_file}')
    in_seq_file = genome_file
    in_seq_handle = open(in_seq_file)
    seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
    in_seq_handle.close()

    run_command(
        f"{sfw.prodigal} -i {genome_file} -o {temp_out}/{strain_name}.prodigal")
    gene_dict = extract_attr_from_prodigal(
        f'{temp_out}/{strain_name}.prodigal',)

    prot_file = f'{temp_out}/{strain_name}.prot'
    annot_file = f'{temp_out}/{strain_name}.annot'

    bad_gene_num = 0
    good_gene_num = []
    contig_name_out = None
    contig_name_map = {}
    with open(annot_file, 'w') as annot_fh, open(prot_file, 'w') as prot_fh:
        for contig_name in gene_dict:
            if contig_name_out is None:
                contig_name_out = contig_name
                contig_index = 0
                gene_index = 0
                good_gene_num.append(0)
            else:
                if contig_name != contig_name_out:
                    contig_index += 1
                    contig_name_out = contig_name
                    gene_index = 0
                    good_gene_num.append(0)
            contig_name_map[contig_name] = contig_index
            for feature in gene_dict[contig_name]:
                if re.match(r'pgap2_\S+_\d+', feature.id):
                    logger.warning(
                        f"Cannot find right Gene ID in {feature.id} through feature 'Parent=XXX', replace it with CDS id {feature.id}, that may cause some problem when that gene have multiple exons.")

                id_name = feature.id
                nucl_fa = feature.extract(seq_dict[contig_name].seq)
                gene_name = feature.qualifiers.get('gene', '')
                product_name = feature.qualifiers.get('product', '')
                try:
                    prot_fa = nucl_fa.translate(table=gcode, cds=True)
                    gene_len = len(prot_fa)
                    if gene_len < falen:
                        raise ValueError(
                            f'Gene {id_name} translated length is too short. Length: {gene_len} < {falen}')
                    location = feature.location
                    gene_name_index = f'{strain_index}:{contig_index}:{gene_index}'

                    prot_fh.write(f'>{gene_name_index}\n{prot_fa}\n')
                    annot_fh.write(
                        f'{gene_name_index}\t{strain_name}\t{contig_name}\t{location}\t{gene_len}\t{id_name}\t{gene_name}\t{product_name}\t{nucl_fa}\t{prot_fa}\n')
                    gene_index += 1
                    good_gene_num[-1] += 1
                except Exception as e:
                    bad_gene_num += 1
                    logger.debug(
                        f'[Skip unregular gene] {id_name} from {strain_name}: {e}')

    if os.path.exists(f'{temp_out}/../../genome_index/'):  # retrieve mode
        dir_index = strain_index // 1000
        strain_index_path = f'{temp_out}/../../genome_index/{dir_index}/{strain_index}'
        os.makedirs(strain_index_path, exist_ok=True)
        records = SeqIO.parse(genome_file, 'fasta')
        for record in records:
            record.id = str(contig_name_map[record.id])
        genome_file = f'{strain_index_path}/ref.fa'
        SeqIO.write(records, genome_file, 'fasta')

    if retrieve:
        run_command(
            f"{sfw.miniprot} -t 1 -d {strain_index_path}/ref.mpi {genome_file}")

    return good_gene_num, bad_gene_num, annot_file, prot_file


def gbf_parser(gbf_file, strain_name, temp_out, strain_index: int, annot: bool, partial: bool = False, retrieve: bool = False, falen: int = 11, gcode: int = 11, read_type: str = 'CDS', read_attr: str = 'gene'):
    logger.debug(f'Reading genome file: {gbf_file}')

    if annot:
        output_file = f'{temp_out}/{strain_name}.genome.fa'
        records = SeqIO.parse(gbf_file, 'genbank')
        SeqIO.write(records, output_file, 'fasta')
        good_gene_num, bad_gene_num, annot_file, prot_fh = fa_parser(genome_file=output_file, strain_name=strain_name,
                                                                     temp_out=temp_out, strain_index=strain_index, annot=True, partial=partial, falen=falen)
    else:
        in_seq_handle = open(gbf_file)
        seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "genbank"))

        annot_file = f'{temp_out}/{strain_name}.annot'
        prot_file = f'{temp_out}/{strain_name}.prot'

        bad_gene_num = 0
        good_gene_num = []
        contig_name_out = None
        contig_name_map = {}
        with open(annot_file, 'w') as annot_fh, open(prot_file, 'w') as prot_fh:
            for contig_name in seq_dict:
                if contig_name_out is None:
                    contig_name_out = contig_name
                    contig_index = 0
                    gene_index = 0
                    good_gene_num.append(0)
                else:
                    if contig_name != contig_name_out:
                        contig_index += 1
                        contig_name_out = contig_name
                        gene_index = 0
                        good_gene_num.append(0)
                contig_name_map[contig_name] = contig_index
                rec = seq_dict[contig_name]
                gene_order_map = {}
                for per in rec.features:
                    if per.type == read_type:
                        start = per.location.start
                        if start in gene_order_map:
                            logger.debug(
                                f"Duplicate gene start position [{start}] in [{contig_name}].")
                            start += 1
                        gene_order_map[start] = per
                for per_pos in sorted(gene_order_map.keys()):
                    per = gene_order_map[per_pos]
                    gene_name = per.qualifiers.get('gene', '')
                    product_name = per.qualifiers.get('product', '')
                    id_name = per.qualifiers.get(read_attr, '')
                    locus_tag = per.qualifiers.get('locus_tag', '')
                    if not id_name:
                        if locus_tag:
                            id_name = locus_tag[0]
                        elif gene_name:
                            id_name = gene_name
                        else:
                            id_name = per.location
                    nucl_fa = per.extract(rec.seq)
                    try:
                        prot_fa = nucl_fa.translate(table=gcode, cds=True)
                        gene_len = len(prot_fa)
                        if gene_len < falen:
                            raise ValueError(
                                f'Gene {id_name} translated length is too short. Length: {gene_len} < {falen}')
                        location = per.location
                        gene_name_index = f'{strain_index}:{contig_index}:{gene_index}'
                        prot_fh.write(f'>{gene_name_index}\n{prot_fa}\n')
                        annot_fh.write(
                            f'{gene_name_index}\t{strain_name}\t{contig_name}\t{location}\t{gene_len}\t{id_name}\t{gene_name}\t{product_name}\t{nucl_fa}\t{prot_fa}\n')
                        gene_index += 1
                        good_gene_num[-1] += 1
                    except Exception as e:
                        bad_gene_num += 1
                        logger.debug(
                            f'[Skip unregular gene] {id_name} from {strain_name}: {e}')
        in_seq_handle.close()

        if os.path.exists(f'{temp_out}/../../genome_index/'):
            dir_index = strain_index // 1000
            strain_index_path = f'{temp_out}/../../genome_index/{dir_index}/{strain_index}'
            os.makedirs(strain_index_path, exist_ok=True)
            records = list(SeqIO.parse(gbf_file, 'genbank'))
            for record in records:
                record.id = str(contig_name_map[record.id])
            genome_file = f'{strain_index_path}/ref.fa'
            SeqIO.write(records, genome_file, 'fasta')

        if retrieve:
            run_command(
                f"{sfw.miniprot} -t 1 -d {strain_index_path}/ref.mpi {genome_file}")

    return good_gene_num, bad_gene_num, annot_file, prot_file


def open_file(file_name) -> IO:
    compression_formats = {".gz": gzip.open, ".zip": zipfile.ZipFile}

    file_format = None
    for format, opener in compression_formats.items():
        if file_name.endswith(format):
            file_format = format
            break

    if file_format:
        if file_format == ".gz":
            return opener(file_name, 'rt')
        elif file_format == ".zip":
            with opener(file_name, 'r') as zip_file:
                # assuming the zip file contains only one file
                return zip_file.open(zip_file.namelist()[0], 'r')
    else:
        # if no compression, just open the file normally
        return open(file_name, 'r')


def dict_to_fasta(dictionary, name_map, output_file):
    records = []
    for key, value in dictionary.items():
        if name_map:
            value.id = str(name_map[key])
        records.append(value)

    with open(output_file, 'w') as handle:
        SeqIO.write(records, handle, 'fasta')

    return output_file


def gffa_parser(gffa_file, fa_file, strain_name, temp_out, strain_index: int, annot: bool, partial: bool = False, retrieve: bool = False, falen: int = 10, gcode: int = 11, read_type: str = 'CDS', read_attr: str = 'ID'):
    logger.debug(f'Reading gffa file: {gffa_file}')
    in_file = gffa_file

    if fa_file and os.path.exists(fa_file):
        seq_dict = SeqIO.to_dict(SeqIO.parse(fa_file, "fasta"))
    else:
        seq_dict = SeqIO.to_dict(SeqIO.parse(gffa_file, "fasta"))

    if not seq_dict:
        logger.error(
            f'Cannot read the genome file {fa_file} or any sequence in {gffa_file}')
        raise ValueError(
            f'Cannot read the genome file {fa_file} or any sequence in {gffa_file}')

    gene_dict = defaultdict(dict)
    # ATTR_REGEX = re.compile(r'\W')

    def parse_attributes(attributes):
        """
        解析 GFF 第九列的属性字段，返回一个字典。
        """
        attr_dict = {}
        for attr in attributes.split(';'):
            key_value = attr.split('=', 1)
            if len(key_value) == 2:
                key, value = key_value
                attr_dict[key.lower()] = re.sub(
                    ',', '_', value)  # replace illegal characters with underscores
        return attr_dict

    with open_file(in_file) as fh:
        count = 0
        for line in fh:
            if line.startswith('#'):
                continue
            line_list = line.strip('\n').split('\t')
            if (len(line_list) == 9) and (line_list[2].lower() == read_type.lower()):
                count += 1
                contig_name = line_list[0]
                location = SimpleLocation(
                    int(line_list[3])-1,  # 0-based start
                    int(line_list[4]),  # [)
                    strand=1 if str(line_list[6]) in ('+', '.') else -1)

                attrs = parse_attributes(line_list[8])
                geneid = attrs.get(read_attr.lower(), None)
                parent = attrs.get('parent', geneid)
                name = attrs.get('name', '')
                gene = attrs.get('gene', '')
                product = attrs.get('product', '')
                if geneid is None:
                    logger.warning(
                        f'Cannot find {read_attr} in {line_list[8]}')

                # if parent is None:
                if read_attr.lower() != 'cds':  # for test simpan
                    parent = geneid

                if parent in gene_dict[contig_name]:
                    logger.debug(f"Duplicate parent ID [{parent}] in contig [{contig_name}], "
                                 f"seems to be a gene with multiple exons.")
                    gene_dict[contig_name][parent].location += location
                else:
                    # add new gene feature
                    feature = SeqFeature(
                        location=location,
                        type="CDS",
                        id=geneid,
                        qualifiers={'product': product,
                                    'gene': gene, 'name': name}
                    )
                    gene_dict[contig_name][parent] = feature

    keys_to_keep = set(seq_dict.keys()) & set(gene_dict.keys())
    seq_dict = {key: seq_dict[key] for key in keys_to_keep}

    if annot:
        genome_fa = dict_to_fasta(
            seq_dict, None, f'{temp_out}/{strain_name}.genome.fa')
        good_gene_num, bad_gene_num, annot_file, prot_file = fa_parser(genome_file=genome_fa, strain_name=strain_name, temp_out=temp_out,
                                                                       strain_index=strain_index, annot=True, partial=partial, falen=falen)
    else:
        annot_file = f'{temp_out}/{strain_name}.annot'
        prot_file = f'{temp_out}/{strain_name}.prot'
        bad_gene_num = 0
        good_gene_num = []
        contig_name_out = None
        contig_name_map = {}
        with open(annot_file, 'w') as annot_fh, open(prot_file, 'w') as prot_fh:
            for contig_name in gene_dict:
                if contig_name_out is None:
                    contig_name_out = contig_name
                    contig_index = 0
                    gene_index = 0
                    good_gene_num.append(0)
                else:
                    if contig_name != contig_name_out:
                        contig_index += 1
                        contig_name_out = contig_name
                        gene_index = 0
                        good_gene_num.append(0)
                contig_name_map[contig_name] = contig_index
                gene_order_map = {}
                used_starts = set()
                for gene, feature in gene_dict[contig_name].items():
                    start = gene_dict[contig_name][gene].location.start
                    while start in used_starts:
                        logger.debug(
                            f"Duplicate gene start position [{start+1}] in contig [{contig_name}].")
                        start += 1
                    used_starts.add(start)
                    feature.id = gene
                    gene_order_map[start] = feature
                for start in sorted(gene_order_map.keys()):
                    feature = gene_order_map[start]
                    id_name = feature.id
                    nucl_fa = feature.extract(seq_dict[contig_name].seq)
                    gene_name = feature.qualifiers.get('gene', '')
                    product_name = feature.qualifiers.get('product', '')
                    try:
                        prot_fa = nucl_fa.translate(table=gcode, cds=True)
                        gene_len = len(prot_fa)
                        if gene_len < falen:
                            raise ValueError(
                                f'Gene {id_name} translated length is too short. Length: {gene_len} < {falen}')
                        location = feature.location
                        gene_name_index = f'{strain_index}:{contig_index}:{gene_index}'
                        prot_fh.write(
                            f'>{gene_name_index} {location}\n{prot_fa}*\n')
                        annot_fh.write(
                            f'{gene_name_index}\t{strain_name}\t{contig_name}\t{location}\t{gene_len}\t{id_name}\t{gene_name}\t{product_name}\t{nucl_fa}\t{prot_fa}\n')
                        gene_index += 1
                        good_gene_num[-1] += 1
                    except Exception as e:
                        bad_gene_num += 1
                        logger.debug(
                            f'[Skip unregular gene] {id_name} from {strain_name}: {e}')

        if os.path.exists(f'{temp_out}/../../genome_index/'):
            dir_index = strain_index // 1000
            strain_index_path = f'{temp_out}/../../genome_index/{dir_index}/{strain_index}'
            os.makedirs(strain_index_path, exist_ok=True)
            genome_file = dict_to_fasta(
                seq_dict, contig_name_map, f'{strain_index_path}/ref.fa')

        if retrieve:
            run_command(
                f"{sfw.miniprot} -t 1 -d {strain_index_path}/ref.mpi {genome_file}")
    return good_gene_num, bad_gene_num, annot_file, prot_file


def pool_file_parser(file_dict_with_index, falen, retrieve, annot, temp_out, gcode, id_attr_key, type_filter):
    strain_index, file_list = file_dict_with_index
    strain_name, file_dict = file_list
    dir_index = strain_index//1000
    temp_out = f'{temp_out}/{dir_index}'
    os.makedirs(temp_out, exist_ok=True)

    if 'gbf' in file_dict:
        gbf_file = file_dict['gbf']
        good_gene_num, bad_gene_num, annot_file, prot_file = gbf_parser(gbf_file, strain_name,
                                                                        temp_out, strain_index, annot, retrieve=retrieve, falen=falen, gcode=gcode, read_type=type_filter, read_attr=id_attr_key)

    elif 'gff' in file_dict:
        gffa_file = file_dict['gff']
        fa_file = file_dict.get('fa', None)
        good_gene_num, bad_gene_num, annot_file, prot_file = gffa_parser(gffa_file, fa_file, strain_name,
                                                                         temp_out, strain_index, annot, retrieve=retrieve, falen=falen, gcode=gcode, read_type=type_filter, read_attr=id_attr_key)

    elif ('fa' in file_dict) and annot:
        good_gene_num, bad_gene_num, annot_file, prot_file = fa_parser(gffa_file, strain_name,
                                                                       temp_out, strain_index, annot, retrieve=retrieve, falen=falen, gcode=gcode)
    else:
        logger.warning(
            f'Skip because [{strain_name}] is not recognized.')
        logger.error(
            f'PGAP2 can only understand gff, fa, gbf file')
        return None, None, None, None, None, None

    return good_gene_num, bad_gene_num, strain_name, strain_index, annot_file, prot_file


def process_file(file_pair):
    """
    read protein and annotation files from a pair of file paths.
    This function is used in multiprocessing to read files concurrently.
    :param file_pair: A tuple containing the annotation file path and protein file path.
    :return: A tuple containing the protein lines and annotation lines.
    """
    annot_file, prot_file = file_pair
    prot_lines = []
    annot_lines = []

    # read protein sequences
    with open(prot_file) as pf:
        prot_lines = pf.readlines()

    # read annotation information
    with open(annot_file) as af:
        annot_lines = af.readlines()

    return prot_lines, annot_lines


def file_parser(indir, outdir, annot, threads: int,  disable: bool = False, id_attr_key: str = 'ID', type_filter: str = 'CDS', retrieve: bool = False, falen: int = 11, gcode: int = 11, prefix='partition') -> Pangenome:
    temp_out = tempfile.mkdtemp(dir=outdir)
    if retrieve or prefix == 'preprocess':
        genome_index_path = f'{outdir}/genome_index'
        if not os.path.exists(genome_index_path):
            try:
                os.makedirs(genome_index_path)
            except OSError as e:
                logger.info(
                    f"Error creating folder '{genome_index_path}': {e}")
        else:
            try:
                shutil.rmtree(genome_index_path)
                os.makedirs(genome_index_path)
            except OSError as e:
                logger.info(
                    f"Error creating folder '{genome_index_path}': {e}")

    file_dict = get_file_dict(indir)

    pg = Pangenome(outdir=outdir, threads=threads,
                   gcode=gcode, disable=disable)
    logger.debug(f'Sequence extraction in {temp_out}')

    bar = tqdm(range(len(file_dict)),
               unit=" strain", disable=disable, desc=tqdm_.step(1))

    total_gene_num = 0
    per_file_list = []

    with get_context('fork').Pool(processes=threads, initializer=set_logger, initargs=(logger,)) as p:
        total_bad_gene_num = 0
        for good_gene_num, bad_gene_num, strain_name, strain_index, annot_file, prot_file in p.imap_unordered(partial(pool_file_parser, falen=falen, retrieve=retrieve, annot=annot, temp_out=temp_out, gcode=gcode, id_attr_key=id_attr_key, type_filter=type_filter), enumerate(file_dict.items())):
            if bad_gene_num is None:
                continue
            if bad_gene_num > 0:
                logger.warning(
                    f'{strain_name} invalid gene count: {bad_gene_num}')
            if sum(good_gene_num) == 0:
                logger.warning(
                    f'{strain_name} has no valid gene. Try using -r to retrieve it.')
            per_file_list.append((annot_file, prot_file))
            total_bad_gene_num += bad_gene_num
            strain = Strain(strain_name=strain_name,
                            strain_index=strain_index,
                            bed_gene_num=bad_gene_num,
                            gene_num=good_gene_num)
            total_gene_num += sum(good_gene_num)
            pg.load_strain(
                strain=strain)
            bar.update()
        bar.close()

        logger.info(f'Total loaded gene count: {total_gene_num}')
        pg.total_gene_num = total_gene_num
        if total_bad_gene_num > 0:
            logger.info(
                f'Total invalid gene count: {total_bad_gene_num}')
            logger.info(f'Check all record in log file: {outdir}/{prefix}.log')

    logger.info(f'Writing the total involved gene and annotation...')
    batch_size = 1000
    for i in range(0, len(per_file_list), batch_size):
        batch = per_file_list[i:i + batch_size]

        with open(f'{outdir}/total.involved_prot.fa', 'a') as prot_fh, open(f'{outdir}/total.involved_annot.tsv', 'a') as annot_fh:
            if i == 0:  # Write header only once
                annot_fh.write(
                    f'#Gene_index\tStrain\tContig\tLocation\tLength\tGene_ID\tGene_name\tProduct_name\tNucleotide_sequence\tProtein_sequence\n')

            for annot_file, prot_file in tqdm(batch, unit=' strain', desc=tqdm_.step(1), disable=disable):
                with open(prot_file) as fh:
                    # prot_lines.extend(fh.readlines())
                    for line in fh:
                        prot_fh.write(line)
                with open(annot_file) as fh:
                    # annot_lines.extend(fh.readlines())
                    for line in fh:
                        annot_fh.write(line)

    logger.info(
        f'Check the total involved protein sequence in {outdir}/total.involved_prot.fa')
    logger.info(
        f'Check the total annotation in {outdir}/total.involved_annot.tsv')
    # shutil.rmtree(temp_out)  # clean up the temporary directory

    return pg


def get_file_dict(indir):
    support_format = {'gffa': 'gff', 'gff': 'gff',  'gff3': 'gff',
                      'fasta': 'fa', 'fa': 'fa', 'fsa': 'fa', 'fna': 'fa',
                      'gbk': 'gbf', 'gbff': 'gbf', 'gbf': 'gbf'}

    file_dict = {}

    for filename in os.listdir(indir):
        filename = os.path.abspath(f'{indir}/{filename}')
        basename = os.path.basename(filename)

        while basename.endswith(('.gz', '.zip')):
            basename = os.path.splitext(basename)[0]

        suffix = os.path.splitext(basename)[-1][1:]

        file_format = support_format.get(suffix, None)
        if file_format is None:
            logger.warning(
                f'Skip because [{basename}] is not recognized.')
            continue

        # get the strain name from the file name
        strain_name = os.path.splitext(basename)[0]

        if strain_name not in file_dict:
            file_dict[strain_name] = {file_format: filename}
        else:
            file_dict[strain_name].update({file_format: filename})

    if len(file_dict) == 0:
        logger.error(f'No file found in {indir}')
        logger.error(
            f'PGAP2 can only understand {list(support_format.keys())} file, please check the input directory {indir}')
        exit(1)

    # sort and put it into the OrderedDict
    file_dict = OrderedDict(sorted(file_dict.items(), key=lambda x: x[0]))
    return file_dict
