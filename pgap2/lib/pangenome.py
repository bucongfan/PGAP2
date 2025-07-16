
import os
import re
import networkx as nx

from tqdm import tqdm
from loguru import logger
from collections import defaultdict, OrderedDict

from pgap2.utils.supply import tqdm_

"""
This module defines the core class of PGAP2: Pangenome.

The Pangenome class serves as the top-level structure in PGAP2’s analysis pipeline, sitting above 
Species and Strain in the hierarchy. It records all global information and configuration parameters 
used throughout the pan-genome analysis process.

This includes attributes such as the output directory, number of threads, total number of genomes, 
and references to genome-level data structures (e.g., genome dictionaries). In addition to serving 
as a central data container, the class also supports functionality for data export, loading, 
and checkpointing of intermediate results.

"""


def pan_judger(ortho_num: int, total_num: int):
    freq = float(ortho_num/total_num)
    tag = 'ERROR_WHEN_JUDGE_PAN'

    if freq == 1:
        tag = 'Strict_core'
    elif freq >= 0.99:
        tag = 'Core'
    elif freq >= 0.95:
        tag = 'Soft_core'
    elif freq >= 0.15:
        tag = 'Shell'
    elif freq > 0:
        tag = 'Cloud'
    else:
        raise ValueError(logger.error(
            f'Cannot find right group of frequency {freq} which ortho_num is {ortho_num} and total_num is {total_num}. Contact me at github please.'))
    return tag


class Pangenome():
    def __init__(self, outdir, threads, gcode, disable) -> None:

        self.outdir = outdir
        self.disable_tqdm = disable
        self.threads: int = threads
        self.gcode = gcode
        self.orth_id: float = 0
        self.para_id: float = 0
        self.dup_id: float = 0
        self.accurate: bool = False
        self.exhaust_orth: bool = False
        self.retrieve: bool = False
        self.evalue: float = 1e-5
        self.aligner: str = 'diamond'
        self.LD: int = 0
        self.AL: int = 0
        self.AS: int = 0
        self.annot_file = None
        self.total_gene_num = 0

        self.strain_dict = {}
        self.pan_clust = ''
        self.pan_array = []
        self.pan_array_symbol = []
        self.pav_array = []
        self.pan_attr = []
        self.not_outlier = []
        self.strain_index = {}

    @property
    def strain_num(self):
        return len(self.strain_dict)

    @property
    def strain_name(self):
        return [_ for _ in self.strain_dict]

    def init_pan_temp(self):
        self.one_pan = [""] * self.strain_num
        self.one_pan_symbol = [""] * self.strain_num
        self.one_pav = [0] * self.strain_num

    def outlier(self, file_dict):
        self.not_outlier = list(file_dict.keys())

    def load_strain(self, strain):
        self.strain_dict[strain.strain_index] = strain

    def load_annot_file(self, annot_file):
        if os.path.exists(annot_file):
            self.annot_file = annot_file
        else:
            logger.error(
                f'Cannot find the annot file {annot_file}. Maybe you should check the path and run the program again')

    def load_prot_file(self, prot_file):
        if os.path.exists(prot_file):
            self.prot_file = prot_file
        else:
            logger.error(
                f'Cannot find the annot file {prot_file}. Maybe you should check the path and run the program again')

    def _get_strain_name_list(self):
        strain_name_list = [str() for _ in range(self.strain_num)]
        for strain_index, strain in self.strain_dict.items():
            strain_name = strain.strain_name
            strain_name_list[strain_index] = strain_name
        return strain_name_list

    def dump_csv(self, outdir='', prefix='pgap2.partition'):
        logger.info(f'Dump csv matrix to {outdir}/{prefix}.gene_content.csv')

        strain_name_list = self._get_strain_name_list()
        header = '#Clust,'+','.join(strain_name_list)
        header2 = '#Clust\t'+'\t'.join(strain_name_list)
        header3 = '#Clust\tgene_name\tproduct\tgroup\trepre_gene\tmin\tmean\tvar\tuni\tinvolved_strain\tpara_strain\tinvolved_gene\tpara_gene\t' + \
            ','.join(strain_name_list)
        statistic_dict = OrderedDict({'Strict_core': 0, 'Core': 0,
                                     'Soft_core': 0, 'Shell': 0, 'Cloud': 0, 'Total': 0})
        total_strain_num = self.strain_num
        with open(f'{outdir}/{prefix}.gene_content.csv', 'w') as fh, open(f'{outdir}/{prefix}.gene_content.pav', 'w') as fh2, open(f'{outdir}/{prefix}.gene_content.detail.tsv', 'w') as fh3, open(f'{outdir}/{prefix}.summary_statistics.txt', 'w') as fh4:
            fh.write(f'{header}\n')
            fh2.write(f'{header2}\n')
            fh3.write(f'{header3}\n')
            for i, one_pan in enumerate(self.pan_array):
                one_pav = self.pav_array[i]
                one_pan_symbol = self.pan_array_symbol[i]
                para_strain = 0
                para_gene_num = 0
                involved_strain = 0
                involved_gene = 0
                for each_gene in one_pav:
                    if each_gene > 0:
                        involved_strain += 1
                        involved_gene += each_gene
                    if each_gene > 1:
                        para_strain += 1
                        para_gene_num += each_gene-1

                group = pan_judger(
                    ortho_num=involved_strain, total_num=total_strain_num)
                statistic_dict[group] += 1
                statistic_dict['Total'] += 1

                min = self.pan_attr[i]['min']
                uni = self.pan_attr[i]['uni']
                mean = self.pan_attr[i]['mean']
                var = self.pan_attr[i]['var']

                repre_node = self.pan_attr[i]['repre_node']
                gene_name = self.pan_attr[i]['gene_name']
                gene_name = '[]' if not gene_name else gene_name
                gene_product = self.pan_attr[i]['gene_product']
                gene_product = '[]' if not gene_product else gene_product

                fh.write('clust_{},{}\n'.format(
                    i, ','.join(one_pan)))
                fh2.write('clust_{}\t{}\n'.format(
                    i, '\t'.join([str(_) for _ in one_pav])))
                fh3.write('clust_{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    i, gene_name, gene_product, group, repre_node,
                    min, mean, var, uni, involved_strain, para_strain, involved_gene, para_gene_num, ','.join(one_pan_symbol)))

            for k, v in statistic_dict.items():
                freq = round(int(v)*100/int(statistic_dict['Total']), 2)
                limit = ''
                if k == 'Strict_core':
                    limit = '(strains = 100%)'
                    k = 'Strict core genes'
                elif k == 'Core':
                    limit = '(99% <= strains < 100%)'
                    k = 'Core genes'
                elif k == 'Soft_core':
                    limit = '(95% <= strains < 99%)'
                    k = 'Soft core genes'
                elif k == 'Shell':
                    limit = '(15% <= strains < 95%)'
                    k = 'Shell genes'
                elif k == 'Cloud':
                    limit = '(0% <= strains < 15%)'
                    k = 'Cloud genes'
                elif k == 'Total':
                    limit = '(0% <= strains <= 100%)'
                    k = 'Total genes'
                else:
                    logger.error('ERROR: tell me at github')
                fh4.write('{}\t{}\t{}\t{}%\n'.format(k, limit, v, freq))

        return (f'{outdir}/{prefix}.gene_content.csv', f'{outdir}/{prefix}.gene_content.pav', f'{outdir}/{prefix}.summary_statistics.txt')

    def reload_nucl_file(self, tree):
        self.falen = {}
        if self.accurate:
            self.nucl_fa = {}
            para_ready_members = set()

            for nodes in nx.connected_components(tree.distance_graph):
                strains = set()
                has_para = False
                for node in nodes:
                    if tree.orth_identity_tree.nodes[node]['has_para']:
                        para_ready_members.update(
                            tree.orth_identity_tree.nodes[node]['members'])

                    this_strains = tree.orth_identity_tree.nodes[node]['strains']
                    if strains & this_strains:
                        has_para = True
                        break
                    else:
                        strains |= this_strains
                if has_para:
                    for node in nodes:
                        para_ready_members.update(
                            tree.orth_identity_tree.nodes[node]['members'])

            with open(self.annot_file) as fh:
                for line in fh:
                    if line.startswith('#'):  # skip the header and retrieved gene
                        continue
                    lines = line.strip().split('\t')
                    gene_index = lines[0]
                    seq = lines[8]
                    falen = int(lines[4])
                    self.falen.update({gene_index: falen})
                    if gene_index in para_ready_members:
                        self.nucl_fa.update({gene_index: seq})
            logger.info(
                f'Load {len(self.nucl_fa)} nucl sequences for bidirection best check')
        elif self.retrieve:
            self.nucl_fa = {}
            with open(self.annot_file) as fh:
                for line in fh:
                    if line.startswith('#'):
                        continue
                    lines = line.strip().split('\t')
                    gene_index = lines[0]
                    falen = int(lines[4])
                    seq = lines[8]
                    self.nucl_fa.update({gene_index: seq})
                    self.falen.update({gene_index: falen})
        else:
            with open(self.annot_file) as fh:
                for line in fh:
                    if line.startswith('#'):
                        continue
                    lines = line.strip().split('\t')
                    gene_index = lines[0]
                    falen = int(lines[4])
                    self.falen.update({gene_index: falen})

    def reload_annot_file(self, retrieve=False):
        flat_annot_file = self.annot_file
        self.annot = {}
        self.annot_contig_map = {}
        self.gene_rank = defaultdict(lambda: [[], []])  # contig: [plus, minus]
        if os.path.exists(flat_annot_file):
            bar = tqdm(total=self.total_gene_num,
                       unit=' Gene', desc=tqdm_.step(6), disable=self.disable_tqdm)
            loaded_contig_name = set()
            if retrieve:
                coord_pat = re.compile(r'\[(\d+):\d+\]\(([+-])\)')

            with open(flat_annot_file) as fh:
                for line in fh:
                    if line.startswith('#'):
                        continue
                    lines = line.strip().split('\t')
                    bar.update()
                    gene_index = lines[0]
                    contig_name = lines[2]
                    if contig_name not in loaded_contig_name:
                        loaded_contig_name.add(contig_name)
                        self.annot_contig_map.update(
                            {':'.join(gene_index.split(':')[:2]): contig_name})
                    gene_len = int(lines[4])
                    gene_id = lines[5]
                    gene_name = lines[6]
                    gene_product = lines[7]
                    if retrieve:
                        # start = int(re.search(r'\[(\d+):', lines[3]).group(1))
                        # strand = re.search(r'\[(\+|-)\d+:\d+\]', lines[3]).group(1)
                        m = coord_pat.search(lines[3])
                        if m:
                            start = int(m.group(1))   # 起始坐标
                            strand = m.group(2)        # '+' 或 '-'
                        contig = ":".join(lines[0].split(":")[:2])
                        if strand == '+':
                            self.gene_rank[contig][0].append(start)
                        else:
                            self.gene_rank[contig][1].append(start)
                    self.annot.update({gene_index: {'len': gene_len, 'id': gene_id,
                                                    'name': gene_name, 'product': gene_product}})
            bar.close()
        else:
            logger.error(f'Cannot find the annot file {flat_annot_file}')
            raise FileNotFoundError

    # load pangeenome result
    def load_one_pan(self, pan_clust):
        one_pan = self.one_pan[:]
        one_pan_symbol = self.one_pan_symbol[:]
        one_pav = self.one_pav[:]
        one_attr = {}
        gene_featrue_stat = {'name': set(), 'product': set()}
        for gene in pan_clust.gene_clust_list:
            strain_index = int(gene.split(':', maxsplit=1)[0])
            gene_id = self.annot[gene]['id']
            gene_name = self.annot[gene]['name']
            gene_product = self.annot[gene]['product']
            del self.annot[gene]
            gene_featrue_stat['name'].add(gene_name)
            gene_featrue_stat['product'].add(gene_product)

            if not one_pan[strain_index]:
                one_pan[strain_index] = gene_id
                one_pan_symbol[strain_index] = gene
                one_pav[strain_index] = 1
            elif one_pan[strain_index]:
                one_pan[strain_index] += f';{gene_id}'
                one_pan_symbol[strain_index] += f';{gene}'
                one_pav[strain_index] += 1

        gene_name = ';'.join(gene_featrue_stat['name'])
        gene_product = ';'.join(gene_featrue_stat['product'])

        one_attr.update({'min': pan_clust.min})
        one_attr.update({'var': pan_clust.var})
        one_attr.update({'uni': pan_clust.uni})
        one_attr.update({'mean': pan_clust.mean})
        one_attr.update({'repre_node': pan_clust.repre_node})
        one_attr.update({'gene_name': gene_name})
        one_attr.update(
            {'gene_product': gene_product})
        self.pan_array.append(one_pan)
        self.pan_array_symbol.append(one_pan_symbol)
        self.pav_array.append(one_pav)
        self.pan_attr.append(one_attr)

    def get_feature(self, gene, feature, level='gene'):
        strain = self.strain_dict[int(gene.split(':', maxsplit=1)[0])]
        assert strain is not None, f"Cannot find the strain name of {gene}"
        if feature == 'strain':
            return strain
        if feature == 'gene_rank':
            raise NotImplementedError(
                'The gene rank feature is not implemented yet. Please tell me in github.')
        return strain.get_gene_feature(gene, feature)

    def get_flank_gene(self, gene: str, flank: int):
        gene_list = []
        strain_symbol, contig_symbol, gene_symbol = gene.split(':')
        strain = self.get_feature(gene, feature='strain')
        gene_symbol = int(gene_symbol)
        for i in range(flank):
            i += 1
            if gene_symbol+i <= strain.ori_gene_num[int(contig_symbol)]:
                flank_gene = '{}:{}:{}'.format(
                    strain_symbol, contig_symbol, gene_symbol+i)
                gene_list.append(flank_gene)
            else:
                break
        for i in range(flank):
            i += 1
            if int(gene_symbol-i) >= 0:
                flank_gene = '{}:{}:{}'.format(
                    strain_symbol, contig_symbol, gene_symbol-i)
                gene_list.append(flank_gene)
            else:
                break
        return set(gene_list)
