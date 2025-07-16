import os
import pandas as pd
from loguru import logger
from collections import defaultdict

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pgap2.lib.pangenome import Pangenome

"""
This class is used to store and manage basic pangenome information extracted from the main function,
which is then used for subsequent post-processing analyses.
It includes methods for loading used clusters, phylogeny information, and handling PAV files.
It does not contain any analysis methods itself, serving only as a data carrier.
"""


class Basic():
    def __init__(self, pg: Pangenome) -> None:
        self.outdir = pg.outdir
        self.strain_num = pg.strain_num
        self.strain_dict = pg.strain_dict

        self.gcode = pg.gcode
        self.orth_id: float = pg.orth_id
        self.para_id: float = pg.para_id
        self.dup_id: float = pg.dup_id
        if os.path.exists(f'{pg.outdir}/pgap2.partition.gene_content.csv'):
            self.pan_array = f'{pg.outdir}/pgap2.partition.gene_content.csv'
        else:
            raise ValueError(logger.error(
                f'BUG: {pg.outdir}/pgap2.partition.gene_content.csv not found. Tell me in github'))
        if os.path.exists(f'{pg.outdir}/pgap2.partition.gene_content.pav'):
            self.pav_array = f'{pg.outdir}/pgap2.partition.gene_content.pav'
        else:
            raise ValueError(logger.error(
                f'BUG: {pg.outdir}/pgap2.partition.gene_content.pav not found. Tell me in github'))
        self.pan_attr = pg.pan_attr
        self.forbidden_chars = set(" ;:,()'")

    def dumper(self):
        attr = {
            'Project location': self.outdir,
            'Genome count': self.strain_num,
            'Genetic code': self.gcode,
            'Ortholog threshold': self.orth_id,
            'Paralog threshold': self.para_id,
            'Duplication threshold': self.dup_id,
            'Profile file': self.pan_array,
            'PAV file': self.pav_array,
            'Cluster number': len(self.pan_attr),
        }
        max_key_length = max(len(key) for key in attr.keys())

        for key, value in attr.items():
            yield f'{key:{max_key_length}}: {value}'

    @staticmethod
    def get_gene_from_comma(gene_str):
        gene_list = []
        for each_genes in gene_str.split(','):
            if ';' in each_genes:
                gene_list.extend(each_genes.split(';'))
            else:
                gene_list.append(each_genes)
        return gene_list

    def load_used_clusters(self, clusts, file):
        self.customized_cluster = set()
        if clusts and os.path.exists(clusts):
            with open(clusts) as fh:
                for line in fh:
                    if line.startswith('#'):
                        continue
                    clust_name = line.strip()
                    self.customized_cluster.add(clust_name)
            logger.info(
                f'Total {len(self.customized_cluster)} customized clusters were loaded.')
        else:
            with open(file, 'r') as fh:
                for line in fh:
                    if line.startswith('#'):
                        continue
                    else:
                        clust_name = line.strip().split('\t')[0]
                        self.customized_cluster.add(clust_name)

    def load_pav(self, file):
        self.pav = pd.read_csv(file, sep='\t', index_col=0)

    def phylogeny_from_detail_pav(self, file, core_thre, also_pan, para_strategy):
        self.core_thre = core_thre
        self.phylogeny_dict = {}
        core_clust_num = 0
        pan_clust_num = 0
        with open(file, 'r') as fh:
            detail_header = {}
            for line in fh:
                line = line.strip()
                if line.startswith('#'):
                    for i, col_name in enumerate(line.split('\t')):
                        if i == 0:  # remove the first '#'
                            col_name = col_name[1:]
                        if i == 13:  # add the last column
                            ...
                        else:
                            detail_header[col_name] = i
                    continue
                line = line.strip().split('\t')

                cluster = line[detail_header.get('Clust')]
                para_gene = int(line[detail_header.get('para_gene')])
                group = line[detail_header.get('group')]
                if para_gene > 0 and para_strategy == 'drop':
                    continue
                involved_strain = int(
                    line[detail_header.get('involved_strain')])
                if involved_strain == 1:
                    continue
                if cluster not in self.customized_cluster:
                    continue
                if involved_strain/self.strain_num >= core_thre:
                    core_clust_num += 1
                    clust_gene_list = self.get_gene_from_comma(line[13])
                    self.phylogeny_dict.update(
                        {cluster: {'Type': 'Core', 'Gene': clust_gene_list, 'Group': group}})
                elif also_pan:
                    pan_clust_num += 1
                    clust_gene_list = self.get_gene_from_comma(line[13])
                    self.phylogeny_dict.update(
                        {cluster: {'Type': 'Pan', 'Gene': clust_gene_list, 'Group': group}})

        logger.info(
            f'Read total cluster number: {core_clust_num + pan_clust_num} (single clouds will be discarded as it is meaningless for phylogenetic analysis)')
        logger.info(
            f'Core cluster number: {core_clust_num}, Pan cluster number: {pan_clust_num}')

    def phylogeny_from_id(self, file):
        self.used_cluster = defaultdict(list)

        gene_count = 0
        involved_genes = {}
        for cluster, cluster_info in self.phylogeny_dict.items():
            for gene in cluster_info['Gene']:
                involved_genes[gene] = cluster
        with open(file, 'r') as fh:
            id_header = {}
            for line in fh:
                line = line.strip()
                if line.startswith('#'):
                    for i, group in enumerate(line.split('\t')):
                        if i == 0:
                            group = group[1:]
                        id_header[group] = i
                    continue
                line = line.split('\t')
                gene_index = line[id_header.get('Gene_index')]
                if gene_index not in involved_genes:
                    continue
                cluster = involved_genes[gene_index]
                cds = line[id_header.get('Nucleotide_sequence')]
                description = '{}:{}:{}'.format(line[id_header.get(
                    'Strain')], line[id_header.get('Contig')], line[id_header.get('Gene_ID')])
                record = SeqRecord(Seq(cds), id=gene_index,
                                   description=description)
                self.used_cluster[cluster].append(record)
                gene_count += 1
        logger.info(
            f'Get total gene number: {gene_count} from {len(self.used_cluster)} clusters')

    def get_real_strain_name(self, symbol):
        strain_name = self.strain_dict[symbol].strain_name
        clean_strain = ''.join(
            char if char not in self.forbidden_chars else '_' for char in strain_name)
        if strain_name != clean_strain:
            logger.warning(
                f'Forbidden char found in {strain_name}, will be cleaned to {clean_strain}')
        return clean_strain, strain_name

    def stat_from_pav(self, file):
        self.para_dict = {'Strict_core': [[], []], 'Core': [[], []],
                          'Soft_core': [[], []], 'Shell': [[], []], 'Cloud': [[], []]}
        with open(file, 'r') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                line = line.strip().split('\t')
                group = str(line[3])

    def stat_from_detail_pav(self, file):
        self.min_dict = {'Strict_core': [], 'Core': [],
                         'Soft_core': [], 'Shell': [], 'Cloud': []}
        self.uni_dict = {'Strict_core': [], 'Core': [],
                         'Soft_core': [], 'Shell': [], 'Cloud': []}
        self.mean_dict = {'Strict_core': [], 'Core': [],
                          'Soft_core': [], 'Shell': [], 'Cloud': []}
        self.var_dict = {'Strict_core': [], 'Core': [],
                         'Soft_core': [], 'Shell': [], 'Cloud': []}
        with open(file, 'r') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                line = line.strip().split('\t')
                group = str(line[3])

                self.min_dict[group].append(float(line[5]))
                self.mean_dict[group].append(float(line[6]))
                self.var_dict[group].append(float(line[7]))
                self.uni_dict[group].append(float(line[8]))
