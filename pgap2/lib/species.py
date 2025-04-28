import os
import pyfastani
from Bio import SeqIO
from tqdm import tqdm
from loguru import logger
from collections import defaultdict

from pgap2.utils.supply import tqdm_


class Species():
    def __init__(self, marker_file, strain_dict, ani, outdir) -> None:
        self.species_tag = {'darb_strain': '', 'outgroup_strain_list': []}
        self.gene_code = {}
        self.gene_len = defaultdict(list)
        self.strain_dict = strain_dict
        self.expect_ani = ani
        self.outdir = outdir
        self.outlier_dict = {}
        if marker_file and os.path.exists(marker_file):
            darb_strain, outgroup_strain_list = self._get_marker(marker_file)
            if darb_strain:
                logger.info(
                    f"Darb strain: {darb_strain} assigned from {marker_file}")
                found_darb = False
                for strain in self.strain_dict.values():
                    if strain.strain_name == darb_strain:
                        self.load_darb(darb_strain)
                        found_darb = True
                        break
                assert found_darb, f"Darb strain {darb_strain} not found in the strain dict"
            if outgroup_strain_list:
                self.species_tag['outgroup_strain_list'] = outgroup_strain_list
                logger.info(
                    f"outgroup strain list: {outgroup_strain_list} assigned from {marker_file}")
        self._ani_index_strain()

    def find_outlier(self, threads):
        outlier_list = []
        self.ani_dict = {}
        darb_record = self._get_darb_record()
        query = (bytes(record.seq) for record in darb_record.values())
        hits = self.mapper.query_draft(query, threads=threads)
        logger.info(f"Finding outlier with ANI < {self.expect_ani}")
        for hit in hits:
            logger.debug(
                f"Strain {self.get_darb()} vs {hit.name} with identity={hit.identity} matches={hit.matches} fragments={hit.fragments}")
            self.ani_dict[hit.name] = hit.identity
            if hit.identity < self.expect_ani:
                outlier_list.append(hit.name)
                logger.warning(
                    f"Outlier found: {hit.name} with ANI {hit.identity}")
        self.load_gene_outlier(outlier_list, 'ani')

    def load_genmoe_attr(self, genome_attr):
        self.genome_attr = genome_attr

    def load_gene_outlier(self, gene_outlier, method):
        self.outlier_dict.update({method: gene_outlier})

    def get_outlier(self):
        for method, outliers in self.outlier_dict.items():
            for strain in outliers:
                strain_name = self.strain_dict[strain].strain_name
                logger.warning(f"Outlier found by {method}: {strain_name}")
        return self.outlier_dict

    def get_ani(self, strain):
        return self.ani_dict[strain]

    def _ani_index_strain(self):
        sketch = pyfastani.Sketch()
        for strain in self.strain_dict:
            genome_file = f'{self.outdir}/genome_index/{strain}/ref.fa'
            if not os.path.exists(f"{self.outdir}/genome_index"):
                logger.error(
                    f"Genome index directory {self.outdir}/genome_index not exists")
                logger.error(
                    f"Please try the pipeline again in the new output directory")
                exit()

            records = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
            self.strain_dict[strain].genome = records
            sketch.add_draft(strain, (bytes(record.seq)
                             for record in records.values()))
        mapper = sketch.index()
        self.mapper = mapper

    def load_darb(self, darb_strain):
        assert darb_strain in self.strain_dict, f"Darb strain {darb_strain} not in the strain dict"
        self.species_tag['darb_strain'] = darb_strain
        self._index_darb()

    def _index_darb(self):
        darb_record = self._get_darb_record()
        sketch = pyfastani.Sketch()
        sketch.add_draft("darb", (bytes(record.seq)
                         for record in darb_record.values()))
        mapper = sketch.index()
        self.darb_mapper = mapper

    def get_total_query(self):
        return self.strain_dict.keys()

    def has_darb(self):
        return bool(self.species_tag['darb_strain'])

    def has_outgroup(self):
        return bool(self.species_tag['outgroup_strain_list'])

    def get_darb(self):
        return self.species_tag['darb_strain']

    def _get_darb_record(self):
        return self.strain_dict[self.species_tag['darb_strain']].genome

    def get_outgroup(self):
        return self.species_tag['outgroup_strain_list']

    def stat_gene_code(self, pg):
        with open(pg.annot_file, 'r') as fh:
            for line in tqdm(fh, desc=tqdm_.step(5), unit=' line', disable=pg.disable_tqdm, total=pg.total_gene_num+1):
                line = line.strip()
                if line.startswith("#"):
                    continue
                lines = line.split("\t")
                strain_index = int(lines[0].split(":")[0])
                nucl_fa = lines[-2]
                start_end = f'{nucl_fa[:3]}|{nucl_fa[-3:]}'
                if start_end not in self.gene_code:
                    self.gene_code[start_end] = defaultdict(int)
                self.gene_code[start_end][strain_index] += 1
                gene_len = len(lines[-1])
                self.gene_len[strain_index].append(gene_len)

    @staticmethod
    def _get_marker(marker_file):
        darb_strain = None
        outgroup_strain_list = []
        with open(marker_file, 'r') as fh:
            # 标志找到darb strain和outgroup strain list
            find_outgroup = False

            for line in fh:
                line = line.strip()

                if line.startswith("#"):
                    continue

                if line.lower().startswith("darb strain:"):
                    find_outgroup = False
                    darb_strain = line.split(":")[1].strip()
                    continue

                if line.lower().startswith("outgroup strain list:"):
                    find_outgroup = True
                    continue

                if find_outgroup and line:
                    outgroup_file = line.strip()
                    if not outgroup_file:
                        continue
                    assert os.path.exists(
                        outgroup_file), f"outgroup file {outgroup_file} not exists"
                    outgroup_strain_list.append(outgroup_file)

        return darb_strain, outgroup_strain_list
