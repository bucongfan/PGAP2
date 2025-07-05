import os
import re
import subprocess

from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from tajimas_d import tajimas_d, watterson_estimator, pi_estimator

from tqdm import tqdm
from loguru import logger
from datetime import datetime
from multiprocessing import get_context

# from Bio.Align.Applications import MafftCommandline
# from Bio.Align.Applications import MuscleCommandline
# from Bio.Align.Applications import TCoffeeCommandline

from pgap2.lib.basic import Basic
from pgap2.utils.supply import sfw, tqdm_
from pgap2.utils.supply import set_golbal


def set_logger(logger_):
    global logger
    logger = logger_


class Phylogeny():
    def __init__(self, basic, outdir, threads, disable, msa_method, tree_method, fastbaps_levels=None, fastbaps_prior=None, add_paras=[]) -> None:
        self.basic: Basic = basic
        self.outdir = outdir
        self.threads = threads
        self.disable = disable
        self.msa_method = msa_method
        self.tree_method = tree_method
        self.fastbaps_levels = fastbaps_levels
        self.fastbaps_prior = fastbaps_prior
        self.add_paras_dict = {2: '', 4: '',
                               6: '', 7: '', 8: '', 9: '', 10: ''}
        self.load_additional_paras(add_paras)

    def load_additional_paras(self, add_paras):
        if not add_paras:
            return
        for this_str in add_paras:
            step = this_str[0]
            seprator = this_str[1]
            paras = this_str[2:]
            if seprator != ':':
                logger.error(
                    f"Seprator should be ':', but got {seprator} in {this_str} when loading additional parameters")
                raise ValueError(
                    logger.error(f"Seprator should be ':', but got {seprator} in {this_str} when loading additional parameters"))
            try:
                # step should be int
                step = int(step)
            except ValueError:
                logger.error(
                    f"Step should be int, but got {step} in {this_str} when loading additional parameters")
                raise ValueError(
                    logger.error(f"Step should be int, but got {step} in {this_str} when loading additional parameters"))

            if step not in self.add_paras_dict:
                logger.error(
                    f"Step should be in [2,4,6,7,8,9,10] because only extra software used in these steps, but got {step} in {this_str} when loading additional parameters")
                raise ValueError(
                    logger.error(f"Step should be in [2,4,6,7,8,9,10] because only extra software used in these steps, but got {step} in {this_str} when loading additional parameters"))
            logger.info(
                f"Additional parameters loaded for step [{step}]: {paras}")
            self.add_paras_dict[step] = paras

    def check_before_run(self, cmd, outdir):
        worksh = os.path.join(outdir, 'work.sh')
        if os.path.exists(worksh):
            with open(worksh, 'r') as f:
                if cmd in f.read():
                    for result in self.results_file:
                        if not os.path.exists(result) or os.path.getsize(result) == 0:
                            return False
                    return True
        return False

    def dump_sequences(self, wd):
        self.results_file = []
        outdir_cds = os.path.join(wd, '01.gene_cds')
        os.makedirs(outdir_cds, exist_ok=True)
        outdir_prot = os.path.join(wd, '01.gene_prot')
        os.makedirs(outdir_prot, exist_ok=True)
        for cluster, seqs in tqdm(self.basic.used_cluster.items(), desc=tqdm_.step(1)):
            path_cds = os.path.join(outdir_cds, f'{cluster}.fa')
            path_prot = os.path.join(outdir_prot, f'{cluster}.fa')
            if os.path.exists(path_cds) and os.path.exists(path_prot) and os.path.getsize(path_cds) > 0 and os.path.getsize(path_prot) > 0:
                logger.debug(
                    f"File already exists {path_cds} and {path_prot}. Skip.")
            else:
                with open(path_cds, 'w') as fh_cds, open(path_prot, 'w') as fh_prot:
                    for seq in seqs:
                        SeqIO.write(seq, fh_cds, 'fasta')
                        protein_seq = seq.translate(to_stop=True)
                        protein_record = SeqRecord(
                            protein_seq.seq, id=seq.id, description=seq.description)
                        SeqIO.write(protein_record, fh_prot, 'fasta')
        self.sequence_cds_path = outdir_cds
        self.sequence_prot_path = outdir_prot
        self.results_file.append(outdir_cds)
        self.results_file.append(outdir_prot)
        logger.info(f'cds: {outdir_cds}')
        logger.info(f'prot: {outdir_prot}')

    @staticmethod
    def _msa(params):
        msa_method, cluster, output, cmd = params
        if os.path.exists(output):
            logger.debug(f"File already exists {output}. Skip.")
        else:
            process = subprocess.run(
                cmd, shell=True, capture_output=True, text=True)
            if process.returncode != 0:
                logger.error(
                    f"Error occurred while running MSA: {process.stderr}")
                raise Exception(
                    f"Error occurred while running MSA: {process.stderr}")
            if msa_method == 'mafft' and process.returncode == 0:
                with open(output, "w") as handle:
                    handle.write(process.stdout)
        return ((cluster, output))

    def _generate_msa_commands(self, prot_align_outdir):
        add_paras = self.add_paras_dict[2]
        commands = []
        for cluster in self.basic.used_cluster:
            expected_output = os.path.join(
                prot_align_outdir, cluster + '.aln.fa')
            input_prot = os.path.join(self.sequence_prot_path, cluster + '.fa')
            if self.msa_method == 'mafft':
                # cline = MafftCommandline(
                #     input=input_prot, auto=True, thread=1, quiet=True)
                cline = f'{sfw.mafft} --auto --thread 1 {input_prot} {add_paras}'
            elif self.msa_method == 'muscle':
                # cline = MuscleCommandline(
                #     input=input_prot, out=expected_output, quiet=True)
                cline = f'{sfw.muscle} -in {input_prot} -out {expected_output} -quiet {add_paras}'
            elif self.msa_method == 'tcoffee':
                # cline = TCoffeeCommandline(
                #     infile=input_prot, output='fasta', outfile=expected_output, quiet=True)
                # print(str(cline))
                cline = f'{sfw.tcoffee} {input_prot} -output fasta -outfile {expected_output} -quiet {add_paras}'
            commands.append((self.msa_method, cluster,
                            expected_output, str(cline)))
        return commands

    def msa(self, wd):
        prot_align_outdir = os.path.join(wd, '02.msa_prot')
        os.makedirs(prot_align_outdir, exist_ok=True)
        commands = self._generate_msa_commands(prot_align_outdir)
        results = {}
        self.results_file = []
        with get_context('fork').Pool(self.threads, initializer=set_golbal, initargs=(logger,)) as p:
            for cluster, output in p.imap_unordered(self._msa, tqdm(commands, desc=tqdm_.step(2))):
                results[cluster] = output
        self.results_file.append(prot_align_outdir)
        return results

    @staticmethod
    def _codon_alignment(params):
        cluster, protein_alignment_file, original_dna_file, expected_output = params
        if os.path.exists(expected_output) and os.path.getsize(expected_output) > 0:
            logger.debug(f"File already exists {expected_output}. Skip.")
        else:
            original_dna = list(SeqIO.parse(original_dna_file, "fasta"))
            protein_alignment = list(SeqIO.parse(
                protein_alignment_file, "fasta"))
            dna_dict = {record.id: record.seq for record in original_dna}
            aligned_dna = []
            for prot_record in protein_alignment:
                dna_seq = dna_dict[prot_record.id]
                codon_aligned_seq = []
                codon_index = 0
                for aa in prot_record.seq:
                    if aa != "-":
                        codon = dna_seq[codon_index:codon_index+3]
                        codon_index += 3
                    else:
                        codon = "---"
                    codon_aligned_seq.append(str(codon))
                aligned_dna_seq = Seq("".join(codon_aligned_seq))
                aligned_dna.append(SeqRecord(
                    aligned_dna_seq, id=prot_record.id, description=prot_record.description))
            with open(expected_output, "w") as handle:
                for record in aligned_dna:
                    SeqIO.write(record, handle, "fasta")
        return ((cluster, expected_output))

    def _generate_codon_alignment_commands(self, codon_align_outdir):
        commands = []
        for cluster in self.msa_results:
            expected_output = os.path.join(
                codon_align_outdir, cluster + '.aln.codon.fa')
            original_dna_file = os.path.join(
                self.sequence_cds_path, cluster + '.fa')
            protein_alignment_file = self.msa_results[cluster]
            commands.append((cluster, protein_alignment_file,
                            original_dna_file, expected_output))
        return commands

    def codon_alignment(self, wd):
        self.results_file = []
        results = {}
        codon_align_outdir = os.path.join(wd, '03.codon_alignment')
        os.makedirs(codon_align_outdir, exist_ok=True)
        cmds = list(self._generate_codon_alignment_commands(
            codon_align_outdir))
        with get_context('fork').Pool(self.threads, initializer=set_golbal, initargs=(logger,)) as p:
            for cluster, exp_output in p.imap_unordered(self._codon_alignment, tqdm(cmds, desc=tqdm_.step(3))):
                results[cluster] = exp_output

        self.results_file.append(codon_align_outdir)
        return results

    @staticmethod
    def _trim_alignment(params):
        cluster, codon_aln, trim_align_outdir, trimmer, add_paras = params
        output = os.path.join(
            trim_align_outdir, cluster + '.aln.codon.trimmed.fa')
        if os.path.exists(output) and os.path.getsize(output) > 0:
            logger.debug(f"File already exists {output}. Skip.")
        else:
            cline = '{} {} --output {} --codon --quiet {}'.format(
                    trimmer, codon_aln, output, add_paras)
            process = subprocess.run(
                cline, shell=True, capture_output=True, text=True)
            if process.returncode != 0:
                logger.error(
                    f"Error occurred while trimming alignment: {process.stderr}")
                raise Exception(
                    f"Error occurred while trimming alignment: {process.stderr}")
        return ((cluster, output))

    def _generate_trim_alignment_commands(self, trim_align_outdir):
        commands = []
        add_paras = self.add_paras_dict[4]
        for cluster, codon_aln in self.codon_results.items():
            commands.append(
                (cluster, codon_aln, trim_align_outdir, sfw.clipkit, add_paras))
        return commands

    def trim_alignment(self, wd):
        self.results_file = []
        results = {}
        trim_align_outdir = os.path.join(wd, '04.trim_alignment')
        os.makedirs(trim_align_outdir, exist_ok=True)
        # cmds=list(self._generate_codon_alignment_commands(codon_align_outdir))
        # with get_context('fork').Pool(self.threads,initializer=set_golbal,initargs=(logger,)) as p:
        #     for cluster, exp_output in p.imap_unordered(self._codon_alignment, tqdm(cmds,desc=tqdm_.step(3))):
        #         results[cluster] = exp_output

        # self.results_file.append(codon_align_outdir)
        cmds = list(self._generate_trim_alignment_commands(trim_align_outdir))
        with get_context('fork').Pool(self.threads, initializer=set_golbal, initargs=(logger,)) as p:
            for cluster, output in p.imap_unordered(self._trim_alignment, tqdm(cmds, desc=tqdm_.step(4))):
                results[cluster] = output
        self.results_file.append(trim_align_outdir)
        return results

    def _tajimas_d_test(self, aln_file):
        # load alignment file with fasta format using biopython
        alignment = AlignIO.read(aln_file, "fasta")
        sequences = [str(rec.seq) for rec in alignment]
        cluster_name = os.path.basename(aln_file).split('.')[0]
        err = False
        try:
            # calculate tajima's D
            tajimas_d_value = tajimas_d(sequences)
            # calculate watterson estimator
            watterson_estimator_value = watterson_estimator(sequences)
            # calculate pi estimator
            pi_estimator_value = pi_estimator(sequences)
            err = False
        except:
            tajimas_d_value = watterson_estimator_value = pi_estimator_value = 0
            err = True

        return (cluster_name, (tajimas_d_value, watterson_estimator_value, pi_estimator_value), err)

    def tajimas_d_test(self, wd):
        # This is the last step of the pipeline, so we don't need to return the results and check before run
        self.results_file = []
        results = {}
        tajimas_d_outdir = os.path.join(wd, '05.tajimas_d')
        os.makedirs(tajimas_d_outdir, exist_ok=True)
        all_files = []
        logger.info('Reading all alignment files...')
        for each_file in os.listdir(os.path.join(wd, '04.trim_alignment')):
            each_file = os.path.join(wd, '04.trim_alignment', each_file)
            if os.path.isfile(each_file) and each_file.endswith('.fa'):
                all_files.append(each_file)
        logger.info('Calculating Tajima\'s D...')
        with get_context('fork').Pool(self.threads, initializer=set_golbal, initargs=(logger,)) as p:
            for cluster, exp_output, err in p.imap_unordered(self._tajimas_d_test, tqdm(all_files, desc=tqdm_.step(5), disable=self.disable)):
                if err:
                    logger.error(
                        f"Error occurred while calculating Tajima's D for {cluster}. Use 0 instead.")
                results[cluster] = exp_output
        logger.info('Writing results...')
        with open(os.path.join(tajimas_d_outdir, 'tajimas_d.txt'), 'w') as fh:
            fh.write('Cluster\tTajimas_D\tWatterson_estimator\tPi_estimator\n')
            for cluster, values in results.items():
                fh.write(f'{cluster}\t{values[0]}\t{values[1]}\t{values[2]}\n')

        self.results_file.append(os.path.join(
            tajimas_d_outdir, 'tajimas_d.txt'))

    def concatenate_alignment(self, wd):
        self.results_file = []
        concatenated_sequences = {}
        cat_align_outdir = os.path.join(wd, '05.concatenate_alignment')
        os.makedirs(cat_align_outdir, exist_ok=True)
        concatenated_gb = SeqRecord(Seq(''),
                                    id='concatenated', description='postprocess phylogeny concatenated alignment')
        concatenated_gb.annotations['division'] = 'BAC'
        concatenated_gb.annotations['source'] = 'pgap2'
        concatenated_gb.annotations['organism'] = 'Bacteria'
        concatenated_gb.annotations['taxonomy'] = ['Bacteria']
        concatenated_gb.annotations['date'] = datetime.now().strftime(
            "%d-%b-%Y").upper()
        concatenated_gb.annotations['comment'] = 'Concatenated alignment of core genes'
        concatenated_gb.annotations['molecule_type'] = 'DNA'

        position_offset = 0  # Track the concatenated sequence length

        real_involved_clusters = {
            strain: 0 for strain in self.basic.strain_dict}
        for core_cluster in tqdm(self.basic.phylogeny_dict, desc=tqdm_.step(5)):
            if self.basic.phylogeny_dict[core_cluster]['Type'] != 'Core':
                continue
            core_cluster_path = self.trim_results[core_cluster]
            core_records = list(AlignIO.read(core_cluster_path, "fasta"))
            example_length = len(core_records[0].seq)

            tmp_dict = {}  # Temporary dictionary to hold sequences for this cluster
            for record in core_records:
                sequence = record.seq
                strain = int(record.id.split(':')[0])

                if strain in tmp_dict:
                    # That would happen if the para_strategy is best and the cluster has paralogs
                    if tmp_dict[strain].count('-') > sequence.count('-'):
                        tmp_dict[strain] = sequence
                else:
                    tmp_dict[strain] = sequence

            # Add sequences or gaps to concatenated_sequences
            for strain in self.basic.strain_dict:
                if strain not in tmp_dict:
                    seq = '-' * example_length
                else:
                    seq = str(tmp_dict[strain])
                    real_involved_clusters[strain] += 1

                if strain not in concatenated_sequences:
                    concatenated_sequences[strain] = seq
                else:
                    concatenated_sequences[strain] += seq

            # Add feature for this cluster
            feature = SeqFeature(FeatureLocation(start=position_offset,
                                                 end=position_offset + example_length),
                                 type="region",
                                 qualifiers={"note": core_cluster_path, "locus_tag": core_cluster})
            concatenated_gb.features.append(feature)
            position_offset += example_length  # Update position for the next cluster

        # Save the concatenated record to a GenBank file
        SeqIO.write(concatenated_gb, os.path.join(
            cat_align_outdir, 'core_gene_alignment.gb'), 'genbank')

        with open(os.path.join(cat_align_outdir, 'core_gene_alignment.aln'), 'w') as fh, open(os.path.join(cat_align_outdir, 'core_gene_alignment.txt'), 'w') as fh2:
            fh2.write('#Clean_strain\tStrain\tSymbol\tInvolved_genes\n')
            for strain, seq in concatenated_sequences.items():
                real_strain, original_strain_name = self.basic.get_real_strain_name(
                    strain)
                fh2.write(
                    f'{real_strain}\t{original_strain_name}\t{strain}\t{real_involved_clusters[strain]}\n')
                involved = real_involved_clusters[strain]
                if involved == 0:
                    raise ValueError(
                        logger.error('BUG: No sequences for strain {real_strain} in core genes. Tell me in github'))
                seq = re.sub(r'[^ATCG-]', 'N', seq.upper())
                record = SeqRecord(Seq(
                    seq), id=real_strain, description='')
                SeqIO.write(record, fh, 'fasta')
        self.results_file.append(f'{cat_align_outdir}/core_gene_alignment.aln')
        self.results_file.append(f'{cat_align_outdir}/core_gene_alignment.gb')
        self.results_file.append(f'{cat_align_outdir}/core_gene_alignment.txt')
        return f'{cat_align_outdir}/core_gene_alignment.aln'

    def construct_tree(self, wd):
        _tqdm = tqdm(total=1, desc=tqdm_.step(6))
        self.results_file = []
        tree_outdir = os.path.join(wd, '06.core_alignment_tree')
        os.makedirs(tree_outdir, exist_ok=True)
        tree_path = os.path.join(tree_outdir, 'core_gene_alignment')
        add_paras = self.add_paras_dict[6]
        if self.tree_method == 'raxml':
            cmd = f'{sfw.raxml} --msa {self.first_aln} --model GTR+G --prefix {tree_path} --all --bs-trees 100 --redo {add_paras}'
            best_tree = os.path.join(
                tree_outdir, 'core_gene_alignment'+'.raxml.bestTree')
        elif self.tree_method == 'fasttree':
            cmd = f'{sfw.fasttree} -gtr -gamma -nt -out {tree_path}.treefile {self.first_aln} {add_paras} 2> {tree_path}.log'
            best_tree = os.path.join(
                tree_outdir, 'core_gene_alignment'+'.treefile')
        elif self.tree_method == 'iqtree':
            cmd = f'{sfw.iqtree} -s {self.first_aln} -m GTR+G -pre {tree_path} -bb 1000 -redo -nt AUTO {add_paras}'
            best_tree = os.path.join(
                tree_outdir, 'core_gene_alignment'+'.treefile')

        self.results_file.append(best_tree)

        finished_status = self.check_before_run(cmd, tree_outdir)
        if finished_status:
            logger.warning(f"Step 6 already finished. Skip.")
        else:
            logger.info(f'Run [{cmd}]')
            with open(os.path.join(tree_outdir, 'work.sh'), 'w') as f:
                f.write(cmd)
            subprocess.run(cmd, shell=True, check=True)

        _tqdm.update(1)
        _tqdm.close()
        return best_tree

    def recombination_inference(self, wd):
        _tqdm = tqdm(total=1, desc=tqdm_.step(7))
        self.results_file = []
        recombination_outdir = os.path.join(wd, '07.recombination')
        tree_path = self.first_tree
        aln_path = self.first_aln
        os.makedirs(recombination_outdir, exist_ok=True)
        add_paras = self.add_paras_dict[7]
        cmd = f'{sfw.cfml} {add_paras} {tree_path} {aln_path} {recombination_outdir}/recombination'

        self.results_file.append(
            f'{recombination_outdir}/recombination.ML_sequence.fasta')
        self.results_file.append(
            f'{recombination_outdir}/recombination.position_cross_reference.txt')
        self.results_file.append(
            f'{recombination_outdir}/recombination.em.txt')
        self.results_file.append(
            f'{recombination_outdir}/recombination.importation_status.txt')
        self.results_file.append(
            f'{recombination_outdir}/recombination.labelled_tree.newick')

        finished_status = self.check_before_run(cmd, recombination_outdir)
        if finished_status:
            logger.warning(f"Step 7 already finished. Skip.")
        else:
            logger.info(f'Run [{cmd}]')
            with open(os.path.join(recombination_outdir, 'work.sh'), 'w') as f:
                f.write(cmd)
            subprocess.run(cmd, shell=True, check=True)

        _tqdm.update(1)
        _tqdm.close()
        return f'{recombination_outdir}/recombination'

    def mask_recombination(self, wd):
        _tqdm = tqdm(total=1, desc=tqdm_.step(8))
        self.results_file = []
        recombination_outdir = os.path.join(wd, '08.mask_recombination')
        os.makedirs(recombination_outdir, exist_ok=True)
        add_paras = self.add_paras_dict[8]
        cmd = f'{sfw.maskrc} {add_paras} --aln {self.first_aln} --out {recombination_outdir}/maskrc.aln --regions {recombination_outdir}/recombinant_regions.txt --svg {recombination_outdir}/recombinant_regions.svg {self.recombination_prefix}'

        self.results_file.append(f'{recombination_outdir}/maskrc.aln')
        self.results_file.append(
            f'{recombination_outdir}/recombinant_regions.txt')
        self.results_file.append(
            f'{recombination_outdir}/recombinant_regions.svg')

        finished_status = self.check_before_run(cmd, recombination_outdir)
        if finished_status:
            logger.warning(f"Step 8 already finished. Skip.")
        else:
            logger.info(f'Run [{cmd}]')
            with open(os.path.join(recombination_outdir, 'work.sh'), 'w') as f:
                f.write(cmd)
            subprocess.run(cmd, shell=True, check=True)

        _tqdm.update(1)
        _tqdm.close()
        return f'{recombination_outdir}/maskrc.aln'

    def reconstruct_tree(self, wd):
        _tqdm = tqdm(total=1, desc=tqdm_.step(9))
        self.results_file = []
        tree_outdir = os.path.join(wd, '09.reconstructed_tree')
        os.makedirs(tree_outdir, exist_ok=True)
        tree_path = os.path.join(tree_outdir, 'reconstructed_tree')
        add_paras = self.add_paras_dict[9]

        if self.tree_method == 'raxml':
            cmd = f'{sfw.raxml} --msa {self.maskrc_aln} --model GTR+G --prefix {tree_path} --all --bs-trees 100 --redo {add_paras}'
            best_tree = os.path.join(
                tree_outdir, 'core_gene_alignment'+'.raxml.bestTree')
        elif self.tree_method == 'fasttree':
            cmd = f'{sfw.fasttree} -gtr -gamma -nt -out {tree_path}.treefile {self.maskrc_aln} {add_paras} 2> {tree_path}.log'
            best_tree = os.path.join(
                tree_outdir, 'core_gene_alignment'+'.treefile')
        elif self.tree_method == 'iqtree':
            cmd = f'{sfw.iqtree} -s {self.maskrc_aln} -m GTR+G -pre {tree_path} -bb 1000 -redo -nt AUTO {add_paras}'
            best_tree = os.path.join(
                tree_outdir, 'core_gene_alignment'+'.treefile')

        self.results_file.append(best_tree)
        finished_status = self.check_before_run(cmd, tree_outdir)
        if finished_status:
            logger.warning(f"File already exists {best_tree}. Skip.")
        else:
            logger.info(f'Run [{cmd}]')
            with open(os.path.join(tree_outdir, 'work.sh'), 'w') as f:
                f.write(cmd)
            subprocess.run(cmd, shell=True, check=True)

        _tqdm.update(1)
        _tqdm.close()
        return best_tree

    def genetic_clustering(self, wd):
        # This is the last step of baps, so we don't need to check before run
        self.results_file = []
        clustering_outdir = os.path.join(wd, '10.genetic_clustering')
        os.makedirs(clustering_outdir, exist_ok=True)
        add_paras = self.add_paras_dict[10]
        cmd = f'python {sfw.fastbaps} {add_paras} --input {self.maskrc_aln} --phylogeny {self.masked_tree} --out {clustering_outdir}/fastbaps_clusters.csv --threads {self.threads} --levels {self.fastbaps_levels} --prior {self.fastbaps_prior}'
        subprocess.run(cmd, shell=True, check=True)
        logger.info(f'Run [{cmd}]')

        self.results_file.append(f'{clustering_outdir}/fastbaps')
        return f'{clustering_outdir}/fastbaps'

    def start_at(self, step: int):
        wd = f'{self.outdir}/postprocess_phylogeny'
        if step == 1:
            logger.info(f'01.Extract core cds and proteins')
            self.dump_sequences(wd)
        elif step == 2:
            logger.info(
                f'02.Multiple sequence alignment using {self.msa_method}')
            self.msa_results = self.msa(wd)
        elif step == 3:
            logger.info(f'03.codon alignment')
            self.codon_results = self.codon_alignment(wd)
        elif step == 4:
            logger.info(f'04.Trim alignment')
            self.trim_results = self.trim_alignment(wd)
        elif step == 5:
            logger.info(f'05.Core alignment concatenation')
            self.first_aln = self.concatenate_alignment(wd)
        elif step == 6:
            logger.info(
                f'06.Phylogenetic tree construction using {self.tree_method}')
            self.first_tree = self.construct_tree(wd)
        elif step == 7:
            logger.info(f'07.Recombination inference')
            self.recombination_prefix = self.recombination_inference(wd)
        elif step == 8:
            logger.info(f'08.Musk the recombination region')
            self.maskrc_aln = self.mask_recombination(wd)
        elif step == 9:
            logger.info(f'09.Reconstruct the phylogeny tree')
            self.masked_tree = self.reconstruct_tree(wd)
        elif step == 10:
            logger.info(f'10.Genetic clustering')
            self.genetic_clustering(wd)
            pass
        elif step == 11:
            logger.info(f"5.Apply Tajima's D Test")
            self.tajimas_d_test(wd)
            pass

    def dump_results(self):
        for i, source_path in enumerate(self.results_file):
            source_file_name = os.path.basename(source_path)
            target_path = os.path.join(self.basic.outdir, source_file_name)
            if not os.path.exists(target_path):
                os.symlink(source_path, target_path)
                logger.success(f"result [{i+1}]: {target_path}")
            else:
                logger.warning(
                    f"File already exists {target_path}. Skip symlink.")
                logger.success(f"result [{i+1}]: {source_path}")
