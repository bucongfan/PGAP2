import bisect
import networkx as nx

from loguru import logger
from tqdm import tqdm
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

from pgap2.lib.pangenome import Pangenome
from pgap2.lib.tree import Tree
from pgap2.utils.tools import shortest_path_length_with_max_length, insert_node, merge_node
from pgap2.utils.supply import run_command
from pgap2.utils.supply import tqdm_, sfw

from collections import defaultdict
import tempfile
from Bio import BiopythonWarning
import warnings


def find_nearest_numbers(array, target):

    index = bisect.bisect_left(array, target)
    if index == 0:
        return array[0], array[1]
    if index == len(array):
        return array[-2], array[-1]
    before = array[index - 1]
    after = array[index]

    return (before, after) if target - before < after - target else (after, before)


def retrieve_from_paf(fpaf):
    need_retrieve_node = {}
    with open(fpaf) as fh:
        for line in fh:
            split_list = line.strip().split('\t')
            identity = round(int(
                split_list[9])/min(int(split_list[1])*3, int(split_list[10])), 5)
            cigar = split_list[19]
            contig = split_list[5]
            if contig not in need_retrieve_node:
                need_retrieve_node[contig] = {}
            need_retrieve_node[contig].update({split_list[0]: {'cigar': cigar,
                                                               'identity': identity,
                                                               'loc': {'start': int(split_list[7]), 'end': int(split_list[8]), 'strand': 1 if split_list[4] == '+' else -1}
                                                               }})
    return need_retrieve_node


def retrieve_gene(G: nx.Graph, pg: Pangenome, outdir, radius, threads):
    fh_annot = open(pg.annot_file, 'a')
    member2leaf = {}
    for node in G.nodes():
        for member in G.nodes[node]['members']:
            member2leaf[member] = node
    need_retrieve = defaultdict(list)
    need_retrieve_count = 0
    for node in G.nodes():
        for neigh in G.neighbors(node):
            # node need to retrieve
            for need_retrieve_strain in (G.nodes[neigh]['strains']-G.nodes[node]['strains']):
                need_retrieve_count += 1
                need_retrieve[need_retrieve_strain].append(node)

    logger.info(
        f'There are a total of {need_retrieve_count} nodes with possibilities to retrieve.')
    i = 0
    bar = tqdm(total=len(need_retrieve), desc=tqdm_.step(
        7), unit='strain ', disable=pg.disable_tqdm)
    maxlen_map = {}
    for strain in need_retrieve:
        bar.update()
        strain_seq_dict = {}
        # the retrieved node want to merge to this mode
        for node in need_retrieve[strain]:
            if not G.has_node(node):
                continue
            max_len = 0
            repre_member = None
            if node in maxlen_map:
                repre_member = maxlen_map[node]
            else:
                for member in G.nodes[node]['members']:
                    if pg.annot[member]['len'] > max_len:
                        max_len = pg.annot[member]['len']
                        repre_member = member
                maxlen_map[node] = repre_member
        repre_member2node = {v: k for k, v in maxlen_map.items()}

        with tempfile.NamedTemporaryFile(mode='w') as temp_file:
            for repre_member in maxlen_map.values():
                temp_file.write(f'{repre_member}\n')
            temp_file.flush()
            run_command(
                f'{sfw.seqtk} subseq {pg.prot_file} {temp_file.name} >{outdir}/genome_index/{strain}/query.fa')

        run_command(
            f'{sfw.miniprot} -t {threads} -N 1 -S {outdir}/genome_index/{strain}/ref.mpi {outdir}/genome_index/{strain}/query.fa >{outdir}/genome_index/{strain}/out.paf')
        need_retrieve_node = retrieve_from_paf(
            fpaf=f'{outdir}/genome_index/{strain}/out.paf')

        for contig in need_retrieve_node:
            retrieved_node = set()
            gene_rank = pg.gene_rank[f'{strain}:{contig}']
            gene_num = len(gene_rank)
            for j, repre_member in enumerate(need_retrieve_node[contig]):
                expect_node = repre_member2node[repre_member]
                # may happen when the node is already merged
                if not G.has_node(expect_node):
                    continue
                if expect_node in retrieved_node:
                    continue
                identity = need_retrieve_node[contig][repre_member]['identity']
                cigar = need_retrieve_node[contig][repre_member]['cigar']
                loc = need_retrieve_node[contig][repre_member]['loc']
                start = loc['start']
                end = loc['end']

                if identity < pg.para_id:
                    continue

                if len(gene_rank) <= 1:
                    continue

                before_s, after_s = find_nearest_numbers(gene_rank, start)
                before_s2, after_s2 = find_nearest_numbers(gene_rank, end)
                if (before_s != before_s2) or (after_s != after_s2):
                    logger.trace(
                        f'Retrived gene ({loc})) has overlap with exists coding gene, skip it.')
                    continue
                else:
                    before_e = before_s + \
                        pg.annot[f'{strain}:{contig}:{gene_rank.index(before_s)}']['len']*3
                    if before_e >= start:
                        logger.trace(
                            f'Retrived gene ({loc})) has overlap with exists coding gene, skip it.')
                        continue

                assumed_gene = f'{strain}:{contig}:{gene_num+j}'
                before_node = member2leaf[f'{strain}:{contig}:{gene_rank.index(before_s)}']
                after_node = member2leaf[f'{strain}:{contig}:{gene_rank.index(after_s)}']

                if not G.has_node(before_node) or not G.has_node(after_node):
                    # this gene has been merged baited by other strain's gene
                    continue

                G = insert_node(G, before_node,
                                after_node,
                                (assumed_gene,
                                 {
                                     'mci': identity,
                                     'uni': identity,
                                     'length': 1,
                                     'members': [assumed_gene],
                                     'strains': set([strain]), 'has_para': False}
                                 ))

                _, path = shortest_path_length_with_max_length(
                    G, expect_node, assumed_gene, {}, radius)
                if not path:
                    need_merge = False
                else:
                    need_merge = True

                if need_merge:
                    # No more need other method to judge such as synteny, tolerate diff.
                    # Because retrieve method is only used in retrieve gene to given expected node but not found new gene
                    i += 1
                    retrieve_gene_name = f'retrieve_{i}'
                    retrieved_node.add(expect_node)
                    if '_' in expect_node:
                        expect_node_name = pg.annot[expect_node.split('_')[
                            1]]['id']
                    else:
                        expect_node_name = pg.annot[expect_node]['id']

                    G = merge_node(
                        G, pg, None, [assumed_gene, expect_node], target=expect_node)
                    # get seqeuence of this retrieved gene using biopython
                    if not strain_seq_dict:
                        strain_seq_dict = SeqIO.to_dict(
                            SeqIO.parse(f'{outdir}/genome_index/{strain}/ref.fa', 'fasta'))
                    feature = SeqFeature(FeatureLocation(
                        start=start-1, end=end, strand=loc['strand']), type="CDS", id=retrieve_gene_name)
                    seq = feature.extract(strain_seq_dict[contig].seq)
                    with warnings.catch_warnings(record=True) as w:
                        # 只捕获 BiopythonWarning
                        warnings.simplefilter("always", BiopythonWarning)

                        prot = seq.translate()

                        # 检查是否捕获到任何警告
                        if len(w) > 0 and issubclass(w[0].category, BiopythonWarning):
                            pseudo_gene_product = f'pseudo|{identity}|{cigar}'
                        else:
                            pseudo_gene_product = f'coding|{identity}|{cigar}'

                    contig_name = pg.annot_contig_map[f'{strain}:{contig}']
                    strain_name = pg.strain_dict[strain].strain_name
                    fh_annot.write(
                        f'#{assumed_gene}\t{strain_name}\t{contig_name}\t{feature.location}\t{len(prot)}\t{retrieve_gene_name}\t{expect_node_name}\t{pseudo_gene_product}\t{seq}\t{prot}\n')
                    pg.annot[assumed_gene] = {'id': retrieve_gene_name, 'name': expect_node_name, 'product': pseudo_gene_product,
                                              'len': len(prot)}
                    logger.debug(
                        f'Retrieved {retrieve_gene_name} in {loc} with identity {identity}, and merged to {expect_node_name}.')
                else:
                    G.remove_node(assumed_gene)
    fh_annot.close()
    logger.info(f'A total of {i} genes were found.')
    return G, pg
