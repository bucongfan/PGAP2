import logging
import itertools
import os
import bisect
import tempfile
import warnings
import networkx as nx
from collections import defaultdict

from loguru import logger
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import BiopythonWarning
from Bio.SeqFeature import SeqFeature, FeatureLocation

from pgap2.lib.pangenome import Pangenome
from pgap2.lib.tree import Tree
from pgap2.utils.generate_tree import run_alignment, process_lines
from pgap2.utils.tools import insert_node, merge_node
from pgap2.utils.supply import run_command
from pgap2.utils.supply import tqdm_, sfw


def set_logger(logger_):
    global logger
    logger = logger_


def find_nearest_numbers(array, target):
    if len(array) < 2:
        return (array[0], array[0])

    index = bisect.bisect_left(array, target)
    if index == 0:
        return array[0], array[1]
    if index == len(array):
        return array[-2], array[-1]
    before = array[index - 1]
    after = array[index]

    return (before, after)


def retrieve_from_paf(fpaf, min_overlap_bp=1):
    """
    Read miniprot PAF and return dict[contig][query] -> {...}

    filter out overlapping sections on the same contig and retain the one with the higher score

    Parameters
    ----------
    fpaf : str
        PAF file path

    min_overlap_bp : int, optional
        Minimum overlap in base pairs to consider for deduplication, by default 1

    Returns
    -------
    dict
        {contig: {query: {...}}}
    """

    keep = {}  # {contig: [{record dict}, ...]}
    with open(fpaf) as fh:
        for line in fh:
            if not line.strip():
                continue
            split_list = line.rstrip('\n').split('\t')

            query = split_list[0]
            qlen = int(split_list[1]) * 3
            strand = 1 if split_list[4] == '+' else -1
            contig = split_list[5]
            t_start = int(split_list[7])
            t_end = int(split_list[8])
            matches = int(split_list[9])
            block = int(split_list[10])

            # get score and identity
            score = None
            for f in split_list[12:]:
                if f.startswith('AS:i:'):
                    score = int(f.split(':')[-1])
                    break
            if score is None:
                raise ValueError(
                    f"PAF line {line.strip()} does not contain AS:i: score, using matches {matches} as fallback."
                )

            identity = round(
                qlen / block, 5) if block >= qlen else round(matches / block, 5)

            rec = dict(q=query, start=t_start, end=t_end,
                       score=score, identity=identity,
                       cigar=split_list[19] if len(split_list) > 19 else '',
                       strand=strand)

            # deduplication
            if contig not in keep:
                keep[contig] = [rec]
                continue

            overlap_found = False
            for idx, old in enumerate(keep[contig]):
                ov_len = min(t_end, old['end']) - max(t_start, old['start'])
                if ov_len >= min_overlap_bp:
                    overlap_found = True
                    logger.debug(
                        f"Found overlap on {contig}:{t_start}-{t_end} with existing {old['q']} ({old['start']}-{old['end']})"
                    )
                    better = False
                    if score > old['score']:
                        better = True
                    elif score == old['score']:
                        if identity > old['identity']:
                            better = True
                        elif identity == old['identity']:
                            if block > (old['end'] - old['start']):
                                better = True

                    if better:
                        logging.info(
                            f"Replace overlap on {contig}:{old['start']}-{old['end']} "
                            f"(query {old['q']}) with better hit {query}"
                        )
                        keep[contig][idx] = rec
                    break

            if not overlap_found:
                keep[contig].append(rec)

    result = {}
    for ctg, lst in keep.items():
        result[ctg] = {}
        for r in lst:
            result[ctg][r['q']] = {
                'cigar': r['cigar'],
                'identity': r['identity'],
                'score': r['score'],
                'loc': {'start': r['start'],
                        'end': r['end'],
                        'strand': r['strand']}
            }
    return result


def merge_sorted_lists(a, b):
    result = []
    i = j = 0
    while i < len(a) and j < len(b):
        if a[i] <= b[j]:
            result.append(a[i])
            i += 1
        else:
            result.append(b[j])
            j += 1
    result.extend(a[i:])
    result.extend(b[j:])
    return result


def realign_parser(outdir, evalue, aligner, mcxdeblast, mcl, ID, LD, AL, AS, threads=1, attribution: dict = {}):

    alignment_result = run_alignment(
        input_file=f'{outdir}/seq.fa', outdir=outdir, threads=threads, evalue=evalue, aligner=aligner)

    filtered_alignment_result = f'{outdir}/{aligner}.filtered.result'

    edges = []
    with open(alignment_result, 'r') as f, open(filtered_alignment_result, 'w') as out_f:
        lines = []
        for line in f:
            if not line.strip():
                continue
            query, subject = line.strip().split('\t')[:2]
            if attribution[query] != attribution[subject]:
                continue
            lines.append(line)

        filtered_lines = process_lines(ID=ID, LD=LD,
                                       AL=AL, AS=AS, lines=lines)
        out_f.writelines(filtered_lines[0])
        edges = filtered_lines[1]
    if edges:
        edges_dict = dict(((x, y), v) for x, y, v in edges)
        mcl_result = f'{outdir}/mcl.result'
        filtered_edges = set()
        run_command(
            f"{mcxdeblast} -m9 --score r --line-mode=abc {filtered_alignment_result} 2> /dev/null | {mcl} - --abc -I 1.5 -te 1 -o {mcl_result} &>/dev/null")
        with open(mcl_result, 'r') as mcl_f:
            for mcl_line in mcl_f:
                mcl_line = mcl_line.strip()
                if not mcl_line:
                    continue
                nodes = mcl_line.split()
                for a, b in itertools.combinations(nodes, 2):
                    if (a, b) in edges_dict:
                        filtered_edges.add((a, b, edges_dict[(a, b)]))
                        del edges_dict[(a, b)]
                    elif (b, a) in edges_dict:
                        filtered_edges.add((b, a, edges_dict[(b, a)]))
                        del edges_dict[(b, a)]

        return filtered_edges, edges_dict
    else:
        logger.debug(
            f'No edges found in {outdir}, skipping re-alignment for this cluster.')
        return None, None


def realign(pg: Pangenome, tree: Tree, outdir, need_realign_root: set):
    logger.info('Dumping query sequences...')
    exitst_file_flag = False
    with open(f'{outdir}/seq.fa', 'w') as fh:
        for root in tqdm(list(need_realign_root), desc=tqdm_.step(8), unit='group', disable=pg.disable_tqdm):
            nodes = tree.root_leaf[root]
            if len(nodes) <= 2:
                logger.debug(
                    f'Only one node ({nodes}) in the cluster {root}, skipping re-alignment.')
                need_realign_root.remove(root)
                continue
            for node in nodes:
                if node in pg.nucl_fa:
                    exitst_file_flag = True
                    nucl_seq = Seq(pg.nucl_fa[node])
                    fh.write(f'>{node}\n{nucl_seq.translate()}\n')
    if exitst_file_flag:
        logger.info('Re-aligning sequences to update the identity network...')
        edges, notused_edges = realign_parser(outdir=outdir, evalue=pg.evalue, aligner=pg.aligner, mcxdeblast=sfw.mcxdeblast,
                                              mcl=sfw.mcl, ID=pg.para_id, LD=pg.LD, AL=pg.AL, AS=pg.AS, threads=pg.threads, attribution=tree.leaf_root)
        if edges is not None:
            for node_a, node_b, weight in edges:
                assert tree.distance_graph.has_node(
                    node_a), f'Node {node_a} not found in the distance graph.'
                assert tree.distance_graph.has_node(
                    node_b), f'Node {node_b} not found in the distance graph.'
                tree.distance_graph.add_edge(node_a, node_b, weight=weight)
                tree.raw_distance_graph.add_edge(node_a, node_b, weight=weight)
            for (node_a, node_b), weight in notused_edges.items():
                assert tree.raw_distance_graph.has_node(
                    node_a), f'Node {node_a} not found in the distance graph.'
                assert tree.raw_distance_graph.has_node(
                    node_b), f'Node {node_b} not found in the distance graph.'
                tree.raw_distance_graph.add_edge(node_a, node_b, weight=weight)
    else:
        logger.warning(
            f'There is no gene alignment available, skipping re-alignment.')
    return pg, tree


def retrieve_gene(G: nx.Graph, pg: Pangenome, tree: Tree):
    G, pg, tree, need_realign_root = retrieve_gene_from_genome(
        G, pg, tree)
    outdir = f"{pg.outdir}/retrieve_gene"
    os.makedirs(outdir, exist_ok=True)
    pg, tree = realign(pg, tree, outdir, need_realign_root)

    tree.update_distance_matrix()

    return G, pg, tree


def retrieve_gene_from_genome(G: nx.Graph, pg: Pangenome, tree: Tree):

    fh_annot = open(f'{pg.outdir}/total.involved_retrieved_annot.tsv', 'w')
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
    need_realign_root = set()
    logger.info(
        f'Extract sequences from {len(need_retrieve)} strains for retrieval.')
    for strain in need_retrieve:
        logger.debug(
            f'Retrieving {len(need_retrieve[strain])} genes for strain {strain}...')
        bar.update()
        strain_seq_dict = {}
        dir_index = int(strain)//1000
        query_fa = f'{pg.outdir}/genome_index/{dir_index}/{strain}/query.fa'
        ref_fa = f'{pg.outdir}/genome_index/{dir_index}/{strain}/ref.fa'
        ref_mpi = f'{pg.outdir}/genome_index/{dir_index}/{strain}/ref.mpi'
        out_paf = f'{pg.outdir}/genome_index/{dir_index}/{strain}/out.paf'
        # the retrieved node want to merge to this node
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
                f'{sfw.seqtk} subseq {pg.prot_file} {temp_file.name} >{query_fa}')
        logger.debug(f'Running and parsing miniprot for strain {strain}...')
        run_command(
            f'{sfw.miniprot} -t {pg.threads} -N 1 -S {ref_mpi} {query_fa} >{out_paf}')
        need_retrieve_node = retrieve_from_paf(fpaf=out_paf)

        for contig in need_retrieve_node:
            retrieved_node = set()
            gene_rank_plus = pg.gene_rank[f'{strain}:{contig}'][0]
            gene_rank_minus = pg.gene_rank[f'{strain}:{contig}'][1]
            gene_rank_all = merge_sorted_lists(gene_rank_plus, gene_rank_minus)
            gene_num = len(gene_rank_plus) + len(gene_rank_minus)
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
                strand = loc['strand']

                if identity < pg.para_id:
                    continue

                if gene_num <= 1:
                    continue

                # check if the gene has been retrieved
                gene_rank = gene_rank_plus if strand == 1 else gene_rank_minus
                if not gene_rank:
                    logger.warning(
                        f'PGAP2 cannot retrieve genes from {strain}:{contig} because no useful synteny information is available.')
                    continue
                if (start in gene_rank) or (end in gene_rank):
                    logger.trace(
                        f'Retrieved gene ({loc}) overlaps existing CDS, skip.')
                    continue
                else:
                    before_s, after_s = find_nearest_numbers(
                        gene_rank, start)
                    before_s2, after_s2 = find_nearest_numbers(gene_rank, end)

                    if (before_s, after_s) != (before_s2, after_s2):
                        logger.trace(
                            f'Retrieved gene ({loc}) overlaps existing CDS, skip.')
                        continue
                before_e = before_s + \
                    pg.annot[f'{strain}:{contig}:{gene_rank_all.index(before_s)}']['len']*3+3
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
                                        'length': 1,  # never used retrieved gene to be the representative
                                        'members': [assumed_gene],
                                        'strains': set([strain]), 'has_para': False,
                                        'repre_nodes': [assumed_gene],
                                    }
                                 ))

                if identity > pg.orth_id and strain not in G.nodes[expect_node]['strains']:
                    G = merge_node(
                        G, pg, tree, [expect_node, assumed_gene], target=expect_node)
                    tree.update_removed_nodes(assumed_gene)

                i += 1
                retrieve_gene_name = f'retrieve_{i}'
                retrieved_node.add(expect_node)
                if '_' in expect_node:
                    expect_node_name = pg.annot[expect_node.split('_')[
                        1]]['id']
                else:
                    expect_node_name = pg.annot[expect_node]['id']
                if not strain_seq_dict:
                    strain_seq_dict = SeqIO.to_dict(
                        SeqIO.parse(ref_fa, 'fasta'))
                feature = SeqFeature(FeatureLocation(
                    start=start, end=end, strand=loc['strand']), type="CDS", id=retrieve_gene_name)
                seq = feature.extract(strain_seq_dict[contig].seq)
                with warnings.catch_warnings(record=True) as w:
                    # only capture BiopythonWarning
                    warnings.simplefilter("always", BiopythonWarning)

                    prot = seq.translate()

                    # check if any warnings were captured
                    if len(w) > 0 and issubclass(w[0].category, BiopythonWarning):
                        pseudo_gene_product = f'pseudo|{identity}|{cigar}'
                        tree.pseudo_root_leaf[tree.leaf_root[expect_node]].add(
                            assumed_gene)
                    else:
                        pseudo_gene_product = f'coding|{identity}|{cigar}'
                        pg.nucl_fa[assumed_gene] = seq
                        need_realign_root.add(tree.leaf_root[expect_node])
                        tree.root_leaf[tree.leaf_root[expect_node]].add(
                            assumed_gene)

                contig_name = pg.annot_contig_map[f'{strain}:{contig}']
                strain_name = pg.strain_dict[strain].strain_name
                location = f"[{start+1}:{end}]({'+' if strand == 1 else '-'})"
                fh_annot.write(
                    f'#{assumed_gene}\t{strain_name}\t{contig_name}\t{location}\t{len(prot)}\t{retrieve_gene_name}\t{expect_node_name}\t{pseudo_gene_product}\t{seq}\t{prot}\n')
                pg.annot[assumed_gene] = {'id': retrieve_gene_name, 'name': expect_node_name, 'product': pseudo_gene_product,
                                          'len': len(prot)}
                logger.debug(
                    f'Retrieved {retrieve_gene_name} in {loc} with identity {identity}, and merged to {expect_node_name}.')
                tree.distance_graph.add_node(assumed_gene)
                tree.distance_graph.add_edge(
                    assumed_gene, expect_node, weight=identity)
                tree.raw_distance_graph.add_node(assumed_gene)
                tree.raw_distance_graph.add_edge(
                    assumed_gene, expect_node, weight=identity)
                tree.leaf_member[assumed_gene] = set([assumed_gene])
                tree.member_leaf[assumed_gene] = assumed_gene
                tree.leaf_root[assumed_gene] = tree.leaf_root[expect_node]
                tree.leaf_member_strains[assumed_gene] = set([strain])

    fh_annot.close()
    logger.info(
        f'A total of {i} genes were found, involving {len(need_realign_root)} clusters that need to be realigned.')
    return G, pg, tree, need_realign_root
