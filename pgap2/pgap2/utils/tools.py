import edlib
import argparse
import itertools
import networkx as nx
import numpy as np
import pandas as pd

from math import ceil
from Bio.Seq import Seq
from loguru import logger
from collections import defaultdict

from pgap2.lib.tree import Tree
from pgap2.lib.pangenome import Pangenome
from pgap2.utils.supply import run_command


def check_min_falen(value):
    ivalue = int(value)
    if ivalue < 11:
        raise argparse.ArgumentTypeError(
            f"Minimum value for --min_falen is 11, but got {ivalue}")
    return ivalue


def is_numeric_pd(pav: pd.DataFrame):
    for value in pav.values.flatten():
        if pd.isna(value):
            continue
        try:
            pd.to_numeric(value)
            return True
        except ValueError:
            return False


def detect_separator(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        possible_seps = [',', '\t']
        sep_counts = {sep: lines[0].count(
            sep) for sep in possible_seps}
        separator = max(sep_counts, key=sep_counts.get)
        if separator == '\t':
            logger.info(f'The separator of the file is [table].')
        elif separator == ',':
            logger.info(f'The separator of the file is [comma].')
        if sep_counts[separator] == 0:
            logger.error('The separator of the file is not found.')
            raise ValueError('The separator of the file is not found.')
    return separator


def check_gcode(value):
    VALID_GCODES_MAP = {
        1: "The Standard Code",
        2: "The Vertebrate Mitochondrial Code",
        3: "The Yeast Mitochondrial Code",
        4: "The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code",
        5: "The Invertebrate Mitochondrial Code",
        6: "The Ciliate, Dasycladacean and Hexamita Nuclear Code",
        9: "The Echinoderm and Flatworm Mitochondrial Code",
        10: "The Euplotid Nuclear Code",
        11: "The Bacterial, Archaeal and Plant Plastid Code",
        12: "The Alternative Yeast Nuclear Code",
        13: "The Ascidian Mitochondrial Code",
        14: "The Alternative Flatworm Mitochondrial Code",
        15: "Blepharisma Nuclear Code",
        16: "Chlorophycean Mitochondrial Code",
        21: "Trematode Mitochondrial Code",
        22: "Scenedesmus obliquus Mitochondrial Code",
        23: "Thraustochytrium Mitochondrial Code",
        24: "Rhabdopleuridae Mitochondrial Code",
        25: "Candidate Division SR1 and Gracilibacteria Code",
        26: "Pachysolen tannophilus Nuclear Code",
        27: "Karyorelict Nuclear Code",
        28: "Condylostoma Nuclear Code",
        29: "Mesodinium Nuclear Code",
        30: "Peritrich Nuclear Code",
        31: "Blastocrithidia Nuclear Code",
        32: "Balanophoraceae Plastid Code",
        33: "Cephalodiscidae Mitochondrial UAA-Tyr Code"
    }

    ivalue = int(value)

    if ivalue not in VALID_GCODES_MAP:
        valid_gcodes_str = '\n'.join(
            f"\t\t{code}: {name}" for code, name in VALID_GCODES_MAP.items())
        raise argparse.ArgumentTypeError(
            f"Invalid genetic code {ivalue}. Valid codes are:\n{valid_gcodes_str}"
        )
    return ivalue


def bbh_check(pg: Pangenome, tree: Tree, para_clust_a, para_clust_b):
    need_merge = False
    clust_a_mem = tree.leaf_member[para_clust_a]
    clust_b_mem = tree.leaf_member[para_clust_b]
    union_mem = clust_a_mem.union(clust_b_mem)
    clust_a_strains = tree.leaf_member_strains[para_clust_a]
    clust_b_strains = tree.leaf_member_strains[para_clust_b]
    overlap_strains = clust_a_strains.intersection(clust_b_strains)

    overlap_repre_strain_dict = defaultdict(list)
    for each_gene in union_mem:
        if (strain_index := int(each_gene.split(':')[0])) in overlap_strains:
            overlap_repre_strain_dict[strain_index].append(
                each_gene)

    best_para_id = 0
    best_overlap_pair = []
    # stop_event = mp.Event()
    overlap_repre_gen = (pair for value_list in overlap_repre_strain_dict.values(
    ) for pair in itertools.combinations(value_list, 2))
    for a, b in overlap_repre_gen:
        a, b, this_pwid = get_identity_with_name(
            a, b, pg.nucl_fa[a], pg.nucl_fa[b])
        if this_pwid >= pg.dup_id:
            need_merge = True
            best_para_id = this_pwid
            break
        else:
            if this_pwid > best_para_id:
                best_para_id = this_pwid
                best_overlap_pair = [a, b]
    if best_para_id < pg.orth_id:
        need_merge = False
    else:
        if not need_merge:
            bbh_dict = {key: [item for item in (
                clust_b_mem if key in clust_b_mem
                else clust_a_mem) if item != key]
                for key in best_overlap_pair}
            bbh_dict_gen = (
                (key, value) for key, value_list in bbh_dict.items() for value in value_list)
            event_flag = False
            for a, b in bbh_dict_gen:
                this_pwid = get_identity(pg.nucl_fa[a], pg.nucl_fa[b])
                if this_pwid > best_para_id:
                    event_flag = True
                    break
            if event_flag:
                need_merge = False
            else:
                need_merge = True

    return need_merge


def calculate_pwid(sA, sB, nucl=True):
    if nucl:
        additional_equalities = [
            ('A', 'N'), ('C', 'N'), ('G', 'N'), ('T', 'N')]
    else:
        additional_equalities = [('*', 'X'), ('A', 'X'), ('C', 'X'), ('B', 'X'),
                                 ('E', 'X'), ('D', 'X'), ('G', 'X'), ('F', 'X'),
                                 ('I', 'X'), ('H', 'X'), ('K', 'X'), ('M', 'X'),
                                 ('L', 'X'), ('N', 'X'), ('Q', 'X'), ('P', 'X'),
                                 ('S', 'X'), ('R', 'X'), ('T', 'X'), ('W', 'X'),
                                 ('V', 'X'), ('Y', 'X'), ('X', 'X'), ('Z', 'X'),
                                 ('D', 'B'), ('N', 'B'), ('E', 'Z'), ('Q', 'Z')]
    aln = edlib.align(sA, sB, mode="NW", task='distance', k=0.5 *
                      len(sA), additionalEqualities=additional_equalities)
    if aln['editDistance'] == -1:
        return 0.0
    return 1.0 - aln['editDistance'] / float(len(sB))


def get_identity(seqA, seqB, nucl=True):
    if len(seqA) > len(seqB):
        seqA, seqB = seqB, seqA
    if nucl:
        pwid = max(calculate_pwid(sA, seqB, nucl=True)
                   for sA in [seqA, str(Seq(seqA).reverse_complement())])
    else:
        pwid = calculate_pwid(seqA, seqB, nucl=False)
    return pwid


def get_identity_with_name(conA, conB, seqA, seqB, nucl=True):
    if len(seqA) > len(seqB):
        seqA, seqB = seqB, seqA
    if nucl:
        pwid = max(calculate_pwid(sA, seqB, nucl=True)
                   for sA in [seqA, str(Seq(seqA).reverse_complement())])
    else:
        pwid = calculate_pwid(seqA, seqB, nucl=False)
    return ((conA, conB, pwid))


def get_similarity(conA, conB):
    if 0 in (len(conA), len(conB)):
        return 0
    similarity = len(set(conA).intersection(
        set(conB)))/min(len(conA), len(conB))
    return similarity


def run_mmseq2(data, data_type: str, id: float, coverage: float, outdir: str, threads: int = 8, verbose: bool = False):
    '''
    Deprecated method
    '''

    if data_type == "fasta":
        fa = data
        run_command(
            f'mmseqs createdb --shuffle 1 {fa} {outdir}/seq.db -v 0')
        data_index = f'{outdir}/seq.db'
    elif data_type == "index":
        data_index = data
    else:
        raise (f'ERROR. Tell me on github')

    run_command(
        f'mmseqs linclust {data_index} {outdir}/seq.clst {outdir}/tmp --min-seq-id {id} -c 0 --threads {threads} --kmer-per-seq 20 -v 0')
    run_command(
        f'mmseqs createtsv --first-seq-as-repr 1 {data_index} {data_index} {outdir}/seq.clst {outdir}/this_clust.tab')
    run_command(
        f'mmseqs createsubdb {outdir}/seq.clst {data_index} {outdir}/seq.clst.rep')

    mydict = {}
    with open(f'{outdir}/this_clust.tab') as fh:
        for line in fh:
            group, header = line.strip().split('\t')
            if group not in mydict:
                mydict[group] = []
            mydict[group].append(header)
    return mydict, f'{outdir}/seq.clst.rep'


def find_mci(tree: Tree, G: nx.Graph, node):
    this_mci = None
    try:
        this_mci = G.nodes[node]['mci']
    except:
        if tree.has_node(node):
            this_mci = tree.nodes[node]['mci']
    assert this_mci is not None, logger.error(f'No MCI for {node}')
    return this_mci


def insert_node(G, before, after, node):
    assert G.has_node(before), f"Error, could not find {before} in G"
    assert G.has_node(after), f"Error, could not find {after} in G"
    G.add_nodes_from([node])
    G.add_edges_from([(before, node[0]), (after, node[0])])
    return G


def find_duplicates(input_list):
    # Find duplicates in a list
    seen = set()
    duplicates = set()

    for item in input_list:
        if item in seen:
            duplicates.add(item)
        else:
            seen.add(item)

    return set(duplicates)


def gen_node_iterables(G, nodes, attr):
    for n in nodes:
        yield G.nodes[n][attr]


def merge_value_from_nodes(G, nodes, attr, dedup=True):
    value = set()
    for f in itertools.chain.from_iterable(gen_node_iterables(G, nodes, attr)):
        value.add(f)
    return value


def test_expect_identity(tree: Tree, G: nx.Graph, u, v, identity):

    max_u = G.nodes[u].get('max_id', None)
    max_v = G.nodes[v].get('max_id', None)
    need_merge = None
    if max_u is None and max_v is None:
        G.nodes[u]['max_id'] = identity
        G.nodes[v]['max_id'] = identity
        need_merge = True
    elif max_u is None:
        if abs(identity - max_v) <= tree.expect_identity:
            G.nodes[u]['max_id'] = max(max_v, identity)
            need_merge = True
        else:
            need_merge = False
    elif max_v is None:
        if abs(identity - max_u) <= tree.expect_identity:
            G.nodes[v]['max_id'] = max(max_u, identity)
            need_merge = True
        else:
            need_merge = False
    else:
        max_id = max(max_u, max_v)
        if abs(identity - max_id) <= tree.expect_identity:
            G.nodes[u]['max_id'] = max_id
            G.nodes[v]['max_id'] = max_id
            need_merge = True
        else:
            need_merge = False
    if not need_merge:
        logger.debug(
            '[Diversity] reject merge {} and {} -> {}/{}/{}'.format(u, v, identity, max_u, max_v))
    return need_merge


def merge_node(G: nx.Graph, pg, tree: Tree, sources: list, target: str):
    '''
    merge node, from source to target, and its attributions

    Parameters:
    G: network
    pg: Pangenome [used to judge the length of source and target node]
    tree: In order to get the corresponding uni and mci in tree
    sources: node that need to merge
    target: target nood

    Return:
    G
    '''

    sources = set(sources)
    nodes = sources.copy()
    sources.remove(target)

    if tree is None:
        G.nodes[target]['repre_nodes'] = [target]
    else:
        G.nodes[target]['repre_nodes'] = merge_value_from_nodes(
            G, nodes, 'repre_nodes')

    members = merge_value_from_nodes(G, nodes, 'members')
    strains = merge_value_from_nodes(G, nodes, 'strains')

    G.nodes[target]['members'] = members
    G.nodes[target]['strains'] = set(strains)
    G.nodes[target]['has_para'] = True if len(
        members) != len(strains) else False

    edges = []
    for source in sources:
        for neighbour in list(G.neighbors(source)):
            if neighbour in nodes:
                continue
            if not G.has_edge(target, neighbour):
                edges.append((target, neighbour))
    G.remove_nodes_from(sources)
    G.add_edges_from(edges)
    return G


def shortest_path_length_with_max_length(G, source, target, source_record, depth_limit=20):
    """Bidirectional shortest path helper with depth limit.

    Returns the shortest path and its length within the specified depth limit.
    """
    if target == source:
        # Source and target are the same node, return the path with one node
        return source_record, [source]

    Gpred = G.adj  # For undirected graph, both Gpred and Gsucc are G.adj
    Gsucc = G.adj

    # Predecessors and successors in search
    pred = {source: None}  # Track predecessors for the forward search
    succ = {target: None}  # Track successors for the reverse search

    # Initialize fringes
    forward_fringe = [source]
    reverse_fringe = [target]
    forward_depth = {source: 0}
    reverse_depth = {target: 0}

    while forward_fringe and reverse_fringe:
        # Expand the smaller fringe first
        if len(forward_fringe) <= len(reverse_fringe):
            this_level = forward_fringe
            forward_fringe = []

            for v in this_level:
                # Check the cumulative depth limit (forward + reverse)
                if forward_depth[v] + min(reverse_depth.values(), default=0) + 1 > depth_limit:
                    return source_record, []  # Return empty path if cumulative depth limit exceeded

                # Check if v has already been cached
                if v in source_record:
                    source_adjs = source_record[v]
                else:
                    # Cache the neighbors of v for future use
                    source_adjs = set(Gsucc[v])  # Store as a set directly
                    source_record[v] = source_adjs

                for w in source_adjs:
                    if w not in pred:
                        forward_fringe.append(w)
                        pred[w] = v  # Track the predecessor
                        forward_depth[w] = forward_depth[v] + 1

                    # Check if the forward and reverse fringes meet
                    if w in succ:
                        # Build the path when forward and reverse meet
                        forward_path = []
                        reverse_path = []

                        # Trace the path from source to w
                        current = w
                        while current is not None:
                            forward_path.append(current)
                            current = pred[current]
                        forward_path.reverse()  # We built it backwards, so reverse it

                        # Trace the path from w to target
                        current = succ[w]
                        while current is not None:
                            reverse_path.append(current)
                            current = succ[current]

                        # Combine forward and reverse paths
                        return source_record, forward_path + reverse_path

        else:
            this_level = reverse_fringe
            reverse_fringe = []

            for v in this_level:
                # Check the cumulative depth limit (forward + reverse)
                if reverse_depth[v] + min(forward_depth.values(), default=0) + 1 > depth_limit:
                    return source_record, []  # Return empty path if cumulative depth limit exceeded

                # Check if v has already been cached
                if v in source_record:
                    source_adjs = source_record[v]
                else:
                    # Cache the neighbors of v for future use
                    source_adjs = set(Gpred[v])  # Store as a set directly
                    source_record[v] = source_adjs

                for w in source_adjs:
                    if w not in succ:
                        reverse_fringe.append(w)
                        succ[w] = v  # Track the successor
                        reverse_depth[w] = reverse_depth[v] + 1

                    # Check if the forward and reverse fringes meet
                    if w in pred:
                        # Build the path when forward and reverse meet
                        forward_path = []
                        reverse_path = []

                        # Trace the path from source to w
                        current = w
                        while current is not None:
                            forward_path.append(current)
                            current = pred[current]
                        forward_path.reverse()  # We built it backwards, so reverse it

                        # Trace the path from w to target
                        current = succ[w]
                        while current is not None:
                            reverse_path.append(current)
                            current = succ[current]

                        # Combine forward and reverse paths
                        return source_record, forward_path + reverse_path

    return source_record, []  # Return empty list if no path found


def test_connectedness(tree: Tree, G: nx.Graph, u, v, sensitivity='moderate'):
    flag = False
    if sensitivity == 'soft':
        flag = True
    else:
        u_repre_nodes = list(G.nodes[u]['repre_nodes'])
        v_repre_nodes = list(G.nodes[v]['repre_nodes'])
        len_u = len(u_repre_nodes)
        len_v = len(v_repre_nodes)

        if max(len_u, len_v) == 1:
            if tree.distance_graph.has_edge(u_repre_nodes[0], v_repre_nodes[0]):
                flag = True
            else:
                flag = False
        # elif sensitivity == 'strict' or max(len_u, len_v) == 2:
        elif sensitivity == 'strict':
            u_indices = [int(tree.matrix_node_map[node])
                         for node in u_repre_nodes]
            v_indices = [int(tree.matrix_node_map[node])
                         for node in v_repre_nodes]
            submatrix = tree.distance_matrix[u_indices, :][:, v_indices]

            # if np.any(submatrix == 0):
            #     flag = False
            # else:
            #     flag = True

            full_size = len(u_repre_nodes) * len(v_repre_nodes)
            if submatrix.size == full_size:
                flag = True
            else:
                flag = False

        elif sensitivity == 'moderate':
            u_indices = [tree.matrix_node_map[node] for node in u_repre_nodes]
            v_indices = [tree.matrix_node_map[node] for node in v_repre_nodes]
            submatrix = tree.distance_matrix[u_indices, :][:, v_indices]
            needed_a = ceil(len(v_repre_nodes) / 2)
            needed_b = ceil(len(u_repre_nodes) / 2)

            u_connections = submatrix.sum(axis=1)
            v_connections = submatrix.sum(axis=0)

            sufficient_u = np.sum(u_connections >= needed_a)
            sufficient_v = np.sum(v_connections >= needed_b)

            if sufficient_u >= needed_b and sufficient_v >= needed_a:
                flag = True
            else:
                flag = False

        else:
            raise ValueError(
                "Invalid sensitivity level. Choose from 'soft', 'moderate', or 'strict'.")
    if not flag:
        logger.debug(
            '[Connectedness] reject merge {} and {} under the {} sensitivity mode'.format(u, v, sensitivity))
    return flag


def merge_judge(tree: Tree, G: nx.Graph, pg: Pangenome, u, v, identity, context_sim=0, flank=5, sensitivity='moderate'):
    need_merge = False
    need_merge = test_connectedness(tree, G, u, v, sensitivity)
    if need_merge:
        if G.nodes[u]['strains'].intersection(G.nodes[v]['strains']):
            need_merge = test_paralog_bbh(
                G, pg, tree, u, v, context_sim, flank)
    if need_merge:
        # must be the last step of merge judge, because it will change the graph attr max_id
        need_merge = test_expect_identity(tree, G, u, v, identity)
    return need_merge


def get_orth_mci(G: nx.Graph, tree: Tree, clust_a, a):
    lca_inner_a = 0.9
    if len(G.nodes[clust_a]['repre_nodes']) == 1:
        if len(G.nodes[clust_a]['members']) > 1:
            level = set([tree.ortho_para[member]
                        for member in G.nodes[clust_a]['members']])
            if len(level) == 1:
                lca_inner_a = tree.dup_id
            elif len(level) > 1:
                lca_inner_a = tree.orth_id
    else:
        for each_repre in G.nodes[clust_a]['repre_nodes']:
            if each_repre == a:
                continue
            if not tree.leaf_member_strains[each_repre] & tree.leaf_member_strains[a]:
                if tree.distance_graph.has_edge(each_repre, a):
                    this_lca_inner = tree.distance_graph[each_repre][a]['weight']
                    if this_lca_inner > lca_inner_a:
                        lca_inner_a = this_lca_inner
    return lca_inner_a


def test_paralog_bbh(G: nx.Graph, pg: Pangenome, tree: Tree, clust_a, clust_b, context_sim, flank):

    need_merge = False
    if not pg.exhaust_orth:
        for a, b in itertools.product(G.nodes[clust_a]['repre_nodes'], G.nodes[clust_b]['repre_nodes']):
            need_merge = False
            if tree.leaf_member_strains[a] & tree.leaf_member_strains[b]:
                para_id = None
                if tree.distance_graph.has_edge(a, b):
                    para_id = tree.distance_graph[a][b]['weight']
                if para_id is None:
                    continue
                # if para_id >= pg.dup_id:
                #     need_merge = True
                # elif para_id >= get_orth_mci(G, tree, clust_a, a) or para_id >= get_orth_mci(G, tree, clust_b, b):
                #     need_merge = True
                clust_a_orth_mci = get_orth_mci(G, tree, clust_a, a)
                clust_b_orth_mci = get_orth_mci(G, tree, clust_b, b)
                if para_id >= clust_a_orth_mci or para_id >= clust_b_orth_mci:
                    need_merge = True
                if pg.accurate and need_merge:
                    if (a, b) in tree.pan_judge_dict:
                        need_merge = tree.pan_judge_dict[(
                            a, b)]
                    elif (b, a) in tree.pan_judge_dict:
                        need_merge = tree.pan_judge_dict[(
                            b, a)]
                    else:
                        need_merge = bbh_check(
                            pg, tree, a, b)
                        tree.pan_judge_dict[(
                            a, b)] = need_merge
                if not need_merge:
                    logger.debug(
                        '[BBH] reject merge {} and {} -> {}/{}/{}'.format(clust_a, clust_b, para_id, clust_a_orth_mci, clust_b_orth_mci))
                    break

    return need_merge
