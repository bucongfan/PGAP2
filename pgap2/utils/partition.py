import os
import pickle
import argparse
import itertools
import numpy as np
import networkx as nx

from tqdm import tqdm
from loguru import logger
from collections import Counter, defaultdict
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components as ccomps

from argparse import ArgumentParser, _SubParsersAction

from pgap2.lib.tree import Tree
from pgap2.lib.basic import Basic
from pgap2.lib.pklcheck import PklCheck
from pgap2.lib.panclust import Panclust
from pgap2.lib.pangenome import Pangenome

from pgap2.utils.supply import sfw, tqdm_
from pgap2.utils.supply import set_verbosity_level, run_command
from pgap2.utils.generate_tree import generate_tree
from pgap2.utils.gene_retriever import retrieve_gene
from pgap2.utils.arrangement_detector import merge_by_synteny
from pgap2.utils.data_loader import file_parser, get_file_dict
from pgap2.utils.tools import merge_node, shortest_path_length_with_max_length, get_similarity, merge_judge, check_min_falen, check_gcode, find_final_node

"""
main function for partitioning pangenome data into clusters based on identity and synteny.
This module provides functions to generate a network of genes, classify paralogs, and select representative nodes based on gene identity and synteny.

input:
- G: Gene network graph, edge is gene distance.
- tree: identity tree of the genes.
- pg: Pangenome object containing strain and gene information.

output:
- total.involved_annot.tsv: gene annotation information of each cluster
- total.involved_prot.fa: all protein sequences involved in this analysis
- pgap2.partition.summary_statistics.txt: Pan-group statistic result
- pgap2.partition.gene_content.detail.tsv: Partitioning result with annotation information
- pgap2.partition.gene_content.pav: Presence-Absence Variation Matrix
- pgap2.partition.gene_content.csv: gene content of each cluster
- pgap2.partition.map.gml: Graph of pangenome
- basic.pkl: Binary file recorded necessary parameters and file paths used to quick downstream analysis
- partition.log: Running log

"""


def select_repre_node(G: nx.Graph, tree: Tree, each_gene, para_repre, para_context):
    # Select the representative node for each gene based on the context similarity
    child_context = para_context[each_gene]
    max_similarity = 0
    # max_repre = None
    max_pool = []
    for each_repre in para_repre:
        repre_context = para_context[each_repre]
        similarity = get_similarity(
            child_context, repre_context)
        if similarity > max_similarity:
            max_pool = [each_repre]
            max_similarity = similarity
            # max_repre = each_repre
        elif similarity == max_similarity:
            max_pool.append(each_repre)
        if similarity >= 1:
            max_pool = [each_repre]
            break

    if len(max_pool) == 1:
        return max_pool[0]
    elif len(max_pool) > 1:  # multiple representatives has similar context
        clust_099 = tree.ortho_para[each_gene]
        for each_repre in max_pool:
            clust_099_ref = tree.ortho_para[each_repre]
            if clust_099 == clust_099_ref:
                return each_repre
        # Choose the largest one
        repre_with_max_len = max(
            tmp_len := {k: G.nodes[k]['length'] for k in max_pool}, key=tmp_len.get)
        logger.trace(
            f'Conflict in {each_gene} with {max_pool}, choose {repre_with_max_len}')
        return repre_with_max_len
    elif len(max_pool) == 0:
        return None


def classify_paralogs(members):
    class_dict = defaultdict(list)
    for member in members:
        strain_index, _ = member.split(':', 1)
        class_dict[strain_index].append(member)
    para_clusts = []
    other_clusts = []
    for clust in class_dict.values():
        if len(clust) > 1:
            para_clusts.extend(clust)
        else:
            other_clusts.extend(clust)
    return [para_clusts], other_clusts


def get_paralogs_repre(members):
    # Find a paralogous gene that has the most occurrences in the members
    most_common_elements = Counter([member.split(':')[0]
                                    for member in members]).most_common()
    # Filter out
    filtered_elements = [
        element for element in most_common_elements if element[1] >= 2]

    if not filtered_elements:
        return None

    max_count = filtered_elements[0][1]

    # Find the gene that meets the most occurrences
    for element in filtered_elements:
        if element[1] == max_count:
            paralog_repre = defaultdict(list)
            this_count = 0
            for each_gene in members:
                strain_index, _ = each_gene.split(':', 1)
                if strain_index == element[0]:
                    paralog_repre[each_gene] = [each_gene]
                    this_count += 1
                if this_count == max_count:
                    return paralog_repre


def generate_network(pg: Pangenome, tree: Tree):
    '''
    Put the index of the gene into the network,
    If the gene is a paralog, it will be split to the single node.
    Then merge nodes according to the synteny information

    pg: Pangenome object
    tree: nx.DiGraph from trimmed N-rooted fusion tree from iterative cd-hit results

    return:
        G: nx.Graph link all the nodes with its real relative distance in the genome
        tree: nx.DiGraph update the node if the node is a paralog
    '''

    G = nx.Graph()
    para_dict = defaultdict(list)
    nodes = []
    edges = []

    member2node = {}
    orth_identity_tree: nx.DiGraph = tree.orth_identity_tree
    split_result_map = defaultdict(list)

    logger.info(
        f'---- Adding high confidence nodes as the network\'s backbone...')
    for node in tqdm(tree.leaf_root.keys(), unit=' node', desc=tqdm_.step(3), disable=pg.disable_tqdm):
        has_para = orth_identity_tree.nodes[node]['has_para']
        if has_para:
            para_dict[node] = orth_identity_tree.nodes[node]['members']
            for each_gene in para_dict[node]:
                falen = pg.falen[each_gene]
                nodes.append((each_gene, {'length': falen, 'members': set([
                    each_gene]), 'strains': set([int(each_gene.split(':')[0])]), 'has_para': False, 'repre_nodes': [each_gene]}))
                member2node.update({each_gene: each_gene})
        else:
            member2node.update(
                {_: node for _ in orth_identity_tree.nodes[node]['members']})
            falen = pg.falen[node]
            nodes.append(
                (node, {'length': falen, 'members': orth_identity_tree.nodes[node]['members'], 'strains': orth_identity_tree.nodes[node]['strains'], 'has_para': False, 'repre_nodes': [node]}))
            split_result_map[node].append(node)
    G.add_nodes_from(nodes)
    del pg.falen
    logger.info(f'---- Connecting adjacent nodes as the backbone edges...')
    for strain_index in tqdm(pg.strain_dict.keys(), unit=' edge', desc=tqdm_.step(3), disable=pg.disable_tqdm):
        gene_num = pg.strain_dict[strain_index].gene_num
        for contig_index, gene_count in enumerate(gene_num):
            prev_gene = f'{strain_index}:{contig_index}:0'
            if gene_count > 0:
                repre_prev_gene = member2node[prev_gene]
                for gene_index in range(1, gene_count):
                    curr_gene = f'{strain_index}:{contig_index}:{gene_index}'
                    repre_curr_gene = member2node[curr_gene]
                    edges.append((repre_prev_gene, repre_curr_gene))
                    repre_prev_gene = repre_curr_gene
    G.add_edges_from(edges)

    logger.info(
        f'---- Attempting to split potential paralogous genes based on synteny...')
    relabel_dict = {}

    for repre_node, clusts in tqdm(para_dict.items(), unit=' paralog', desc=tqdm_.step(3), disable=pg.disable_tqdm):

        # seperate paralogs and other clusters with signle strain
        para_clusts, other_clusts = classify_paralogs(clusts)
        para_context = {}
        for each_clust in clusts:
            para_context[each_clust] = tree.get_context(each_clust, flank=10)
        logger.debug(
            f'---- Splitting the paralogous clusts of {repre_node} with {len(clusts)} nodes...')
        while True:
            split_clusts = para_clusts
            para_clusts = []
            has_para = False
            for para_clust in split_clusts:
                # get the most strain representative
                para_repre = get_paralogs_repre(para_clust)
                if para_repre:
                    has_para = True
                else:
                    # until all paralogs are processed
                    para_clusts.append(para_clust)
                    continue
                for each_gene in para_clust:
                    if each_gene in para_repre:
                        continue
                    max_repre = select_repre_node(
                        G, tree, each_gene, para_repre, para_context)

                    if not max_repre:
                        logger.debug(
                            f'No repre for {each_gene}, set as itself')
                        para_repre[each_gene] = [each_gene]
                    else:
                        para_repre[max_repre].append(each_gene)
                para_clusts.extend(para_repre.values())
            if not has_para:  # after all paralogs are processed, try to assign the other clusters to the splitted paralogous clusters
                para_repre = defaultdict(list)
                para_repre_map = {}
                for each_clust in para_clusts:
                    repre = each_clust[0]
                    para_repre[repre] = each_clust
                    for each_gene in each_clust:
                        para_repre_map[each_gene] = repre

                for each_other_clust in other_clusts:
                    max_repre = select_repre_node(
                        G, tree, each_other_clust, para_repre_map.keys(), para_context)
                    if not max_repre:
                        logger.debug(
                            f'No repre for {each_other_clust}, set as itself')
                        para_repre[each_other_clust] = [each_other_clust]
                        para_repre_map[each_other_clust] = each_other_clust
                    else:
                        max_repre = para_repre_map[max_repre]
                        para_repre[max_repre].append(each_other_clust)
                para_clusts = para_repre.values()
                break
        logger.debug(
            f'{len(split_clusts)} paralogs split into {len(para_clusts)} clusters')
        for each_clust in para_clusts:
            if len(each_clust) > 1:
                target_node = max(each_clust_len := {
                                  k: G.nodes[k]['length'] for k in each_clust}, key=each_clust_len.get)
                G = merge_node(G, pg, None, sources=each_clust,
                               target=target_node)
            else:
                # node member (through clustering) that cannot be merged through synteny but they always have very high identity
                target_node = each_clust[0]
            relabel_dict.update({target_node: f'{repre_node}_{target_node}'})
            G.nodes[target_node]['repre_nodes'] = [
                f'{repre_node}_{target_node}']
            split_result_map[repre_node].append(f'{repre_node}_{target_node}')

    nx.relabel_nodes(G, relabel_dict, copy=False)
    tree.load_split_result_map(split_result_map)

    logger.info(f'---- Updating the paralogous nodes...')
    update_nodes = defaultdict(list)
    for node in G.nodes():
        if '_' in node:
            father, child = node.split('_')
            update_nodes[father].append(node)
    for father, childs in update_nodes.items():
        root = tree.leaf_root[father]
        del tree.leaf_root[father]
        for child in childs:
            tree.leaf_root[child] = root
            for member in G.nodes[child]['members']:
                tree.member_leaf[member] = child

    root_leaf = defaultdict(set)
    for node in tree.leaf_root:
        root_leaf[tree.leaf_root[node]].add(node)
    tree.root_leaf = root_leaf

    leaf_member = defaultdict(set)
    leaf_member_strains = defaultdict(set)
    for member in tree.member_leaf:
        leaf_member[tree.member_leaf[member]].add(member)
        leaf_member_strains[tree.member_leaf[member]].add(
            int(member.split(':')[0]))
    tree.leaf_member = leaf_member
    tree.leaf_member_strains = leaf_member_strains

    return G, tree


def find_para_node(G, node):
    if G.nodes[node]['has_para']:
        return node
    else:
        predecessors = list(G.predecessors(node))
        if len(predecessors) > 0:
            return find_para_node(G, predecessors[0])
        else:
            return None


def get_pan_clust(G: nx.Graph, pg: Pangenome, tree: Tree, clust):

    clust_nodes = G.nodes[clust]
    nodes = defaultdict(set)
    repre_node = None
    repre_node_len = 0

    for node in clust_nodes['repre_nodes']:
        child_node = node.split('_')[1] if '_' in node else node
        members = tree.leaf_member[node]
        nodes[node].update(members)
        annot = pg.annot[child_node]
        if annot['len'] > repre_node_len:
            repre_node = child_node
            repre_node_len = annot['len']

    subgraph = tree.raw_distance_graph.subgraph(nodes.keys())

    # cache the count of each node
    counts = {node: len(members) for node, members in nodes.items()}

    id_list = np.array([weight * counts[u] * counts[v]
                       for u, v, weight in subgraph.edges(data='weight')])

    minimum_id = np.min(id_list, initial=1) if id_list.size > 0 else 1
    average_id = np.round(np.mean(id_list), 5) if id_list.size > 0 else 1
    var = np.round(np.var(id_list), 5) if id_list.size > 0 else 0

    node_set = set(nodes)

    uni = 0
    for u, v, weight in subgraph.edges(data='weight'):
        if (u in node_set and v not in node_set) or (v in node_set and u not in node_set):
            uni = max(uni, weight)

    my_pan_clust = Panclust(one_pan=clust_nodes['members'], min=minimum_id, repre_node=repre_node,
                            uni=uni, var=var, mean=average_id)
    return my_pan_clust


def similarity_partition(tree: Tree, G: nx.Graph, nodes, search_distance, pre_compute, pre_changed_nodes):
    nodes = list(nodes)
    real_nodes = set()
    node_distances = {}  # record the distance

    for a, b in itertools.combinations(nodes, 2):
        has_changed = False
        if a in pre_compute:
            path = pre_compute[a]
            if set(path).intersection(pre_changed_nodes):
                del pre_compute[a]
                has_changed = True
        else:
            has_changed = True

        if b in pre_compute:
            path = pre_compute[b]
            if set(path).intersection(pre_changed_nodes):
                del pre_compute[b]
                has_changed = True
        else:
            has_changed = True

        if not has_changed:
            continue
        a_adj, path = shortest_path_length_with_max_length(
            G, a, b, {}, search_distance)
        if path:
            for each_node in path:
                if G.has_node(each_node) and G.degree(each_node) > 1:
                    real_nodes.add((a, b))
                    node_distances[(a, b)] = len(path)
                    break
        else:
            pre_compute.update({a: set.union(*a_adj.values())})
    if not real_nodes:
        return [], []
    data = []
    row_ind = []
    col_ind = []
    for (a, b) in real_nodes:
        data.append(1)
        row_ind.append(nodes.index(a))
        col_ind.append(nodes.index(b))
    csr = csr_matrix((data, (row_ind, col_ind)),
                     shape=(len(nodes), len(nodes)))
    num_components, labels = ccomps(csr, directed=False)

    # A dictionary that maps each label to the nodes that belong to it
    components = defaultdict(list)
    for i, label in enumerate(labels):
        components[label].append(nodes[i])

    need_merge_nodes = []
    merge_node_attr = []  # store the similarity and distance of each node pair

    for leaf in components.values():
        if len(leaf) > 1:
            this_nodes = set()
            node_map = {}
            for node in leaf:
                for each_node in G.nodes[node]['repre_nodes']:
                    node_map[each_node] = node
                    this_nodes.add(each_node)
            subgraph = tree.distance_graph.subgraph(this_nodes)
            for sub_components in nx.connected_components(subgraph):
                result = set()
                for component in sub_components:
                    result.add(node_map[component])
                if len(result) > 1:
                    need_merge_nodes.append(list(result))
                    cluster_nodes = list(result)
                    cluster_attr = {}

                    for node1, node2 in itertools.combinations(cluster_nodes, 2):
                        dist = node_distances.get(
                            (node1, node2), node_distances.get((node2, node1), None))
                        if not dist:
                            continue

                        repre_nodes1 = G.nodes[node1]['repre_nodes']
                        repre_nodes2 = G.nodes[node2]['repre_nodes']
                        max_similarity = 0
                        # Traverse all edges between node1 and node2's repre_nodes to find the maximum weight
                        for rn1, rn2 in itertools.product(repre_nodes1, repre_nodes2):
                            if subgraph.has_edge(rn1, rn2):
                                edge_weight = tree.distance_graph[rn1][rn2].get(
                                    'weight', 0)
                                max_similarity = max(
                                    max_similarity, edge_weight)
                                if max_similarity >= tree.dup_id:
                                    break
                        if not max_similarity:
                            continue
                        cluster_attr.update({(node1, node2): (
                            max_similarity, dist)})
                    merge_node_attr.append(cluster_attr)

    return need_merge_nodes, merge_node_attr


def merge_by_similarity(G: nx.Graph, pg: Pangenome, tree: Tree, fast: bool = False, sensitivity: str = 'strict', radius: int = 3, context_sim: float = 0, flank: int = 5, disable: bool = True, step: int = 4):

    search_distance = radius*2+1
    root_leaf = tree.root_leaf

    iter_count = 0
    removed_nodes = tree.get_removed_nodes()

    merge_event = True
    pre_compute = {}
    changed_nodes = set()
    while merge_event:
        iter_count += 1
        merge_event = False
        pre_changed_nodes = changed_nodes
        changed_nodes = set()
        for main_nodes in tqdm(root_leaf.values(), unit=f" Round: {iter_count}", disable=disable, desc=tqdm_.step(step=step)):
            if len(main_nodes) == 1:
                continue
            exists_nodes = set()
            for node in main_nodes:
                if node in removed_nodes:
                    continue
                exists_nodes.add(node)
            if len(exists_nodes) < 2:
                continue

            logger.trace(
                f'Process {len(exists_nodes)} nodes in {main_nodes} with search distance {search_distance}...')
            need_merge_nodes, merge_node_attr = similarity_partition(
                tree, G, exists_nodes, search_distance, pre_compute, pre_changed_nodes)
            logger.trace(f'Found {len(need_merge_nodes)} merge candidates...')

            if not need_merge_nodes:
                continue

            for clusters, cor_attr in zip(need_merge_nodes, merge_node_attr):
                if fast:
                    longest_node = max(
                        clusters, key=lambda x: G.nodes[x]['length'])
                    G = merge_node(G, pg, tree, clusters,
                                   target=longest_node)
                    merge_event = True
                    for v in clusters:
                        if v != longest_node:
                            changed_nodes.update(set([longest_node, v]))
                            removed_nodes.add(v)
                else:
                    split_clust_map = {_: _ for _ in clusters}
                    sorted_edges = sorted(
                        cor_attr.items(), key=lambda x: (-x[1][0], x[1][1]))
                    logger.trace(
                        f'Process {len(sorted_edges)} nodes itertively')

                    changed_result = {}
                    for (u, v), (identity, distance) in sorted_edges:
                        u = find_final_node(u, split_clust_map)
                        v = find_final_node(v, split_clust_map)

                        if u == v:
                            continue
                        u_i = len(G.nodes[u]['repre_nodes'])
                        v_i = len(G.nodes[v]['repre_nodes'])
                        flag = False
                        if (u, v) in changed_result:
                            pre_u_i, pre_v_i, pre_need_merge = changed_result[(
                                u, v)]
                            if u_i == pre_u_i and v_i == pre_v_i:
                                need_merge = pre_need_merge
                                flag = True

                        if not flag:
                            logger.trace(
                                f'[Fine analysis] Checking {u} and {v} with identity {identity}, context_sim {context_sim}, flank {flank}, sensitivity {sensitivity}')
                            # Check if the nodes need
                            need_merge = merge_judge(
                                tree, G, pg, u, v, identity, context_sim, flank, sensitivity)
                            changed_result[(u, v)] = (u_i, v_i, need_merge)
                            changed_result[(v, u)] = (v_i, u_i, need_merge)

                        if need_merge:
                            u, v = (u, v) if G.nodes[u]['length'] > G.nodes[v]['length'] else (
                                v, u)
                            G = merge_node(G, pg, tree, [u, v], target=u)
                            removed_nodes.add(v)
                            split_clust_map[v] = u
                            merge_event = True
                            changed_nodes.update(set([u, v]))
    if pg.retrieve:
        logger.info(f'---- Retrieving genes from the removed nodes...')
        tree.set_removed_nodes(removed_nodes)
    return (G, pg, tree)


def mcl(pg: Pangenome, tree: Tree):
    mcl_result = f'{pg.outdir}/mcl.result'
    run_command(
        f"{sfw.mcxdeblast} -m9 --score r --line-mode=abc {tree.alignment_result} 2> /dev/null | {sfw.mcl} - --abc -I 1.5 -te {pg.threads} -o {mcl_result} &>/dev/null")

    G = nx.Graph()
    G.add_nodes_from(tree.raw_distance_graph.nodes(data=True))
    raw_G = tree.raw_distance_graph
    with open(mcl_result, 'r') as fh:
        for line in fh:
            line = line.rstrip()
            clust = line.split('\t')
            for a, b in itertools.combinations(clust, 2):
                if raw_G.has_edge(a, b):
                    G.add_edge(a, b, weight=raw_G[a][b]['weight'])
    tree.load_distance_graph(G)


def is_complete_graph(G):
    n = len(G.nodes())
    return nx.is_connected(G) and len(G.edges()) == n * (n - 1) / 2


def get_expect_identity(tree: Tree, G: Pangenome, pg: Pangenome):
    all_range_value = []
    for clusts in tree.root_leaf.values():
        strain_all = set()
        need_next = True
        for clust in clusts:
            if strain_all.intersection(tree.leaf_member_strains[clust]):
                need_next = False
                break
            strain_all.update(tree.leaf_member_strains[clust])
        if not need_next:
            continue
        if len(strain_all) >= pg.hconf_count_thre and len(clusts) > 1:
            subgraph = tree.distance_graph.subgraph(clusts)
            if not is_complete_graph(subgraph):
                continue
            weight_set = set()
            need_next = True
            for a, b in itertools.combinations(clusts, 2):
                _, path = shortest_path_length_with_max_length(
                    G, a, b, {})
                if not path:
                    need_next = False
                    break
                if tree.distance_graph.has_edge(a, b):
                    weight_set.add(tree.distance_graph[a][b]['weight'])
            if not need_next:
                continue
            range_value = max(weight_set) - min(weight_set)
            all_range_value.append(range_value)
    if not all_range_value:
        max_in_range = 1 - pg.para_id
        logger.warning(
            f"No valid range values found. returning default value of {max_in_range}")
    else:
        max_in_range = max(all_range_value)
    return round(max_in_range, 5)


def main(indir: str, outdir: str, evalue: float, hconf_thre: float, aligner: str, clust_method: str, falen: int, fast_mode: bool, threads: int, orth_id: float, para_id: float, dup_id: float, id_attr_key: str, type_filter: str, max_targets: int, coverage: float, LD: float, AS: float, AL: float, context_similirity: float, accurate: bool, exhaust_orth: bool, flank: int, disable: bool, annot: bool, gcode: int, retrieve: bool, radius: int, sensitivity: int, ins: bool):

    decode_status = False
    file_dict = get_file_dict(indir)
    if os.path.exists(f'{outdir}/preprocess.pkl'):
        '''
        Found a previous preprocess.pkl file, loading...
        This file contains the parameters and file structure of the previous run.
        If the parameters are not match, it will raise a ValueError.
        If the file structure has changed, it will warn the user and reload the file structure from the current input directory.
        If the file structure is the same, it will load the pangenome and tree from the previous run.
        If the pangenome has invalid genes, it will warn the user and continue to the next step.
        '''
        logger.info(f'Found {outdir}/preprocess.pkl')
        logger.info(f'Loding...')
        with open(f'{outdir}/preprocess.pkl', 'rb') as fh:
            previous: PklCheck = pickle.load(fh)
            logger.info(f'Check the previous file parameters...')
            decode_status = previous.decode(
                orth_id=orth_id, para_id=para_id, dup_id=dup_id, accurate=accurate, coverage=coverage, id_attr_key=id_attr_key, type_filter=type_filter, LD=LD, AS=AS, AL=AL, evalue=evalue, aligner=aligner, clust_method=clust_method, falen=falen, annot=annot, retrieve=retrieve,)
            if decode_status:
                # success
                pg = previous.data_dump('pangenome')
                tree = previous.data_dump('tree')
                previous_file_dict = previous.data_dump('file_dict')

                if previous_file_dict != file_dict:
                    logger.warning(
                        f'File structure has changed')
                    total_name = list(previous_file_dict.keys(
                    )) + [k for k in file_dict.keys() if k not in previous_file_dict]

                    max_width = max(
                        [len(name) for name in total_name+['Previous', 'Current']])+2
                    logger.warning(
                        f'{"Previous":<{max_width}}\t{"Current":<{max_width}}')
                    logger.warning(f'{"-"*(max_width*2)}')
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
                    logger.warning(f'{"-"*(max_width*2)}')
                    len_prev = len(previous_file_dict)
                    len_cur = len(file_dict)
                    logger.warning(
                        f'{len_prev:<{max_width}}\t{len_cur:<{max_width}}')

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

                if decode_status:
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
        '''
        Load strain from input directory
        If the previous preprocess.pkl file is not found or the parameters are not match,
        it will load the strain from the input directory and create a new pangenome object.
        '''
        logger.info(f'Load strain from {indir}')
        pg = file_parser(
            indir=indir, outdir=outdir, annot=annot, threads=threads, disable=disable, retrieve=retrieve, falen=falen, gcode=gcode, id_attr_key=id_attr_key, type_filter=type_filter, prefix='partition')
        file_prot = f'{outdir}/total.involved_prot.fa'
        logger.info(f'Total gene invovled in this project in: {file_prot}')
        pg.load_annot_file(f'{outdir}/total.involved_annot.tsv')
        pg.load_prot_file(file_prot)
        logger.info(
            f'Create distane tree with {pg.strain_num} strains')

        tree = generate_tree(
            input_file=file_prot, orth_list=[dup_id, orth_id], outdir=pg.outdir, evalue=evalue, aligner=aligner, falen=falen, disable=disable, threads=threads, max_targets=max_targets, coverage=coverage, ID=para_id, LD=LD, AS=AS, AL=AL, clust_method=clust_method)

        logger.info(
            f'To save the complete information of this project for breakpoint resume...')
        pickle_preprocess = PklCheck(outdir=outdir, name='preprocess')
        pickle_preprocess.load('file_dict', main_data=file_dict)
        pickle_preprocess.load('pangenome', main_data=pg, parameter={'orth_id': orth_id, 'para_id': para_id, 'dup_id': dup_id, 'accurate': accurate,
                                                                     'id_attr_key': id_attr_key, 'type_filter': type_filter,
                                                                     'coverage': coverage, 'AS': AS, 'AL': AL, 'LD': LD, 'retrieve': retrieve,
                                                                     'evalue': evalue, 'aligner': aligner, 'clust_method': clust_method,
                                                                     'annot': annot, 'falen': falen})
        pickle_preprocess.load('tree', main_data=tree)
        pickle_preprocess.pickle_()

    '''
    load necessary parameters and file paths used to quick downstream analysis
    '''
    pg.orth_id = orth_id
    pg.para_id = para_id
    pg.dup_id = dup_id
    pg.accurate = accurate
    pg.exhaust_orth = exhaust_orth
    pg.retrieve = retrieve
    pg.evalue = evalue
    pg.aligner = aligner
    pg.LD = LD
    pg.AL = AL
    pg.AS = AS
    pg.load_hconf(hconf_thre=hconf_thre)

    # -----------------------------------partition step-----------------------------------#

    logger.info('Get the gene primal clust result by mcl')
    mcl(pg, tree)
    logger.info('Load the gene length information')
    pg.reload_nucl_file(tree)
    logger.info('Create synteny network')
    G, tree = generate_network(pg=pg, tree=tree)

    tree.load_para_id(para_id)
    tree.load_orth_id(orth_id)
    tree.load_dup_id(dup_id)
    logger.info('Build index')
    max_in_range = get_expect_identity(tree, G, pg)

    logger.info(f'Load expect identity: {max_in_range}')
    tree.load_expect_identity(max_in_range)
    logger.info(f'Clean up the distance graph according to paralogous genes')
    tree.update_distance_graph(disable=disable)
    logger.info(f'Merge by gene similarity')

    G, pg, tree = merge_by_similarity(G=G, pg=pg, tree=tree,
                                      sensitivity=sensitivity,
                                      radius=radius, fast=fast_mode,
                                      context_sim=context_similirity,
                                      flank=flank,
                                      disable=disable,)

    logger.info(f'Double check through gene synteny')
    G = merge_by_synteny(G, pg, tree,
                         context_sim=context_similirity,
                         flank=flank,
                         sensitivity=sensitivity,
                         ins=ins,
                         )

    logger.info(f'Reload the gene annotation')
    pg.reload_annot_file(retrieve=retrieve)

    if retrieve:
        logger.info(f'Retrieve gene from missing')
        G, pg, tree = retrieve_gene(G, pg, tree)
        G, pg, tree = merge_by_similarity(G=G, pg=pg, tree=tree,
                                          sensitivity=sensitivity,
                                          radius=radius, fast=fast_mode,
                                          context_sim=context_similirity,
                                          flank=flank,
                                          disable=disable, step=10)
        G = merge_by_synteny(G=G, pg=pg, tree=tree,
                             context_sim=context_similirity,
                             flank=flank,
                             sensitivity=sensitivity,
                             ins=ins, step=11
                             )

    logger.info('Organize the results')
    pg.init_pan_temp()

    bar = tqdm(range(G.number_of_nodes()), unit=f" Organize",
               disable=disable, desc=tqdm_.step(-1))
    for node in G.nodes():
        bar.update()
        my_pan_clust = get_pan_clust(G, pg, tree, node)
        pg.load_one_pan(pan_clust=my_pan_clust)
    bar.close()
    logger.info('Dump the gene content matrix')
    pg.dump_csv(outdir=outdir, prefix='pgap2.partition')
    logger.info('Dump the gene map graph')
    H = nx.Graph()
    H.add_nodes_from(G.nodes())
    H.add_edges_from(G.edges())
    nx.write_gml(H, f"{outdir}/pgap2.partition.map.gml")

    logger.info(
        f'To save the basic results of this project for downstream visulization...')
    pickle_basic = PklCheck(outdir=outdir, name='basic')
    pickle_basic.load('basic', main_data=Basic(pg=pg))
    pickle_basic.pickle_()
    return 0


def launch(args: argparse.Namespace):
    main(indir=os.path.abspath(args.indir), outdir=os.path.abspath(args.outdir),
         falen=args.min_falen, threads=args.threads, evalue=args.evalue,
         aligner=args.aligner, clust_method=args.clust_method,
         orth_id=args.orth_id, para_id=args.para_id, dup_id=args.dup_id,
         id_attr_key=args.id_attr_key, type_filter=args.type_filter,
         coverage=0.98, fast_mode=args.fast_mode, hconf_thre=args.hconf_thre,
         #  coverage=args.coverage,
         LD=args.LD, AS=args.AS, AL=args.AL, max_targets=args.max_targets,
         # notused parameters, set a default value and will discard in the next release if dont use for sure
         context_similirity=0, flank=5,
         accurate=args.accurate,
         exhaust_orth=args.exhaust_orth,
         gcode=args.gcode,
         disable=args.disable, annot=args.annot, retrieve=args.retrieve,
         radius=args.radius, sensitivity=args.sensitivity, ins=args.ins)


def partition_portal(args):
    set_verbosity_level(args.outdir, args.verbose,
                        args.debug, 'partition')

    if args.clust_method == 'mmseqs2':
        sfw.check_dependency("mmseqs2")
    elif args.clust_method == 'cdhit':
        sfw.check_dependency("cdhit")

    sfw.check_dependency("mcl")
    sfw.check_dependency("mcxdeblast")
    if args.aligner == 'diamond':
        sfw.check_dependency("diamond")
    elif args.aligner == 'blastp':
        sfw.check_dependency("blastp")
        sfw.check_dependency("makeblastdb")

    if args.retrieve:
        sfw.check_dependency("miniprot")
        sfw.check_dependency("seqtk")
    if args.annot:
        sfw.check_dependency("prodigal")
    if args.retrieve:
        tqdm_.set_total_step(12)
    else:
        tqdm_.set_total_step(7)
    launch(args)


def partition_cmd(subparser: _SubParsersAction):

    subparser_partition: ArgumentParser = subparser.add_parser(
        'main', help='Core functions of PGAP2', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparser_partition.add_argument(
        '--indir', '-i', required=True, help='Input file contained, same prefix seems as the same strain.',)
    subparser_partition.add_argument(
        '--outdir', '-o', required=False, help='Output directory', default='.',)
    subparser_partition.add_argument('--dup_id', required=False, type=float,
                                     default=0.99, help='The maximum identity between the most recent duplication envent.')
    subparser_partition.add_argument('--orth_id', required=False, type=float,
                                     default=0.98, help='The maximum identity between the most similar panclusters.')
    subparser_partition.add_argument('--para_id', required=False, type=float,
                                     default=0.7, help='Use this identity as the paralogous identity.')
    subparser_partition.add_argument("--type-filter", required=False, type=str,
                                     default='CDS', help="Only for gff file as input, feature type (3rd column) to include, Only lines matching these types will be processed.")
    subparser_partition.add_argument("--id-attr-key", required=False, type=str,
                                     default='ID', help="Only for gff file as input, Attribute key to extract from the 9th column as the record ID (e.g., 'ID', 'gene', 'locus_tag').")
    subparser_partition.add_argument('--hconf_thre', required=False, type=float,
                                     default=1, help='The threshold to define high confidence cluster which is used to evaluate the cluster diversity. Loose this value when your input is too large or too diverse, such as 0.95.')
    subparser_partition.add_argument('--exhaust_orth', '-e', required=False, action='store_true',
                                     default=False, help='Try to split every paralogs gene exhausted')
    subparser_partition.add_argument('--sensitivity', '-s', required=False, type=str,
                                     default='strict', choices=('soft', 'moderate', 'strict'), help='The degree of connectedness between each node of clust.')
    subparser_partition.add_argument('--ins', '-n', required=False,
                                     action='store_true', default=False, help='Ignore the influence of insertion sequence.')
    subparser_partition.add_argument('--fast', '-f', dest='fast_mode', required=False,
                                     action='store_true', default=False, help='Do not apply fine feature analysis just partition according to the gene identity and synteny.')
    subparser_partition.add_argument('--accurate', '-a', required=False,
                                     action='store_true', default=False, help='Apply bidirection check for paralogous gene partition (useless if exhaust_orth asigned).')
    subparser_partition.add_argument('--retrieve', '-r', required=False,
                                     action='store_true', default=False, help='Retrieve gene that may lost with annotations')
    subparser_partition.add_argument(
        '--threads', '-t', required=False, default=1, help='threads used in parallel', type=int)
    # subparser_partition.add_argument('--coverage', required=False, type=float,
    #                                  default=0.98, help='Length difference for cdhit cluster.')
    # subparser_partition.add_argument('--context_similirity', '-s', required=False, type=float,
    #                                  default=0, help='The context similarity threshold of gene synteny.')
    # subparser_partition.add_argument('--flank', '-l', required=False, type=int,
    #                                  default=5, help='The flank region of gene synteny.')
    subparser_partition.add_argument('--max_targets', '-k', required=False, type=int,
                                     default=2000, help='The maximum targets for each query in alignment. Improves accuracy for large-scale analyses, but increases runtime and memory usage.')
    subparser_partition.add_argument('--LD', required=False, type=float,
                                     default=0.6, help='Minimum gene length difference proportion between two genes.')
    subparser_partition.add_argument('--AS', required=False, type=float,
                                     default=0.6, help='Coverage for the shorter sequence.')
    subparser_partition.add_argument('--AL', required=False, type=float,
                                     default=0.6, help='Coverage for the longer sequence.')
    subparser_partition.add_argument('--evalue', required=False, type=float,
                                     default=1E-5, help='The evalue of aligner.')
    subparser_partition.add_argument('--aligner', required=False, type=str,
                                     default='diamond', choices=('diamond', 'blastp'), help='The aligner used to pairwise alignment.')
    subparser_partition.add_argument('--clust_method', required=False, type=str,
                                     default='mmseqs2', choices=('cdhit', 'mmseqs2'), help='The method used to cluster the genes.')
    subparser_partition.add_argument('--radius', required=False, type=int,
                                     default=3, help='The radius of search region.')
    subparser_partition.add_argument('--min_falen', '-m', required=False, type=check_min_falen,
                                     default=20, help='protein length of throw_away_sequences, at least 11')
    subparser_partition.add_argument('--gcode', required=False, type=check_gcode,
                                     default=11, help='The genetic code of your species.')
    subparser_partition.add_argument('--annot', required=False,
                                     action='store_true', default=False, help='Discard original annotation, and re-annote the genome privately using prodigal')

    subparser_partition.add_argument(
        '--verbose', '-v', required=False, action='store_true', default=False, help='Verbose output')
    subparser_partition.add_argument(
        '--debug', '-D', required=False, action='store_true', default=False, help='Debug mode. Note: very verbose')
    subparser_partition.add_argument(
        '--disable', required=False, action='store_true', default=False, help='Disable progress bar')

    subparser_partition.set_defaults(func=partition_portal)
