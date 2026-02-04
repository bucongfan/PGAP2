import os
import pickle
import shlex  # Added shlex
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
from pgap2.utils.generate_tree import merge_tree, generate_tree
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


def merge_pangenome(previous_pg: Pangenome, new_pg: Pangenome) -> Pangenome:
    """
    Merge a new pangenome with the previous pangenome.

    This function combines the strain information from both pangenomes,
    ensuring that the new genome's strain indices continue from where 
    the previous pangenome ended.

    Args:
        previous_pg: The existing pangenome object with previous strains
        new_pg: The new pangenome object containing newly added strain(s)

    Returns:
        Pangenome: A merged pangenome object containing all strains
    """
    logger.info(
        f'Merging pangenomes: {previous_pg.strain_num} previous strains + {new_pg.strain_num} new strain(s)')

    # Create a new merged pangenome with the same configuration as previous
    merged_pg = Pangenome(
        outdir=new_pg.outdir,  # Use new output directory
        threads=new_pg.threads,
        gcode=new_pg.gcode,
        disable=new_pg.disable_tqdm
    )

    # Copy all strains from previous pangenome
    for strain_index, strain in previous_pg.strain_dict.items():
        merged_pg.load_strain(strain)

    # Add all strains from new pangenome (they should already have correct indices)
    for strain_index, strain in new_pg.strain_dict.items():
        merged_pg.load_strain(strain)
        logger.debug(
            f'Added new strain: {strain.strain_name} (index: {strain_index})')

    # Merge total gene numbers
    merged_pg.total_gene_num = previous_pg.total_gene_num + new_pg.total_gene_num

    # Copy other important attributes from previous pangenome
    merged_pg.orth_id = previous_pg.orth_id
    merged_pg.para_id = previous_pg.para_id
    merged_pg.dup_id = previous_pg.dup_id
    merged_pg.accurate = previous_pg.accurate
    merged_pg.exhaust_orth = previous_pg.exhaust_orth
    merged_pg.retrieve = previous_pg.retrieve
    merged_pg.evalue = previous_pg.evalue
    merged_pg.aligner = previous_pg.aligner
    merged_pg.LD = previous_pg.LD
    merged_pg.AL = previous_pg.AL
    merged_pg.AS = previous_pg.AS

    logger.info(
        f'Merged pangenome: {merged_pg.strain_num} total strains, {merged_pg.total_gene_num} total genes')

    return merged_pg


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
        try:
            # Choose the largest one
            repre_with_max_len = max(
                tmp_len := {k: G.nodes[k]['length'] for k in max_pool}, key=tmp_len.get)
            logger.trace(
                f'Conflict in {each_gene} with {max_pool}, choose {repre_with_max_len}')
            return repre_with_max_len
        except:
            logger.trace(f'Would only occur in add mode')
            return None
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
    safe_mcl_result = shlex.quote(mcl_result)
    safe_alignment_result = shlex.quote(tree.alignment_result)
    run_command(
        f"{sfw.mcxdeblast} -m9 --score r --line-mode=abc {safe_alignment_result} 2> /dev/null | {sfw.mcl} - --abc -I 1.5 -te {pg.threads} -o {safe_mcl_result} &>/dev/null")

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


def merge_files(previous_file: str, current_file: str, output_file: str, skip_header: bool = False) -> str:
    """
    Merge the previous file with the current file.

    Args:
        previous_file: Path to the previous file
        current_file: Path to the current (new) file
        output_file: Path to the output merged file
        skip_header: If True, skip the first line (header) of the current file when appending

    Returns:
        str: Path to the merged file
    """
    logger.info(f'Merging files:')
    logger.info(f'  Previous: {previous_file}')
    logger.info(f'  Current:  {current_file}')

    with open(output_file, 'w') as out_fh:
        # Write previous file content first
        with open(previous_file, 'r') as fh:
            out_fh.write(fh.read())

        # Append current file content
        with open(current_file, 'r') as fh:
            if skip_header:
                # Skip header line
                next(fh, None)
            out_fh.write(fh.read())

    logger.info(f'  Merged:   {output_file}')

    return output_file


def build_mcl_member_dict(mcl_result_file: str) -> dict:
    """
    Build a dictionary mapping each gene to its MCL cluster ID.

    Args:
        mcl_result_file: Path to the MCL result file

    Returns:
        dict: {gene_id: mcl_cluster_id}
    """
    member_to_mcl = {}
    with open(mcl_result_file, 'r') as fh:
        for cluster_id, line in enumerate(fh):
            genes = line.strip().split('\t')
            for gene in genes:
                member_to_mcl[gene] = cluster_id
    return member_to_mcl


def formate_tree(tree: Tree, previous_tree: Tree):
    logger.info('---- Formatting tree leaves to match previous run...')
    prev_member_leaf = previous_tree.member_leaf
    if not prev_member_leaf:
        logger.warning(
            'Previous tree has no member_leaf mapping; skip formatting')
        assert False, 'Previous tree has no member_leaf mapping; cannot format current tree.'

    if not tree.leaf_member:
        tmp_leaf_member = defaultdict(set)
        for member, leaf in tree.member_leaf.items():
            tmp_leaf_member[leaf].add(member)
        tree.leaf_member = tmp_leaf_member

    leaf_map = {}
    target_sources = defaultdict(list)
    for leaf, members in tree.leaf_member.items():
        prev_counts = Counter(
            prev_member_leaf[member] for member in members if member in prev_member_leaf)
        if prev_counts:
            target, _ = prev_counts.most_common(1)[0]
            if len(prev_counts) > 1:
                logger.warning(
                    f'[w1] Leaf {leaf} maps to multiple previous leaves {dict(prev_counts)}; using {target}')
            leaf_map[leaf] = target
        else:
            # only occurs when new gene as a singleton
            logger.warning(
                f'[w2] Leaf {leaf} has no members in previous tree; keeping as is')
            leaf_map[leaf] = leaf
        target_sources[leaf_map[leaf]].append(leaf)

    for target, sources in target_sources.items():  # target is a node in previous tree
        if len(sources) > 1:
            logger.warning(
                f'Multiple current leaves {sources} map to previous leaf {target}; merging')

    new_leaf_member = defaultdict(set)
    for leaf, members in tree.leaf_member.items():
        target = leaf_map.get(leaf, leaf)
        new_leaf_member[target].update(members)

    new_member_leaf = {}
    for leaf, members in new_leaf_member.items():
        for member in members:
            new_member_leaf[member] = leaf

    new_leaf_member_strains = defaultdict(set)
    for leaf, members in new_leaf_member.items():
        new_leaf_member_strains[leaf] = {int(m.split(':')[0]) for m in members}

    def merge_graph(graph: nx.Graph, mapping: dict) -> nx.Graph:
        new_graph = nx.Graph()
        for node, data in graph.nodes(data=True):
            target = mapping.get(node, node)
            if target not in new_graph:
                new_graph.add_node(target, **data)
            else:
                for key, value in data.items():
                    if key in ('members', 'strains') and isinstance(value, (set, list, tuple)):
                        new_graph.nodes[target].setdefault(
                            key, set()).update(value)
                    elif key not in new_graph.nodes[target]:
                        new_graph.nodes[target][key] = value
        for u, v, edata in graph.edges(data=True):
            tu, tv = mapping.get(u, u), mapping.get(v, v)
            if tu == tv:
                continue
            if new_graph.has_edge(tu, tv):
                if 'weight' in edata:
                    prev = new_graph[tu][tv].get('weight', edata['weight'])
                    new_graph[tu][tv]['weight'] = max(prev, edata['weight'])
            else:
                new_graph.add_edge(tu, tv, **edata)
        return new_graph

    if tree.distance_graph:
        tree.distance_graph = merge_graph(tree.distance_graph, leaf_map)
    if hasattr(tree, 'raw_distance_graph') and tree.raw_distance_graph:
        tree.raw_distance_graph = merge_graph(
            tree.raw_distance_graph, leaf_map)

    ortho_tree = tree.orth_identity_tree
    if ortho_tree and ortho_tree.number_of_nodes() > 0:
        for leaf, target in list(leaf_map.items()):
            if leaf == target or leaf not in ortho_tree.nodes:
                continue
            if target in ortho_tree.nodes:
                for key in ('members', 'strains'):
                    if key in ortho_tree.nodes[leaf]:
                        ortho_tree.nodes[target].setdefault(key, set()).update(
                            ortho_tree.nodes[leaf][key])
                if 'members' in ortho_tree.nodes[target] and 'strains' in ortho_tree.nodes[target]:
                    ortho_tree.nodes[target]['has_para'] = len(
                        ortho_tree.nodes[target]['members']) != len(ortho_tree.nodes[target]['strains'])
                parents = list(ortho_tree.predecessors(leaf))
                for parent in parents:
                    if not ortho_tree.has_edge(parent, target):
                        ortho_tree.add_edge(parent, target)
                    if ortho_tree.has_edge(parent, leaf):
                        ortho_tree.remove_edge(parent, leaf)
                ortho_tree.remove_node(leaf)
            else:
                nx.relabel_nodes(ortho_tree, {leaf: target}, copy=False)
        tree.orth_identity_tree = ortho_tree

    new_leaf_root = {}
    next_root = max(previous_tree.leaf_root.values(), default=-1) + 1
    for leaf in new_leaf_member:
        if leaf in previous_tree.leaf_root:
            new_leaf_root[leaf] = previous_tree.leaf_root[leaf]
        else:
            new_leaf_root[leaf] = next_root
            next_root += 1

    tree.member_leaf = new_member_leaf
    tree.leaf_member = new_leaf_member
    tree.leaf_member_strains = new_leaf_member_strains
    tree.leaf_root = new_leaf_root

    root_leaf = defaultdict(set)
    for leaf, root in new_leaf_root.items():
        root_leaf[root].add(leaf)
    tree.root_leaf = root_leaf

    tree.ortho_para = {}
    if tree.orth_identity_tree and tree.orth_identity_tree.number_of_nodes() > 0:
        for node in tree.orth_identity_tree.nodes:
            if tree.orth_identity_tree.out_degree(node) == 0 and 'members' in tree.orth_identity_tree.nodes[node]:
                for member in tree.orth_identity_tree.nodes[node]['members']:
                    tree.ortho_para[member] = node


def merge_network(G: nx.Graph, pg: Pangenome, tree: Tree, new_strain_index: int, mcl_result_file: str) -> tuple:
    """
    Merge a new genome into the existing network G.

    This function:
    1. Creates individual nodes for each new gene
    2. Connects new nodes based on gene adjacency
    3. Merges new nodes with existing nodes based on MCL clustering and synteny

    Merge rules:
    - If a MCL cluster has only one new gene and one or more old nodes:
      - If one old node: merge directly
      - If multiple old nodes (split): use synteny to choose the best match
    - If a MCL cluster has multiple new genes: do not merge (keep separate)

    Args:
        G: Existing network graph
        pg: Pangenome object (contains all strains including the new one)
        tree: Tree object with MCL clustering info
        new_strain_index: Index of the newly added strain
        mcl_result_file: Path to the MCL result file

    Returns:
        tuple: (G, tree) - Updated graph and tree
    """
    logger.info(
        f'---- Merging new genome (strain index: {new_strain_index}) into the network...')

    # Get the new strain's gene information
    new_strain = pg.strain_dict[new_strain_index]
    gene_num = new_strain.gene_num

    # Step 1: Create nodes for new genes
    logger.info(f'---- Creating nodes for new genes...')
    new_nodes = []
    new_member2node = {}  # new gene -> its node name

    for contig_index, gene_count in enumerate(gene_num):
        for gene_index in range(gene_count):
            gene_id = f'{new_strain_index}:{contig_index}:{gene_index}'
            falen = pg.falen.get(gene_id, 0)
            new_nodes.append((gene_id, {
                'length': falen,
                'members': {gene_id},
                'strains': {new_strain_index},
                'has_para': False,
                'repre_nodes': [gene_id]
            }))
            new_member2node[gene_id] = gene_id

    G.add_nodes_from(new_nodes)
    logger.info(f'     Added {len(new_nodes)} new gene nodes')

    # Step 2: Add edges between adjacent new genes
    logger.info(f'---- Connecting adjacent new genes...')
    new_edges = []
    for contig_index, gene_count in enumerate(gene_num):
        if gene_count > 1:
            prev_gene = f'{new_strain_index}:{contig_index}:0'
            for gene_index in range(1, gene_count):
                curr_gene = f'{new_strain_index}:{contig_index}:{gene_index}'
                new_edges.append((prev_gene, curr_gene))
                prev_gene = curr_gene
    G.add_edges_from(new_edges)
    logger.info(f'     Added {len(new_edges)} edges between new genes')

    # Step 3: Build mapping from leaf -> old nodes / new genes using tree members
    logger.info(
        f'---- Building leaf-to-node mapping based on tree members...')

    leaf_to_old_nodes = defaultdict(set)
    for node in list(G.nodes()):
        if node in new_member2node:
            continue
        for member in G.nodes[node].get('members', set()):
            if member not in tree.member_leaf:
                continue
            if int(member.split(':')[0]) == new_strain_index:
                continue
            leaf_to_old_nodes[tree.member_leaf[member]].add(node)

    leaf_to_new_genes = defaultdict(list)
    for gene_id in new_member2node:
        leaf = tree.member_leaf.get(gene_id)
        if leaf is None:
            assert False, f'New gene {gene_id} not found in tree.member_leaf'
        leaf_to_new_genes[leaf].append(gene_id)

    # Step 4: Merge new nodes with old nodes based on leaf membership
    logger.info(
        f'---- Merging new nodes with existing nodes based on leaf membership...')

    merged_count = 0
    skipped_multi_new = 0
    skipped_no_old = 0
    relabel_dict = {}
    split_result_map = defaultdict(list)

    for leaf, new_genes in leaf_to_new_genes.items():
        if len(new_genes) == 0:
            assert False, f'Leaf {leaf} has no new genes'
        if len(new_genes) > 1:
            logger.debug(
                f'Leaf {leaf}: {len(new_genes)} new genes, skipping merge')
            skipped_multi_new += 1
            for new_gene in new_genes:
                new_node = f'{leaf}_{new_gene}'
                relabel_dict[new_gene] = new_node
                G.nodes[new_gene]['repre_nodes'] = [new_node]
                split_result_map[leaf].append(new_node)
            continue

        old_nodes = leaf_to_old_nodes.get(leaf, set())
        if not old_nodes:
            logger.debug(
                f'Leaf {leaf}: no old nodes, new gene stays as singleton')
            skipped_no_old += 1
            assert False, f'Leaf {leaf} has no old nodes'

        new_gene = new_genes[0]
        target_node = None
        if len(old_nodes) > 1:
            logger.debug(
                f'Leaf {leaf}: multiple old nodes found, using synteny to select best match')
            rep_to_node = {}
            para_repre = []
            para_context = {}

            for old_node in old_nodes:
                members = list(G.nodes[old_node].get('members', set()))
                repre = members[0] if members else old_node
                rep_to_node[repre] = old_node
                para_repre.append(repre)
                para_context[repre] = tree.get_context(repre, flank=10)

            para_context[new_gene] = tree.get_context(new_gene, flank=10)
            best_repre = select_repre_node(
                G, tree, new_gene, para_repre, para_context)
            if best_repre is not None:
                target_node = rep_to_node.get(best_repre)

        if target_node is None:
            target_node = max(
                old_nodes, key=lambda n: G.nodes[n].get('length', 0))
        G = merge_node(
            G, pg, tree, [target_node, new_gene], target=target_node)
        merged_count += 1
        logger.debug(f'Merged {new_gene} into {target_node}')

    if relabel_dict:
        nx.relabel_nodes(G, relabel_dict, copy=False)
        merged_split_map = defaultdict(list)
        if hasattr(tree, '_split_result_map_reverse'):
            for child, parent in tree._split_result_map_reverse.items():
                merged_split_map[parent].append(child)
        for parent, children in split_result_map.items():
            merged_split_map[parent].extend(children)
        tree.load_split_result_map(merged_split_map)
    logger.info(f'     Merged: {merged_count} new genes')
    logger.info(
        f'     Skipped (multiple new genes in MCL cluster): {skipped_multi_new}')
    logger.info(f'     Skipped (no old nodes): {skipped_no_old}')

    # Step 5: Update tree structures for remaining new singleton nodes
    logger.info(f'---- Updating tree structures for new nodes...')

    # Rebuild member_leaf from final G (members are always base genes)
    member_leaf = {}
    leaf_member = {}
    for node in G.nodes():
        assert 'members' in G.nodes[node], f'Node {node} missing members attribute'
        for member in G.nodes[node].get('members', set()):
            member_leaf[member] = node

    # Rebuild leaf_root from final G and orth_identity_tree
    leaf_root = {}
    for node in G.nodes():
        if node in tree.leaf_root:
            leaf_root[node] = tree.leaf_root[node]
        elif '_' in node:
            father = node.split('_')[0]
            leaf_root[node] = tree.leaf_root.get(father, father)
        else:
            members = G.nodes[node].get('members', set())
            for member in members:
                if member in tree.leaf_root:
                    leaf_root[node] = tree.leaf_root[member]
                    break
            else:
                print(G.nodes[node])
                assert False, f'Node {node} not found in tree.leaf_root'

    tree.member_leaf = member_leaf
    tree.leaf_root = leaf_root

    # Rebuild root_leaf/leaf_member/leaf_member_strains to match generate_network
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


def main(indir: str, outdir: str, previous_dir: str, aligner: str, clust_method: str, falen: int, threads: int, id_attr_key: str, type_filter: str, annot: bool, gcode: int, retrieve: bool, disable: bool):

    # Load previous parameters
    logger.info(f'Loading parameters from {previous_dir}/basic.pkl')
    with open(f'{previous_dir}/basic.pkl', 'rb') as fh:
        previous: PklCheck = pickle.load(fh)
        basic = previous.data_dump('basic')
        params = basic.params

    # Log current arguments
    logger.info('Current arguments:')
    current_args = {
        'aligner': aligner,
        'clust_method': clust_method,
        'falen': falen,
        'threads': threads,
        'id_attr_key': id_attr_key,
        'type_filter': type_filter,
        'annot': annot,
        'gcode': gcode,
        'retrieve': retrieve,
        'disable': disable
    }
    for k, v in current_args.items():
        logger.info(f'{k}: {v}')

    # Check for consistency
    check_keys = ['aligner', 'clust_method', 'falen',
                  'id_attr_key', 'type_filter', 'annot', 'gcode', 'retrieve']
    for key in check_keys:
        if key in params and params[key] != current_args[key]:
            logger.warning(
                f"Parameter mismatch: {key} (Previous: {params[key]}, Current: {current_args[key]})")

    evalue = params['evalue']
    hconf_thre = params['hconf_thre']
    fast_mode = params['fast_mode']
    orth_id = params['orth_id']
    para_id = params['para_id']
    dup_id = params['dup_id']
    max_targets = params['max_targets']
    coverage = params.get('coverage', 0.98)
    LD = params['LD']
    AS = params['AS']
    AL = params['AL']
    context_similarity = params.get('context_similarity', 0)
    accurate = params['accurate']
    exhaust_orth = params['exhaust_orth']
    flank = params.get('flank', 5)
    radius = params.get('radius', 3)
    sensitivity = params.get('sensitivity', 'strict')
    ins = params.get('ins', False)

    logger.info('Inherited parameters:')
    inherited_params = {
        'evalue': evalue, 'hconf_thre': hconf_thre, 'fast_mode': fast_mode,
        'orth_id': orth_id, 'para_id': para_id, 'dup_id': dup_id,
        'max_targets': max_targets, 'coverage': coverage,
        'LD': LD, 'AS': AS, 'AL': AL,
        'context_similarity': context_similarity, 'accurate': accurate,
        'exhaust_orth': exhaust_orth, 'flank': flank,
        'radius': radius, 'sensitivity': sensitivity, 'ins': ins
    }
    for k, v in inherited_params.items():
        logger.info(f'{k}: {v}')

    decode_status = False
    file_dict = get_file_dict(indir)

    # Validate that only one genome is being added
    if len(file_dict) != 1:
        logger.error(
            f'The "add" command is designed to add only ONE genome at a time.')
        logger.error(f'Found {len(file_dict)} genome(s) in {indir}:')
        for strain_name in file_dict.keys():
            logger.error(f'  - {strain_name}')
        logger.error(
            f'Please ensure the input directory contains exactly one genome file.')
        logger.error(
            f'If you want to add multiple genomes, add them one by one, or use the "main" command to start a new analysis.')
        raise ValueError(
            f'Invalid input: expected 1 genome, found {len(file_dict)} genomes in {indir}')

    logger.info(f'Validated input: exactly 1 genome to be added from {indir}')

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
        logger.info(f'Loading...')
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

                    # Get the new strain index (it's the last one added)
                    # Load previous_dir pangenome to get the original strain count
                    with open(f'{previous_dir}/preprocess.pkl', 'rb') as prev_fh:
                        prev_pkl: PklCheck = pickle.load(prev_fh)
                        prev_pg = prev_pkl.data_dump('pangenome')
                        start_strain_index = prev_pg.strain_num
                        logger.info(f'New strain index: {start_strain_index}')
            else:
                logger.warning(
                    f'Previous file parameters is not match, start partition from the begining')

    if decode_status is False:
        '''
        Load strain from input directory
        If the previous preprocess.pkl file is not found or the parameters are not match,
        it will load the strain from the input directory and create a new pangenome object.
        '''
        # load single genome
        logger.info(f'Load strain from {indir}')

        # Load previous pangenome to get the starting strain index
        with open(f'{previous_dir}/preprocess.pkl', 'rb') as fh:
            previous_pkl: PklCheck = pickle.load(fh)
            previous_pg = previous_pkl.data_dump('pangenome')
            previous_tree = previous_pkl.data_dump('tree')
            start_strain_index = previous_pg.strain_num
            logger.info(f'Starting strain index: {start_strain_index}')

        # Parse new genome with correct starting index
        new_pg = file_parser(
            indir=indir, outdir=outdir, annot=annot, threads=1, disable=disable,
            retrieve=retrieve, falen=falen, gcode=gcode, id_attr_key=id_attr_key,
            type_filter=type_filter, prefix='partition', start_index=start_strain_index, run_type='add')
        # Merge previous pangenome with new pangenome
        pg = merge_pangenome(previous_pg, new_pg)

        file_prot = f'{outdir}/add.involved_prot.fa'
        file_annot = f'{outdir}/add.involved_annot.tsv'

        # Merge protein and annotation files
        merged_prot = merge_files(
            previous_file=f'{previous_dir}/total.involved_prot.fa',
            current_file=file_prot,
            output_file=f'{outdir}/total.involved_prot.fa',
            skip_header=False
        )
        merged_annot = merge_files(
            previous_file=f'{previous_dir}/total.involved_annot.tsv',
            current_file=file_annot,
            output_file=f'{outdir}/total.involved_annot.tsv',
            skip_header=True  # Skip header in new annot file
        )

        pg.load_annot_file(merged_annot)
        pg.load_prot_file(merged_prot)
        logger.info(
            f'Create distane tree with {pg.strain_num} strains')
        logger.info(
            f'Clustering with orth_id: {orth_id}, para_id: {para_id}, dup_id: {dup_id}')
        tree = generate_tree(
            input_file=merged_prot, orth_list=[dup_id, orth_id], outdir=pg.outdir, evalue=evalue, aligner=aligner, falen=falen, disable=disable, threads=threads, max_targets=max_targets, coverage=coverage, ID=para_id, LD=LD, AS=AS, AL=AL, clust_method=clust_method)

        logger.info(f'Pangenome and tree loaded successfully.')
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

    with open(f'{previous_dir}/graph.pkl', 'rb') as fh:
        previous_G: PklCheck = pickle.load(fh)
        previous_G = previous_G.data_dump('graph')

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

    with open(f'{previous_dir}/preprocess.pkl', 'rb') as fh:
        previous_pkl: PklCheck = pickle.load(fh)
        previous_tree = previous_pkl.data_dump('tree')
    # -----------------------------------partition step-----------------------------------#
    # logger.info('Get the gene primal clust result by mcl')
    mcl(pg, tree)  # add distance graph to tree
    formate_tree(tree, previous_tree)
    # reformat_tree(tree, previous_tree)
    mcl_result_file = f'{pg.outdir}/mcl.result'
    logger.info('Load the gene length information')
    pg.reload_nucl_file(tree)
    logger.info('Create synteny network')
    G, tree = merge_network(
        G=previous_G,
        pg=pg,
        tree=tree,
        new_strain_index=start_strain_index,
        mcl_result_file=mcl_result_file
    )
    # G, tree = generate_network(pg=pg, tree=tree)

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
                                      context_sim=context_similarity,
                                      flank=flank,
                                      disable=disable,)

    logger.info(f'Double check through gene synteny')
    G = merge_by_synteny(G, pg, tree,
                         context_sim=context_similarity,
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
                                          context_sim=context_similarity,
                                          flank=flank,
                                          disable=disable, step=10)
        G = merge_by_synteny(G=G, pg=pg, tree=tree,
                             context_sim=context_similarity,
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
         previous_dir=os.path.abspath(args.previous),
         aligner=args.aligner, clust_method=args.clust_method,
         falen=args.min_falen, threads=args.threads,
         id_attr_key=args.id_attr_key, type_filter=args.type_filter,
         annot=args.annot, gcode=args.gcode, retrieve=args.retrieve,
         disable=args.disable)


def add_portal(args):
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


def add_cmd(subparser: _SubParsersAction):

    subparser_add: ArgumentParser = subparser.add_parser(
        'add', help='Add a single genome to an existing PGAP2 project', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparser_add.add_argument(
        '--indir', '-i', required=True, help='Input file contained, same prefix seems as the same strain.',)
    subparser_add.add_argument(
        '--outdir', '-o', required=False, help='Output directory', default='.',)
    subparser_add.add_argument(
        '--previous', '-p', required=False, help='The previous PGAP2 result directory, used to resume the partition step quickly.', default=None,)
    subparser_add.add_argument('--aligner', required=False, type=str,
                               default='diamond', choices=('diamond', 'blastp'), help='The aligner used to pairwise alignment.')
    subparser_add.add_argument('--clust_method', required=False, type=str,
                               default='mmseqs2', choices=('cdhit', 'mmseqs2'), help='The method used to cluster the genes.')

    subparser_add.add_argument("--type-filter", required=False, type=str,
                               default='CDS', help="Only for gff file as input, feature type (3rd column) to include, Only lines matching these types will be processed.")
    subparser_add.add_argument("--id-attr-key", required=False, type=str,
                               default='ID', help="Only for gff file as input, Attribute key to extract from the 9th column as the record ID (e.g., 'ID', 'gene', 'locus_tag').")
    subparser_add.add_argument('--retrieve', '-r', required=False,
                               action='store_true', default=False, help='Retrieve gene that may lost with annotations')
    subparser_add.add_argument(
        '--threads', '-t', required=False, default=1, help='threads used in parallel', type=int)
    subparser_add.add_argument('--min_falen', '-m', required=False, type=check_min_falen,
                               default=20, help='protein length of throw_away_sequences, at least 11')
    subparser_add.add_argument('--gcode', required=False, type=check_gcode,
                               default=11, help='The genetic code of your species.')
    subparser_add.add_argument('--annot', required=False,
                               action='store_true', default=False, help='Discard original annotation, and re-annote the genome privately using prodigal')

    subparser_add.add_argument(
        '--verbose', '-v', required=False, action='store_true', default=False, help='Verbose output')
    subparser_add.add_argument(
        '--debug', '-D', required=False, action='store_true', default=False, help='Debug mode. Note: very verbose')
    subparser_add.add_argument(
        '--disable', required=False, action='store_true', default=False, help='Disable progress bar')

    subparser_add.set_defaults(func=add_portal)
