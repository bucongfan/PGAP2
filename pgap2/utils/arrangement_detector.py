import itertools
import networkx as nx

from tqdm import tqdm
from collections import defaultdict

from pgap2.lib.pangenome import Pangenome
from pgap2.lib.tree import Tree
from pgap2.utils.graph_utils import merge_node, get_similarity, merge_judge
from pgap2.utils.supply import tqdm_

"""
This module provides functions to handle rearrangements in a pangenome graph based on context similarity.

After core gene partitioning, this method identifies and resolves recent duplication events in genomes 
by detecting highly similar nodes that may represent redundant gene copies. 
Such nodes are candidates for merging if they exhibit extremely high sequence identity and 
are located in similar genomic contexts. This helps refine the pangenome graph by reducing redundancy 
and improving structural coherence.

input:
- G: Pangenome graph containing gene clusters.
- pg: Pangenome object containing strain and gene information.
- tree: Identity tree of the genes.
- context_sim: Context similarity threshold for merging nodes.
- flank: Flank size for context comparison.
- sensitivity: Sensitivity level for merging decisions.
- ins: Boolean indicating whether to consider insertion events.
output:
- G: Updated Pangenome graph after merging nodes based on context similarity.
"""


def get_flank_clust(pg: Pangenome, tree: Tree, repre_node, flank):
    repre_flank_gene = pg.get_flank_gene(gene=repre_node, flank=flank)
    repre_clust = {}
    for each_flank_gene in repre_flank_gene:
        repre = tree.ancestor(each_flank_gene)
        repre_clust.update({each_flank_gene: repre})
    repre_clust.update({repre_node: tree.ancestor(repre_node)})
    return repre_clust


def get_rearrange_list(G, tree):
    rearrange_list = []
    for root, leaf_nodes in tree.iter_root():
        tmp_node_list = []
        for each_leaf in leaf_nodes:
            if G.has_node(each_leaf):
                tmp_node_list.append(each_leaf)
        if len(tmp_node_list) > 1:
            for a, b in itertools.combinations(tmp_node_list, 2):
                rearrange_list.append((a, b))
    return rearrange_list


def _guide_says_merge(G, guide_dict, node_a, node_b):
    """
    Use the guide_dict to decide whether two nodes should be merged.

    Returns:
        True  - all members are old genes and shared a previous cluster → merge.
        False - all members are old genes but different clusters → don't merge.
        None  - has new genes, guide cannot decide → normal judgment.
    """
    # O(1) short-circuit: if either node has new genes, guide cannot decide
    if G.nodes[node_a].get('has_new', False) or G.nodes[node_b].get('has_new', False):
        return None

    members_a = G.nodes[node_a]['members']
    members_b = G.nodes[node_b]['members']

    clusters_a = set()
    for m in members_a:
        c = guide_dict.get(m)
        if c is None:
            return None
        clusters_a.add(c)

    for m in members_b:
        c = guide_dict.get(m)
        if c is None:
            return None
        if c in clusters_a:
            return True  # found intersection, merge

    return False


def merge_by_synteny(G: nx.Graph, pg: Pangenome, tree: Tree, context_sim: float, flank: int, sensitivity: str, ins: bool = False, step: int = 4, guide_dict: dict = None, new_strain_indices: set = None):

    def find_final_node(node, mapping):
        while mapping[node] != node:
            node = mapping[node]
        return node

    use_guide = guide_dict is not None

    for _, nodes in tqdm(tree.root_leaf.items(), unit=' node', desc=tqdm_.step(step=step), disable=pg.disable_tqdm):
        if len(nodes) == 1:
            continue
        exists_nodes = [_ for _ in nodes if G.has_node(_)]

        # --- Guide fast path: skip pairwise search for all-old-gene groups ---
        if use_guide and all(not G.nodes[n].get('has_new', False) for n in exists_nodes):
            cluster_groups = defaultdict(list)
            for n in exists_nodes:
                for m in G.nodes[n]['members']:
                    c = guide_dict.get(m)
                    if c is not None:
                        cluster_groups[c].append(n)
                        break
            for cluster_id, group_nodes in cluster_groups.items():
                unique_nodes = list(dict.fromkeys(group_nodes))
                alive = [n for n in unique_nodes if G.has_node(n)]
                if len(alive) < 2:
                    continue
                target = max(alive, key=lambda x: G.nodes[x]['length'])
                G = merge_node(G, pg, tree, alive, target=target)
                for v in alive:
                    if v != target and pg.retrieve:
                        tree.update_removed_nodes(v)
            continue

        need_merge_nodes = {}
        for clust_pair in itertools.combinations(exists_nodes, 2):
            need_merge = False
            clust_a, clust_b = clust_pair
            if not G.has_node(clust_a) or not G.has_node(clust_b):
                continue
            inter_lca = 0
            for a, b in itertools.product(G.nodes[clust_a]['repre_nodes'], G.nodes[clust_b]['repre_nodes']):
                if tree.distance_graph.has_edge(a, b):
                    this_lca = tree.distance_graph[a][b]['weight']
                    inter_lca = max(inter_lca, this_lca)

            if inter_lca == 0:  # not connected
                continue
            need_merge_nodes[clust_pair] = inter_lca

        if not need_merge_nodes:
            continue

        sored_need_merge_nodes = sorted(
            need_merge_nodes.items(), key=lambda x: x[1], reverse=True)
        split_clust_map = {_: _ for _ in exists_nodes}
        chacned_result = {}
        for clust_pair, lca in sored_need_merge_nodes:
            clust_a, clust_b = clust_pair
            clust_a = find_final_node(clust_a, split_clust_map)
            clust_b = find_final_node(clust_b, split_clust_map)
            if clust_a == clust_b:
                continue

            u_i = len(G.nodes[clust_a]['repre_nodes'])
            v_i = len(G.nodes[clust_b]['repre_nodes'])
            flag = False
            need_merge = False

            # --- Guide-based shortcut for add mode ---
            if use_guide:
                guide_result = _guide_says_merge(
                    G, guide_dict, clust_a, clust_b)
                if guide_result is not None:
                    need_merge = guide_result
                    flag = True

            if (clust_a, clust_b) in chacned_result and not flag:
                pre_u_i, pre_v_i, pre_need_merge = chacned_result[(
                    clust_a, clust_b)]
                if u_i == pre_u_i and v_i == pre_v_i:
                    need_merge = pre_need_merge
                    flag = True

            if not flag:
                if (not pg.exhaust_orth) and (lca >= tree.dup_id):
                    clust_a_real_repre = tree.get_unsplit_repre(
                        G.nodes[clust_a]['repre_nodes'])
                    clust_b_real_repre = tree.get_unsplit_repre(
                        G.nodes[clust_b]['repre_nodes'])
                    if set(clust_a_real_repre).intersection(set(clust_b_real_repre)):
                        if not ins:
                            need_merge = True
                        else:
                            sim = get_similarity(tree.get_context(clust_a, flank=flank), tree.get_context(
                                clust_b, flank=flank))
                            if sim <= context_sim:
                                need_merge = False
                            else:
                                need_merge = True
                    else:
                        need_merge = False

                else:
                    sim = get_similarity(tree.get_context(clust_a, flank=flank), tree.get_context(
                        clust_b, flank=flank))
                    if sim <= context_sim:
                        need_merge = False
                    else:
                        need_merge = merge_judge(tree, G, pg, clust_a, clust_b, lca,
                                                 context_sim, flank, sensitivity)
                chacned_result[(clust_a, clust_b)] = (u_i, v_i, need_merge)
                chacned_result[(clust_b, clust_a)] = (v_i, u_i, need_merge)

            if need_merge:
                clust_a, clust_b = (clust_a, clust_b) if G.nodes[clust_a][
                    'length'] > G.nodes[clust_b]['length'] else (clust_b, clust_a)
                G = merge_node(G, pg, tree, [clust_a, clust_b], target=clust_a)
                split_clust_map[clust_b] = clust_a
                if pg.retrieve:
                    tree.update_removed_nodes(clust_b)
    return G
