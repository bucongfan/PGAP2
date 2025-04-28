import itertools
import numpy as np
import networkx as nx

from tqdm import tqdm
from collections import defaultdict
from pgap2.utils.supply import tqdm_


class Tree(nx.DiGraph):
    def __init__(self, *args, **kwargs):
        super(Tree, self).__init__(*args, **kwargs)
        self.leaf_root = {}
        self.root_leaf = {}
        self.member_leaf = {}
        self.leaf_member = {}
        self.leaf_member_strains = {}
        self.ortho_para = {}
        self.pan_judge_dict = {}
        self.para_id = 0.7
        self.orth_id = 0.98
        self.dup_id = 0.99
        self.distance_graph = nx.Graph()
        self.orth_identity_tree = nx.DiGraph()

    def load_expect_identity(self, expect_identity):
        self.expect_identity = expect_identity

    def load_split_result_map(self, split_result_map):
        self._split_result_map_reverse = {}
        for key, value in split_result_map.items():
            for v in value:
                self._split_result_map_reverse[v] = key

    def load_alignment_result(self, alignment_result):
        self.alignment_result = alignment_result

    def load_mcl_result(self, mcl_result):
        self.mcl_result = mcl_result

    def load_distance_graph(self, distance_graph: nx.Graph, raw: bool = False):
        if raw is True:
            self.raw_distance_graph = distance_graph
            return
        else:
            self.distance_graph = distance_graph
            for root_i, compnent in enumerate(nx.connected_components(distance_graph)):
                for leaf in compnent:
                    self.leaf_root[leaf] = root_i

    def get_unsplit_repre(self, nodes):
        real_nodes = []
        for node in nodes:
            if node in self._split_result_map_reverse:
                real_nodes.append(self._split_result_map_reverse[node])
        return real_nodes

    def are_in_same_clique(self, node_list):
        for i, node1 in enumerate(node_list):
            for node2 in node_list[i+1:]:
                if not self.distance_graph.has_edge(node1, node2):
                    return False
        return True

    def load_ortho_identity_tree(self, ortho_tree):
        self.orth_identity_tree = ortho_tree

        for node in ortho_tree.nodes:
            if ortho_tree.in_degree(node) == 0:
                self.member_leaf.update(
                    {member: node for member in ortho_tree.nodes[node]['members']})

                descendants = nx.descendants(ortho_tree, node)
                if descendants:
                    for leaf in (n for n in descendants if ortho_tree.out_degree(n) == 0):
                        self.ortho_para.update(
                            {member: leaf for member in ortho_tree.nodes[leaf]['members']})
                else:
                    self.ortho_para.update(
                        {member: node for member in ortho_tree.nodes[node]['members']})

    def has_para(self, clusts):
        clust_a = clusts[0]
        strain_set = set([_.split(':')[0] for _ in clust_a])
        for clust_b in clusts[1:]:
            clust_b_strain_set = set([_.split(':')[0] for _ in clust_b])
            if strain_set & clust_b_strain_set:
                return True
            else:
                strain_set |= clust_b_strain_set
        return False

    def load_para_id(self, para_id):
        self.para_id = para_id

    def load_orth_id(self, orth_id=0.98):
        self.orth_id = orth_id

    def load_dup_id(self, dup_id=0.99):
        self.dup_id = dup_id

    def get_context(self, gene, flank=10):
        '''
        Get the context of a gene, returning 'flank' unique ancestors for both upstream and downstream,
        unless the gene sequence ends.
        '''
        if '_' in gene:
            strain_index, contig_index, gene_index = gene.split('_')[
                1].split(':')
        else:
            strain_index, contig_index, gene_index = gene.split(':')

        # 上下游独特祖先集合
        unique_ancestors_up = set()
        unique_ancestors_down = set()

        context = []

        i = 1  # 索引从1开始，即离目标基因最近的相邻基因
        gene_index = int(gene_index)  # 将基因索引转换为整数
        has_upper = True
        has_lower = True
        # 继续搜索，直到上下游各自收集到足够数量的独特祖先
        while has_upper or has_lower:
            # 计算上下游基因索引
            upper_gene_index = gene_index - i
            lower_gene_index = gene_index + i

            # 生成基因标识符
            upper_g = f'{strain_index}:{contig_index}:{upper_gene_index}'
            lower_g = f'{strain_index}:{contig_index}:{lower_gene_index}'

            # 检查上游基因祖先
            if has_upper and upper_g in self.member_leaf:
                ancestor = self.ancestor(self.member_leaf[upper_g])
                if ancestor not in context:
                    unique_ancestors_up.add(ancestor)
                    context.append(ancestor)
                    if len(unique_ancestors_up) == flank:
                        has_upper = False
            else:
                has_upper = False

            # 检查下游基因祖先
            if has_lower and lower_g in self.member_leaf:
                ancestor = self.ancestor(self.member_leaf[lower_g])
                if ancestor not in context:
                    unique_ancestors_down.add(ancestor)
                    context.append(ancestor)
                    if len(unique_ancestors_down) == flank:
                        has_lower = False
            else:
                has_lower = False

            i += 1  # 增加搜索范围

        return context

    def ancestor(self, node):
        return self.leaf_root[node]

    def build_index(self, disable=False):
        new_distance_graph = self._build_index(
            distance_graph=self.distance_graph, disable=disable)
        self.distance_graph = new_distance_graph
        self.distance_matrix: np.ndarray = nx.adjacency_matrix(
            self.distance_graph)
        self.distance_matrix.data[:] = int(1)
        # self.distance_matrix = self.distance_matrix.toarray()
        self.matrix_node_map = {node: idx for idx,
                                node in enumerate(self.distance_graph.nodes())}
        new_distance_graph = self._build_index(
            distance_graph=self.raw_distance_graph, disable=disable)
        self.raw_distance_graph = new_distance_graph

    def _build_index(self, distance_graph, disable=False):
        # Step 1: 拆分节点并更新图结构
        split_dict = defaultdict(list)
        for node in self.leaf_member.keys():
            if '_' in node:
                fa_node = node.split('_')[0]
                split_dict[fa_node].append(node)

        # 复制distance_graph的节点和边，构建H图
        H = nx.Graph()
        H.add_nodes_from(distance_graph.nodes())
        H.add_edges_from(distance_graph.edges(data=True))

        # 处理分裂节点的邻接关系
        for old_node, new_nodes in tqdm(split_dict.items(), unit=f" node",
                                        disable=disable, desc=tqdm_.step(3)):
            neighbors = list(H.neighbors(old_node))
            # 添加新节点并复制旧节点的邻接关系
            for new_node in new_nodes:
                H.add_node(new_node)
                for neighbor in neighbors:
                    if neighbor not in new_nodes:
                        H.add_edge(new_node, neighbor,
                                   weight=H[old_node][neighbor]['weight'])

            # 在新节点之间添加边
            for a, b in itertools.combinations(new_nodes, 2):
                if self.ortho_para[a.split('_')[1]] == self.ortho_para[b.split('_')[1]]:
                    H.add_edge(a, b, weight=self.dup_id)
                else:
                    H.add_edge(a, b, weight=self.orth_id)

            # 移除旧节点
            H.remove_node(old_node)
        return H
