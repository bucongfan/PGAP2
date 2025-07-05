import itertools
import networkx as nx

from tqdm import tqdm
from collections import defaultdict
from pgap2.utils.supply import tqdm_


class Tree(nx.DiGraph):
    def __init__(self, *args, **kwargs):
        super(Tree, self).__init__(*args, **kwargs)
        self.leaf_root = {}
        self.root_leaf = defaultdict(set)
        self.pseudo_root_leaf = defaultdict(set)
        self.member_leaf = {}
        self.leaf_member = {}
        self.leaf_member_strains = {}
        self.ortho_para = {}
        self.pan_judge_dict = {}
        self.para_id = 0.7
        self.orth_id = 0.98
        self.dup_id = 0.99
        self.removed_nodes = set()
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
        """
        Load the orthologous identity tree.
        """
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

        # unique ancestor sets for upstream and downstream
        unique_ancestors_up = set()
        unique_ancestors_down = set()

        context = []

        i = 1  # begin with the first gene index
        gene_index = int(gene_index)  # convert gene index to integer
        has_upper = True
        has_lower = True
        # continue searching until enough unique ancestors are collected for both upstream and downstream
        while has_upper or has_lower:
            # calculate upstream and downstream gene indices
            upper_gene_index = gene_index - i
            lower_gene_index = gene_index + i

            # generate gene identifiers
            upper_g = f'{strain_index}:{contig_index}:{upper_gene_index}'
            lower_g = f'{strain_index}:{contig_index}:{lower_gene_index}'

            # check upper gene ancestor
            if has_upper and upper_g in self.member_leaf:
                ancestor = self.ancestor(self.member_leaf[upper_g])
                if ancestor not in context:
                    unique_ancestors_up.add(ancestor)
                    context.append(ancestor)
                    if len(unique_ancestors_up) == flank:
                        has_upper = False
            else:
                has_upper = False

            # check lower gene ancestor
            if has_lower and lower_g in self.member_leaf:
                ancestor = self.ancestor(self.member_leaf[lower_g])
                if ancestor not in context:
                    unique_ancestors_down.add(ancestor)
                    context.append(ancestor)
                    if len(unique_ancestors_down) == flank:
                        has_lower = False
            else:
                has_lower = False

            i += 1  # increase the index to check the next gene
        return context

    def ancestor(self, node):
        return self.leaf_root[node]

    def get_removed_nodes(self):
        return self.removed_nodes

    def set_removed_nodes(self, removed_nodes):
        """ Set the removed nodes in the tree.
        This method is used to update the set of removed nodes in the tree.
        """
        self.removed_nodes = removed_nodes

    def update_removed_nodes(self, node):
        """ Update the removed nodes in the tree by adding a new node.
        This method is used to add a new node to the set of removed nodes.
        """
        if node not in self.removed_nodes:
            self.removed_nodes.add(node)
        else:
            raise ValueError(
                f"Node {node} is already in the removed nodes set.")

    def update_distance_graph(self, disable=False):
        """ Update the distance graph by rebuilding the index.
        This method is called after any changes to the distance graph or its related attributes.
        """
        self.distance_graph = self._build_index(
            distance_graph=self.distance_graph, disable=disable)
        self.raw_distance_graph = self._build_index(
            distance_graph=self.raw_distance_graph, disable=disable)
        self.update_distance_matrix()

    def update_distance_matrix(self):
        """
        Update the distance matrix based on the current distance graph.
        This method is called after the distance graph has been modified.
        """
        self.distance_matrix = nx.adjacency_matrix(self.distance_graph)
        self.distance_matrix.data[:] = int(1)
        self.matrix_node_map = {node: idx for idx,
                                node in enumerate(self.distance_graph.nodes())}
        for root, nodes in self.pseudo_root_leaf.items():
            self.root_leaf[root].update(nodes)

    def _build_index(self, distance_graph, disable=False):
        """ Build the index for the distance graph by splitting nodes and updating the graph structure.
        This method is called to create a new graph structure based on the existing distance graph.
        It handles the splitting of nodes that have multiple members and updates the graph accordingly.
        """
        # Step 1: Identify nodes that need to be split
        split_dict = defaultdict(list)
        for node in self.leaf_member.keys():
            if '_' in node:
                fa_node = node.split('_')[0]
                split_dict[fa_node].append(node)

        # Step 2: Create a new graph structure
        H = nx.Graph()
        H.add_nodes_from(distance_graph.nodes())
        H.add_edges_from(distance_graph.edges(data=True))

        # Step 3: Split nodes and update the graph
        for old_node, new_nodes in tqdm(split_dict.items(), unit=f" node",
                                        disable=disable, desc=tqdm_.step(3)):
            neighbors = list(H.neighbors(old_node))
            # add new nodes to the graph
            for new_node in new_nodes:
                H.add_node(new_node)
                for neighbor in neighbors:
                    if neighbor not in new_nodes:
                        H.add_edge(new_node, neighbor,
                                   weight=H[old_node][neighbor]['weight'])

            # Add edges between new nodes
            for a, b in itertools.combinations(new_nodes, 2):
                if self.ortho_para[a.split('_')[1]] == self.ortho_para[b.split('_')[1]]:
                    H.add_edge(a, b, weight=self.dup_id)
                else:
                    H.add_edge(a, b, weight=self.orth_id)

            # Remove the old node from the graph
            H.remove_node(old_node)
        return H
