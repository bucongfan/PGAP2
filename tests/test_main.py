
def test_file_parser():
    from pgap2.utils.data_loader import file_parser
    from pgap2.lib.pangenome import Pangenome
    import os
    import shutil
    input_dir = 'pgap2/toy_input'
    output_dir = 'tests/test_output'
    os.makedirs(output_dir, exist_ok=True)
    pg = file_parser(
        indir=input_dir, outdir=output_dir, annot=False, threads=1, disable=False, retrieve=False, falen=10, gcode=11, id_attr_key='CDS', type_filter='locus_tag', prefix='test')
    shutil.rmtree(output_dir)  # Clean up after test
    assert isinstance(
        pg, Pangenome), "file_parser did not return a Pangenome object correctly"
    pass


def test_get_file_dict():
    from pgap2.utils.data_loader import get_file_dict

    input_dir = 'pgap2/toy_input'
    file_dict = get_file_dict(input_dir)
    assert isinstance(
        file_dict, dict), "get_file_dict did not return a dictionary"
    pass


def test_decode_success():
    from pgap2.lib.pklcheck import PklCheck
    import pickle

    pkl_file = 'tests/output/preprocess.pkl'
    with open(pkl_file, 'rb') as fh:
        previous: PklCheck = pickle.load(fh)
        decode_status = previous.decode(
            orth_id=0.98, para_id=0.7, dup_id=0.99,
            accurate=False, coverage=0.98,
            id_attr_key='ID', type_filter='CDS',
            LD=0.6, AS=0.6, AL=0.6, evalue=1e-5, aligner='diamond',
            clust_method='cdhit', falen=20, annot=False, retrieve=False,)
        assert decode_status is True, "decode did not return True for matching parameters"


def test_decode_failure():
    from pgap2.lib.pklcheck import PklCheck
    import pickle

    pkl_file = 'tests/output/preprocess.pkl'
    with open(pkl_file, 'rb') as fh:
        previous: PklCheck = pickle.load(fh)
        decode_status = previous.decode(
            orth_id=0.98, para_id=0.7, dup_id=0.99,
            accurate=True, coverage=0.98,
            id_attr_key='ID', type_filter='CDS',
            LD=0.6, AS=0.6, AL=0.6, evalue=1e-5, aligner='diamond',
            clust_method='cdhit', falen=20, annot=False, retrieve=False,)
        assert decode_status is False, "decode did not return False for non-matching parameters"


def test_get_similarity():
    from pgap2.utils.partition import get_similarity

    similarity = get_similarity([1, 2, 3, 4, 5, 6], [1, 2, 3])
    assert similarity == 1, "get_similarity did not return the expected similarity value"

    similarity = get_similarity([1, 2, 3], [1, 2, 3, 4, 5, 6])
    assert similarity == 1, "get_similarity did not return the expected similarity value"

    similarity = get_similarity([1, 2, 3], [4, 5, 6])
    assert similarity == 0, "get_similarity did not return the expected similarity value"

    similarity = get_similarity([1, 2, 3, 4, 5], [5, 6])
    assert similarity == 0.5, "get_similarity did not return the expected similarity value"


def test_shortest_path_length_with_max_length():
    from pgap2.utils.partition import shortest_path_length_with_max_length
    import networkx as nx
    G = nx.Graph()
    G.add_edge(1, 2, weight=1)
    G.add_edge(2, 3, weight=1)
    G.add_edge(3, 4, weight=1)
    G.add_edge(4, 5, weight=1)
    a_adj, path = shortest_path_length_with_max_length(
        G, 1, 5, {}, 5)
    assert path == [
        1, 2, 3, 4, 5], "should return the expected path from 1 to 5"

    assert a_adj == {1: {2}, 2: {1, 3}, 3: {2, 4}, 4: {
        3, 5}}, "should return the expected adjacency list"

    a_adj, path = shortest_path_length_with_max_length(
        G, 1, 5, a_adj, 5)
    assert path == [
        1, 2, 3, 4, 5], "should return the same path when max length is not exceeded"

    a_adj, path = shortest_path_length_with_max_length(
        G, 1, 5, a_adj, 3)
    assert path == [], "should return an empty path when max length is exceeded"


def test_similarity_partition():
    from pgap2.utils.partition import similarity_partition
    from pgap2.lib.tree import Tree
    import networkx as nx
    G = nx.Graph()
    G.add_node(1, repre_nodes=[1])
    G.add_node(2, repre_nodes=[2])
    G.add_node(3, repre_nodes=[3])
    G.add_node(4, repre_nodes=[4])
    G.add_node(5, repre_nodes=[5])

    G.add_edge(1, 2, weight=1)
    G.add_edge(2, 3, weight=1)
    G.add_edge(3, 4, weight=1)
    G.add_edge(4, 5, weight=1)
    tree = Tree()
    identity_G = nx.Graph()
    identity_G.add_edge(1, 5, weight=1)
    identity_G.add_edge(2, 4, weight=1)
    tree.distance_graph = identity_G

    need_merge_nodes, merge_node_attr = similarity_partition(
        tree, G, [1, 2, 3, 4, 5], 3, {}, {})
    assert need_merge_nodes == [[1, 5], [2, 4]
                                ], "should return the expected merge nodes"


def test_find_final_node():
    from pgap2.utils.tools import find_final_node

    mapping = {1: 1, 2: 1, 3: 2, 4: 4, 5: 5}
    assert find_final_node(
        1, mapping) == 1, "should return the same node if it points to itself"
    assert find_final_node(
        3, mapping) == 1, "should return the final node after following the mapping"
    assert find_final_node(
        2, mapping) == 1, "should return the final node after following multiple mappings"


def test_test_connectedness():
    from pgap2.utils.tools import test_connectedness
    from pgap2.lib.tree import Tree
    import networkx as nx
    tree = Tree()
    G = nx.Graph()
    # G.add_edge(1, 2)
    # G.add_edge(2, 3)
    # G.add_edge(3, 4)
    # G.add_edge(4, 5)
    G.add_node(1, repre_nodes=[1, 2, 3])
    G.add_node(5, repre_nodes=[4, 5])

    distance_graph = nx.Graph()
    distance_graph.add_node(1)
    distance_graph.add_node(2)
    distance_graph.add_node(3)
    distance_graph.add_node(4)
    distance_graph.add_node(5)
    distance_graph.add_edge(1, 5, weight=1)
    distance_graph.add_edge(1, 4, weight=1)
    distance_graph.add_edge(2, 4, weight=1)
    tree.distance_graph = distance_graph
    tree.update_distance_matrix()

    need_merge = test_connectedness(tree, G, 1, 5, 'moderate')
    assert need_merge == True, "sensitivity:moderate: should return True for connected nodes"

    need_merge = test_connectedness(tree, G, 1, 5, 'soft')
    assert need_merge == True, "sensitivity:soft: should return True for connected nodes with soft sensitivity"

    need_merge = test_connectedness(tree, G, 1, 5, 'strict')
    assert need_merge == False, "sensitivity:strict: should return False for connected nodes with strict sensitivity"


def test_merge_node():
    from pgap2.utils.tools import merge_node
    from pgap2.lib.tree import Tree
    import networkx as nx
    G = nx.Graph()
    tree = Tree()
    pg = None

    G.add_node(1, repre_nodes=[1, 2, 3], members=[0], strains=[0])
    G.add_node(5, repre_nodes=[4, 5], members=[6], strains=[6])

    G = merge_node(G, pg, tree, [1, 5], 5)
    assert not G.has_node(1), "should keep the merged node 5"
    assert G.has_node(5), "should remove the merged node 1"
    assert G.nodes[5]['repre_nodes'] == {1, 2, 3, 4,
                                         5}, "should merge the representative nodes correctly"
    assert G.nodes[5]['members'] == {
        0, 6}, "should merge the members correctly"
    assert G.nodes[5]['strains'] == {
        0, 6}, "should merge the strains correctly"


def test_retrieve_merge_sorted_lists():
    from pgap2.utils.gene_retriever import merge_sorted_lists
    list1 = [1, 4, 6, 7]
    list2 = [2, 3, 5]
    merged_list = merge_sorted_lists(list1, list2)
    assert merged_list == [
        1, 2, 3, 4, 5, 6, 7], "merge_sorted_lists did not return the expected merged list"


def test_retrieve_find_nearest_numbers():
    from pgap2.utils.gene_retriever import find_nearest_numbers
    numbers = [1, 3, 5, 7, 9]
    target = 6
    nearest_numbers = find_nearest_numbers(numbers, target)
    assert nearest_numbers == (
        5, 7), "find_nearest_numbers did not return the expected nearest numbers"
