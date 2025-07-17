import pytest


@pytest.fixture
def decode_pkl():
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
        pg = previous.data_dump('pangenome')
        tree = previous.data_dump('tree')
        return pg, tree


def test_stat_core(decode_pkl):
    from pgap2.utils.preprocess import stat_core
    from collections import defaultdict

    pg, tree = decode_pkl
    core_clusters = stat_core(tree, pg)
    assert isinstance(
        core_clusters, tuple), "stat_core did not return a tuple"
    assert len(core_clusters) == 2, "stat_core returned an empty tuple"


def test_genome_feature_stat(decode_pkl):
    from pgap2.utils.preprocess import genome_feature_stat
    from pgap2.lib.pangenome import Pangenome
    from pgap2.lib.strain import Strain
    from pgap2.lib.species import Species

    pg, tree = decode_pkl
    outdir = 'tests/output'
    sp = Species(marker_file=None,
                 strain_dict=pg.strain_dict, ani=95, outdir=outdir)
    assert isinstance(sp, Species), "Species object was not created correctly"

    for strain_index in pg.strain_dict:
        genome_record = pg.strain_dict[strain_index].genome
        chr_count, genome_len, atcg_dict = genome_feature_stat(
            genome_record)
        assert isinstance(chr_count, int), "chr_count is not an integer"
        assert isinstance(genome_len, int), "genome_len is not an integer"
        assert isinstance(atcg_dict, dict), "atcg_dict is not a dictionary"
        assert all(isinstance(v, int) for v in atcg_dict.values()
                   ), "atcg_dict values are not all integers"
        assert len(
            atcg_dict) == 5, "atcg_dict does not contain exactly 5 keys (A, T, C, G, N)"
