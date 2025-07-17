import pytest


@pytest.fixture
def test_basic():
    basic_file = 'tests/output/basic.pkl'
    from pgap2.lib.pklcheck import PklCheck
    from pgap2.lib.basic import Basic
    import pickle

    with open(f'{basic_file}', 'rb') as fh:
        previous: PklCheck = pickle.load(fh)
        decode_status = previous.decode()
        assert decode_status is True, "decode did not return True for matching parameters"
        basic: Basic = previous.data_dump('basic')
        assert isinstance(
            basic, Basic), "Basic object was not created correctly"
        return basic


def test_get_pan_group(test_basic):
    from pgap2.postprocess.stat import get_pan_group
    basic = test_basic
    outdir = 'tests/output'
    basic.load_pav(file=f'{outdir}/pgap2.partition.gene_content.pav')
    pan_group_freq, pan_para_stat = get_pan_group(basic.pav, outdir)
    assert len(
        pan_group_freq) == 3, "get_pan_group did not return the expected number of groups"
    assert len(
        pan_para_stat) == 5, "get_pan_group did not return the expected number of parameters"


def test_rare(test_basic):
    from pgap2.postprocess.stat import get_rarefaction, fit_rerefaction
    basic = test_basic
    outdir = 'tests/output/'
    basic.load_pav(file=f'{outdir}/pgap2.partition.gene_content.pav')
    pan_profile, new_clusters = get_rarefaction(
        basic.pav, outdir, N=1)
    assert len(
        pan_profile) == 3, "get_rarefaction did not return the expected number of profiles"
    assert len(
        new_clusters) == 2, "get_rarefaction did not return the expected number of clusters"
