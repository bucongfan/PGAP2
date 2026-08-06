import pytest
import numpy as np


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


def test_fit_rarefaction_keeps_running_when_core_fit_fails(tmp_path, monkeypatch):
    from pgap2.postprocess import stat

    pan_profile = [
        [1, 2, 3, 4],
        [[3000], [2500], [2300], [2200]],
        [[4000], [4500], [4800], [5100]],
    ]
    original_curve_fit = stat.curve_fit

    def fail_core_fit(model, *args, **kwargs):
        if model is stat.core_gene_model:
            raise RuntimeError("forced core failure")
        return original_curve_fit(model, *args, **kwargs)

    monkeypatch.setattr(stat, "curve_fit", fail_core_fit)

    stat.fit_rerefaction(pan_profile, tmp_path)

    content = (tmp_path / "postprocess.curve_fit.txt").read_text()
    assert "Optimal parameters for core genes: NA" in content
    assert "Pan Gene Fit: y=" in content


def test_fit_rarefaction_uses_pan_loglog_fallback(tmp_path, monkeypatch):
    from pgap2.postprocess import stat

    pan_profile = [
        [1, 2, 4, 8],
        [[3000], [2500], [2300], [2200]],
        [[1000], [1414], [2000], [2828]],
    ]

    def fail_pan_fit(model, *args, **kwargs):
        if model is stat.pan_gene_model:
            raise RuntimeError("forced pan failure")
        return (
            np.asarray([800.0, 0.1, 2200.0], dtype=float),
            np.eye(3),
        )

    monkeypatch.setattr(stat, "curve_fit", fail_pan_fit)

    stat.fit_rerefaction(pan_profile, tmp_path)

    content = (tmp_path / "postprocess.curve_fit.txt").read_text()
    assert "Optimal parameters for pan genes:" in content
    assert "log-log fallback" in content
    assert "Pan Gene Fit: y=" in content


@pytest.mark.parametrize(
    ("alpha", "expected_state"),
    [
        (0.5, "open"),
        (1.5, "closed"),
    ],
)
def test_fit_rarefaction_classifies_state_from_new_genes(
        tmp_path, alpha, expected_state):
    from pgap2.postprocess import stat

    x_data = np.asarray([1, 2, 4, 8, 16, 32], dtype=float)
    pan_profile = [
        x_data.tolist(),
        [[3000 - 20 * i] for i in range(len(x_data))],
        [[4000 + 300 * np.sqrt(x)] for x in x_data],
    ]
    new_x_data = x_data[1:]
    new_gene_means = 1000 * new_x_data**(-alpha)
    new_clusters = [
        new_x_data.tolist(),
        [[value] for value in new_gene_means],
    ]

    stat.fit_rerefaction(pan_profile, tmp_path, new_clusters)

    content = (tmp_path / "postprocess.curve_fit.txt").read_text()
    assert "New Gene Fit: y=" in content
    assert f"The pan-genome is {expected_state}" in content


def test_fit_rarefaction_reports_na_with_too_few_new_gene_points(tmp_path):
    from pgap2.postprocess import stat

    pan_profile = [
        [1, 2, 3, 4],
        [[3000], [2500], [2300], [2200]],
        [[4000], [4500], [4800], [5100]],
    ]
    new_clusters = [[2, 3], [[500], [200]]]

    stat.fit_rerefaction(pan_profile, tmp_path, new_clusters)

    content = (tmp_path / "postprocess.curve_fit.txt").read_text()
    assert "Optimal parameters for new genes: NA" in content
    assert "Pan-genome state: NA" in content
