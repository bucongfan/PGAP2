from types import SimpleNamespace
from unittest.mock import Mock

import networkx as nx
import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from pgap2.lib.pangenome import Pangenome
from pgap2.lib.tree import Tree
from pgap2.utils import gene_retriever


def _paf_line(query, protein_length, matches, block, strand="+", score=50):
    fields = [
        query, protein_length, 0, min(protein_length, block // 3), strand,
        "0", 1200, 200, 200 + block, matches, block, 60, f"AS:i:{score}",
    ]
    return "\t".join(map(str, fields)) + "\n"


@pytest.mark.parametrize(
    "protein_length,matches,block,identity",
    [
        pytest.param(100, 120, 300, 0.4, id="full-length-low-identity"),
        pytest.param(100, 270, 300, 0.9, id="full-length-mismatches"),
        pytest.param(100, 270, 360, 0.75, id="longer-alignment"),
        pytest.param(100, 180, 240, 0.75, id="partial-alignment"),
        pytest.param(100, 300, 300, 1.0, id="perfect-match"),
        pytest.param(100, 198, 297, 0.66667, id="rounding"),
    ],
)
def test_retrieve_from_paf_identity(tmp_path, protein_length, matches, block,
                                   identity):
    paf_file = tmp_path / "hits.paf"
    paf_file.write_text(_paf_line("query", protein_length, matches, block))

    hits = gene_retriever.retrieve_from_paf(paf_file)

    assert hits["0"]["query"]["identity"] == identity


@pytest.mark.parametrize("high_identity_first", [False, True])
def test_retrieve_from_paf_overlap_prefers_identity(tmp_path, high_identity_first):
    lines = [
        _paf_line("low_identity", 100, 120, 300),
        _paf_line("high_identity", 100, 270, 300),
    ]
    if high_identity_first:
        lines.reverse()
    paf_file = tmp_path / "hits.paf"
    paf_file.write_text("".join(lines))

    hits = gene_retriever.retrieve_from_paf(paf_file)

    assert set(hits["0"]) == {"high_identity"}
    assert hits["0"]["high_identity"]["identity"] == 0.9


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize(
    "protein_length,matches,outcome",
    [
        pytest.param(100, 120, "rejected", id="low-identity-filtered"),
        pytest.param(100, 270, "inserted", id="below-orth-not-merged"),
        pytest.param(120, 270, "inserted", id="mixed-strand-neighbors"),
        pytest.param(100, 300, "merged", id="high-identity-merged"),
    ],
)
def test_retrieval_identity_and_neighbors(tmp_path, monkeypatch, strand,
                                         protein_length, matches, outcome):
    expected_node = "1:0:0"
    retrieved_gene = "0:0:4"
    left_node, right_node = "0:0:1", "0:0:3"
    graph = nx.Graph()
    graph.add_node(expected_node, length=protein_length, members={expected_node},
                   strains={1}, repre_nodes=[expected_node])
    # The same-strand flanks are genes 1 and 3, not strand-local indices 0 and 1.
    positions = [(0, -strand), (60, strand), (90, -strand), (900, strand)]
    rows = [
        f"{expected_node}\treference\tcontig_1\t[0:360](+)\t{protein_length}"
        "\treference_gene\tgene\tproduct\tATG\tM\n"
    ]
    for index, (start, gene_strand) in enumerate(positions):
        gene = f"0:0:{index}"
        graph.add_node(gene, length=5, members={gene, f"1:0:{index + 1}"},
                       strains={0, 1}, repre_nodes=[gene])
        sign = "+" if gene_strand == 1 else "-"
        rows.append(
            f"{gene}\ttarget\tcontig_1\t[{start}:{start + 18}]({sign})"
            f"\t5\tgene_{index}\tgene\tproduct\tATG\tM\n")
    graph.add_edges_from([(expected_node, left_node), (expected_node, right_node)])

    pg = Pangenome(outdir=str(tmp_path), threads=1, gcode=11, disable=True)
    annot_file = tmp_path / "total.involved_annot.tsv"
    annot_file.write_text("".join(rows))
    pg.total_gene_num = len(rows)
    pg.load_annot_file(str(annot_file))
    pg.reload_annot_file(retrieve=True)
    pg.para_id, pg.orth_id = 0.7, 0.98
    pg.nucl_fa = {}
    pg.prot_file = str(tmp_path / "query-source.fa")
    pg.strain_dict = {0: SimpleNamespace(strain_name="target")}

    tree = Tree()
    tree.leaf_root[expected_node] = 0
    tree.root_leaf[0].add(expected_node)
    tree.distance_graph.add_node(expected_node)
    tree.raw_distance_graph = tree.distance_graph.copy()

    genome_dir = tmp_path / "genome_index" / "0" / "0"
    genome_dir.mkdir(parents=True)
    dna = "ATG" + "AAA" * 99
    genomic_dna = dna if strand == 1 else str(Seq(dna).reverse_complement())
    (genome_dir / "ref.fa").write_text(
        ">0\n" + "A" * 200 + genomic_dna + "A" * 700 + "\n")
    (genome_dir / "out.paf").write_text(
        _paf_line(expected_node, protein_length, matches, 300,
                  strand="+" if strand == 1 else "-"))
    monkeypatch.setattr(gene_retriever, "sfw",
                        SimpleNamespace(seqtk="seqtk", miniprot="miniprot"))
    monkeypatch.setattr(gene_retriever, "run_command", Mock())

    # Keep the PAF parser, graph updates, extraction, and annotation writing real.
    graph, pg, tree, realign_roots = gene_retriever.retrieve_gene_from_genome(
        graph, pg, tree)

    annotation = (tmp_path / "total.involved_retrieved_annot.tsv").read_text()
    if outcome == "rejected":
        assert annotation == ""
        assert retrieved_gene not in pg.annot
        assert not graph.has_node(retrieved_gene)
        assert graph.nodes[expected_node]["members"] == {expected_node}
        assert not tree.distance_graph.has_node(retrieved_gene)
        assert realign_roots == set()
    else:
        fields = annotation.strip().split("\t")
        assert fields[0] == "#" + retrieved_gene
        assert fields[7].split("|")[:2] == ["coding", str(matches / 300)]
        assert fields[8] == dna
        assert tree.distance_graph[retrieved_gene][expected_node]["weight"] == matches / 300
        assert realign_roots == {0}
        if outcome == "inserted":
            assert set(graph.neighbors(retrieved_gene)) == {left_node, right_node}
            assert graph.nodes[expected_node]["members"] == {expected_node}
        else:
            assert not graph.has_node(retrieved_gene)
            assert retrieved_gene in graph.nodes[expected_node]["members"]
            assert set(graph.neighbors(expected_node)) == {left_node, right_node}


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize(
    "dna,gcode,protein,status",
    [
        pytest.param("ATGAAATAAGGG", 11, "MK*G", "pseudo", id="internal-taa"),
        pytest.param("ATGAAATAGGGG", 11, "MK*G", "pseudo", id="internal-tag"),
        pytest.param("ATGAAATGAGGG", 11, "MK*G", "pseudo", id="internal-tga"),
        pytest.param("ATGAAAGGGTAA", 11, "MKG*", "coding", id="terminal-stop"),
        pytest.param("ATGAAAGGG", 11, "MKG", "coding", id="no-terminal-stop"),
        pytest.param("AAAGGG", 11, "KG", "coding", id="no-start-codon"),
        pytest.param("ATGAAATAAGGGTAA", 11, "MK*G*", "pseudo", id="both-stops"),
        pytest.param("ATGAAATAATAA", 11, "MK**", "pseudo", id="consecutive-stops"),
        pytest.param("ATGAAAGGGA", 11, "MKG", "pseudo", id="partial-codon"),
        pytest.param("ATGAAATGAGGG", 4, "MKWG", "coding", id="alternative-code"),
    ],
)
def test_retrieved_gene_classification(tmp_path, monkeypatch, strand, dna, gcode,
                                       protein, status):
    expected_node = "1:0:0"
    left_node, right_node = "0:0:0", "0:0:1"
    retrieved_gene = "0:0:2"
    graph = nx.Graph()
    graph.add_node(expected_node, length=10, members={expected_node},
                   strains={1}, repre_nodes=[expected_node])
    for node, other_member in [(left_node, "1:0:1"), (right_node, "1:0:2")]:
        graph.add_node(node, length=5, members={node, other_member},
                       strains={0, 1}, repre_nodes=[node])
        graph.add_edge(expected_node, node)

    pg = Pangenome(outdir=str(tmp_path), threads=1, gcode=gcode, disable=True)
    pg.para_id, pg.orth_id = 0.7, 0.98
    pg.nucl_fa = {}
    pg.prot_file = str(tmp_path / "query-source.fa")
    pg.annot = {
        expected_node: {"len": 10, "id": "reference_gene"},
        left_node: {"len": 5},
        right_node: {"len": 5},
    }
    pg.annot_contig_map = {"0:0": "contig_1"}
    pg.gene_rank = {"0:0": [[0, 300], []] if strand == 1 else [[], [0, 300]]}
    pg.strain_dict = {0: SimpleNamespace(strain_name="target_strain")}

    tree = Tree()
    tree.leaf_root[expected_node] = 0
    tree.root_leaf[0].add(expected_node)
    tree.distance_graph.add_node(expected_node)
    tree.raw_distance_graph = tree.distance_graph.copy()

    genome_dir = tmp_path / "genome_index" / "0" / "0"
    genome_dir.mkdir(parents=True)
    genomic_dna = dna if strand == 1 else str(Seq(dna).reverse_complement())
    (genome_dir / "ref.fa").write_text(
        ">0\n" + "A" * 100 + genomic_dna + "A" * 300 + "\n",
        encoding="utf-8",
    )

    # Supply a controlled hit while exercising the real extraction and output path.
    monkeypatch.setattr(gene_retriever, "sfw",
                        SimpleNamespace(seqtk="seqtk", miniprot="miniprot"))
    monkeypatch.setattr(gene_retriever, "run_command", Mock())
    monkeypatch.setattr(gene_retriever, "retrieve_from_paf", Mock(return_value={
        "0": {expected_node: {
            "identity": 0.99,
            "cigar": "cg:Z:4M",
            "loc": {"start": 100, "end": 100 + len(dna), "strand": strand},
        }},
    }))

    graph, pg, tree, realign_roots = gene_retriever.retrieve_gene_from_genome(
        graph, pg, tree)

    rows = (tmp_path / "total.involved_retrieved_annot.tsv").read_text().splitlines()
    assert len(rows) == 1
    fields = rows[0].split("\t")
    assert fields[0] == "#" + retrieved_gene
    assert fields[7].split("|")[0] == status
    assert fields[8] == dna
    assert fields[9] == protein
    assert int(fields[4]) == len(protein)
    assert pg.annot[retrieved_gene]["product"] == fields[7]
    assert retrieved_gene in graph.nodes[expected_node]["members"]
    if status == "pseudo":
        assert retrieved_gene in tree.pseudo_root_leaf[0]
        assert retrieved_gene not in pg.nucl_fa
        assert retrieved_gene not in tree.root_leaf[0]
        assert realign_roots == set()
    else:
        assert retrieved_gene not in tree.pseudo_root_leaf[0]
        assert str(pg.nucl_fa[retrieved_gene]) == dna
        assert retrieved_gene in tree.root_leaf[0]
        assert realign_roots == {0}


@pytest.mark.parametrize("gcode,protein", [(4, "MKW"), (11, "MK*")])
def test_realign_uses_configured_genetic_code(tmp_path, monkeypatch, gcode, protein):
    pg = Pangenome(outdir=str(tmp_path), threads=1, gcode=gcode, disable=True)
    pg.nucl_fa = {"0:0:0": Seq("ATGAAATGA")}
    tree = Tree()
    tree.root_leaf[0].update(["0:0:0", "1:0:0", "2:0:0"])
    parser = Mock(return_value=(None, None))
    monkeypatch.setattr(gene_retriever, "realign_parser", parser)
    monkeypatch.setattr(gene_retriever, "sfw",
                        SimpleNamespace(mcxdeblast="mcxdeblast", mcl="mcl"))

    gene_retriever.realign(pg, tree, str(tmp_path), {0})

    record = SeqIO.read(tmp_path / "seq.fa", "fasta")
    assert record.id == "0:0:0"
    assert str(record.seq) == protein
    parser.assert_called_once()
