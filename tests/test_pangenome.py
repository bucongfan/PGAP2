from pgap2.lib.pangenome import Pangenome


def test_reload_annot_file_maps_same_contig_name_per_strain(tmp_path):
    annot_file = tmp_path / "total.involved_annot.tsv"
    annot_file.write_text(
        "#Gene_index\tStrain\tContig\tLocation\tLength\tGene_ID\tGene_name\tProduct_name\tNucleotide_sequence\tProtein_sequence\n"
        "0:3:0\tstrain_a\tcontig_1\t[1:90](+)\t30\tgene_a\t\tproduct\tATG\tM\n"
        "62:3:0\tstrain_b\tcontig_1\t[5:94](+)\t30\tgene_b\t\tproduct\tATG\tM\n",
        encoding="utf-8",
    )

    pg = Pangenome(outdir=str(tmp_path), threads=1, gcode=11, disable=True)
    pg.total_gene_num = 2
    pg.load_annot_file(str(annot_file))
    pg.reload_annot_file(retrieve=True)

    assert pg.annot_contig_map["0:3"] == "contig_1"
    assert pg.annot_contig_map["62:3"] == "contig_1"
