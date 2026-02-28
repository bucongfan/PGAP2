import os
import pickle
import networkx as nx

from tqdm import tqdm
from loguru import logger
from collections import OrderedDict

from pgap2.lib.tree import Tree
from pgap2.lib.basic import Basic
from pgap2.lib.strain import Strain
from pgap2.lib.pklcheck import PklCheck
from pgap2.lib.pangenome import Pangenome

from pgap2.utils.supply import tqdm_
from pgap2.utils.generate_tree import generate_tree
from pgap2.utils.arrangement_detector import merge_by_synteny
from pgap2.utils.gene_retriever import retrieve_gene

from pgap2.utils.partition import (
    generate_network, get_pan_clust, merge_by_similarity,
    mcl, get_expect_identity
)

"""
Merge mode for PGAP2.

This module merges two independent PGAP2 results into a single pangenome.

"""


# ─────────────────────── Reindex helpers ──────────────────────────────────────

def _reindex_gene_id(gene_id: str, offset: int) -> str:
    """
    Reindex a gene ID by adding offset to its strain index.

    Gene ID format: 'strain_index:contig_index:gene_index'
    """
    parts = gene_id.split(':')
    parts[0] = str(int(parts[0]) + offset)
    return ':'.join(parts)


def reindex_prot_file(src_path: str, dst_path: str, offset: int):
    """
    Rewrite a protein FASTA file with reindexed gene IDs.
    """
    with open(src_path, 'r') as fin, open(dst_path, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                gene_id = line[1:].strip()
                new_id = _reindex_gene_id(gene_id, offset)
                fout.write(f'>{new_id}\n')
            else:
                fout.write(line)


def reindex_annot_file(src_path: str, dst_path: str, offset: int):
    """
    Rewrite an annotation TSV file with reindexed gene IDs.
    The gene index is in the first column.
    """
    with open(src_path, 'r') as fin, open(dst_path, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue
            parts = line.split('\t', 1)
            if len(parts) < 2:
                fout.write(line)
                continue
            new_id = _reindex_gene_id(parts[0], offset)
            fout.write(f'{new_id}\t{parts[1]}')


def reindex_pangenome(pg_b: Pangenome, offset: int) -> Pangenome:
    """
    Reindex all strain indices in pg_b by adding offset.
    Returns a new Pangenome with shifted strain indices.
    """
    new_pg = Pangenome(
        outdir=pg_b.outdir,
        threads=pg_b.threads,
        gcode=pg_b.gcode,
        disable=pg_b.disable_tqdm
    )
    for old_idx, strain in pg_b.strain_dict.items():
        new_idx = old_idx + offset
        new_strain = Strain(
            strain_name=strain.strain_name,
            strain_index=new_idx,
            bed_gene_num=strain.bad_gene_num,
            gene_num=strain.gene_num
        )
        new_pg.load_strain(new_strain)

    new_pg.total_gene_num = pg_b.total_gene_num
    new_pg.orth_id = pg_b.orth_id
    new_pg.para_id = pg_b.para_id
    new_pg.dup_id = pg_b.dup_id
    new_pg.accurate = pg_b.accurate
    new_pg.exhaust_orth = pg_b.exhaust_orth
    new_pg.retrieve = pg_b.retrieve
    new_pg.evalue = pg_b.evalue
    new_pg.aligner = pg_b.aligner
    new_pg.LD = pg_b.LD
    new_pg.AL = pg_b.AL
    new_pg.AS = pg_b.AS
    return new_pg


def reindex_graph(G: nx.Graph, offset: int) -> nx.Graph:
    """
    Reindex all gene IDs in a pangenome graph (graph.pkl).
    Node names and member/strain sets are all reindexed.
    """
    relabel_map = {}
    for node in G.nodes():
        if '_' in node:
            parts = node.rsplit('_', 1)
            new_child = _reindex_gene_id(parts[1], offset)
            # The node name format after paralog split: 'parent_child'
            # where parent is already a gene_id
            new_parent = _reindex_gene_id(parts[0], offset)
            relabel_map[node] = f'{new_parent}_{new_child}'
        else:
            relabel_map[node] = _reindex_gene_id(node, offset)

    new_G = nx.relabel_nodes(G, relabel_map, copy=True)

    for node in new_G.nodes():
        attrs = new_G.nodes[node]
        if 'members' in attrs:
            attrs['members'] = {_reindex_gene_id(
                m, offset) for m in attrs['members']}
        if 'strains' in attrs:
            attrs['strains'] = {s + offset for s in attrs['strains']}
        if 'repre_nodes' in attrs:
            new_repre = []
            for rn in attrs['repre_nodes']:
                if '_' in rn:
                    parts = rn.rsplit('_', 1)
                    new_repre.append(
                        f'{_reindex_gene_id(parts[0], offset)}_{_reindex_gene_id(parts[1], offset)}')
                else:
                    new_repre.append(_reindex_gene_id(rn, offset))
            attrs['repre_nodes'] = new_repre

    return new_G


# ─────────────────────── Guide dict ──────────────────────────────────────────

def build_guide_dict(G: nx.Graph) -> dict:
    """
    Build a guide dictionary from a final graph.
    Maps each gene member to its node name (cluster ID).
    """
    guide_dict = {}
    for node in G.nodes():
        members = G.nodes[node].get('members', set())
        for member in members:
            guide_dict[member] = node
    logger.info(
        f'Built guide dict: {len(guide_dict)} genes → {G.number_of_nodes()} clusters')
    return guide_dict


def build_merged_guide_dict(G_a: nx.Graph, G_b: nx.Graph) -> dict:
    """
    Build a combined guide dictionary from two previous graphs.
    All genes from both graphs are mapped to their cluster names.
    Since G_b has already been reindexed, there should be no key conflicts.
    """
    guide_dict = {}
    for G in (G_a, G_b):
        for node in G.nodes():
            members = G.nodes[node].get('members', set())
            for member in members:
                guide_dict[member] = node
    logger.info(
        f'Built merged guide dict: {len(guide_dict)} genes from 2 graphs')
    return guide_dict


# ─────────────────────── Merge pangenome ─────────────────────────────────────

def merge_two_pangenomes(pg_a: Pangenome, pg_b: Pangenome, outdir: str,
                         threads: int, gcode: int, disable: bool) -> Pangenome:
    """
    Merge two pangenomes into one.
    pg_b should already be reindexed.
    """
    logger.info(
        f'Merging pangenomes: {pg_a.strain_num} + {pg_b.strain_num} strains')

    merged_pg = Pangenome(
        outdir=outdir,
        threads=threads,
        gcode=gcode,
        disable=disable
    )

    for strain in pg_a.strain_dict.values():
        merged_pg.load_strain(strain)
    for strain in pg_b.strain_dict.values():
        merged_pg.load_strain(strain)

    merged_pg.total_gene_num = pg_a.total_gene_num + pg_b.total_gene_num

    # Copy parameters from set A (they should be consistent)
    for attr in ('orth_id', 'para_id', 'dup_id', 'accurate', 'exhaust_orth',
                 'retrieve', 'evalue', 'aligner', 'LD', 'AL', 'AS'):
        setattr(merged_pg, attr, getattr(pg_a, attr))

    logger.info(
        f'Merged pangenome: {merged_pg.strain_num} strains, {merged_pg.total_gene_num} genes')
    return merged_pg


# ─────────────────── Cross-set alignment + tree building ─────────────────────

def build_merged_tree(outdir: str, orth_list: list, coverage: float,
                      evalue: float, falen: int, threads: int,
                      max_targets: int, para_id: float, LD: float,
                      AS: float, AL: float, aligner: str,
                      clust_method: str, disable: bool) -> Tree:
    """
    Build a merged Tree by running generate_tree from scratch on the
    combined protein file.

    Trying to manually merge two orth_trees and two distance graphs is
    extremely error-prone: the internal data structures (leaf_root,
    member_leaf, ortho_para, etc.) become inconsistent because the merged
    orth_tree roots and the distance graph nodes don't match up properly.

    Instead we simply re-cluster and re-align all representatives from
    scratch.  Clustering (mmseqs2/cd-hit) is fast, and alignment only
    runs on the final representatives, so the cost is acceptable.
    """
    merged_prot = f'{outdir}/total.involved_prot.fa'
    logger.info('Building merged tree from combined protein file...')
    tree = generate_tree(
        input_file=merged_prot, orth_list=orth_list, outdir=outdir,
        evalue=evalue, aligner=aligner, falen=falen, disable=disable,
        threads=threads, max_targets=max_targets, coverage=coverage,
        ID=para_id, LD=LD, AS=AS, AL=AL, clust_method=clust_method)
    return tree


# ─────────────────────── Check duplicate strains ─────────────────────────────

def check_duplicate_strains(pg_a: Pangenome, pg_b: Pangenome):
    """
    Check for duplicate strain names between two pangenomes.
    Returns list of duplicate strain names.
    """
    names_a = {s.strain_name for s in pg_a.strain_dict.values()}
    names_b = {s.strain_name for s in pg_b.strain_dict.values()}
    duplicates = names_a & names_b
    if duplicates:
        logger.warning(
            f'Found {len(duplicates)} overlapping strain(s) between the two sets:')
        for dup in sorted(duplicates):
            logger.warning(f'  - {dup}')
        logger.error(
            'Cannot merge results with overlapping strain names. '
            'Please remove duplicates from one side first.')
        raise ValueError(
            f'Overlapping strain names found: {sorted(duplicates)}')
    return duplicates


# ─────────────────────── Main ─────────────────────────────────────────────────

def main(dir_a: str, dir_b: str, outdir: str, aligner: str, clust_method: str,
         falen: int, threads: int, id_attr_key: str, type_filter: str,
         annot: bool, gcode: int, retrieve: bool, disable: bool):

    # ── Load parameters from set A ──
    logger.info(f'Loading parameters from {dir_a}/basic.pkl')
    with open(f'{dir_a}/basic.pkl', 'rb') as fh:
        pkl_a: PklCheck = pickle.load(fh)
        basic_a = pkl_a.data_dump('basic')
        params = basic_a.params

    # Inherit parameters from set A
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

    # ── Log parameters ──
    logger.info('Inherited parameters from set A:')
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
        logger.info(f'  {k}: {v}')

    # ── Check parameter consistency with set B ──
    logger.info(f'Loading parameters from {dir_b}/basic.pkl')
    with open(f'{dir_b}/basic.pkl', 'rb') as fh:
        pkl_b: PklCheck = pickle.load(fh)
        basic_b = pkl_b.data_dump('basic')
        params_b = basic_b.params

    if params_b:
        critical_keys = ['orth_id', 'para_id',
                         'dup_id', 'evalue', 'LD', 'AS', 'AL']
        for key in critical_keys:
            if key in params_b and key in params and params[key] != params_b[key]:
                logger.warning(
                    f"Parameter mismatch between set A and B: {key} "
                    f"(A: {params[key]}, B: {params_b[key]})")

    # ── Load pangenomes ──
    logger.info('Loading pangenomes from preprocess.pkl...')
    with open(f'{dir_a}/preprocess.pkl', 'rb') as fh:
        pkl_pre_a: PklCheck = pickle.load(fh)
        pg_a = pkl_pre_a.data_dump('pangenome')

    with open(f'{dir_b}/preprocess.pkl', 'rb') as fh:
        pkl_pre_b: PklCheck = pickle.load(fh)
        pg_b = pkl_pre_b.data_dump('pangenome')

    # ── Check for duplicate strain names ──
    check_duplicate_strains(pg_a, pg_b)

    n_a = pg_a.strain_num
    n_b = pg_b.strain_num
    reindex_offset = n_a
    logger.info(f'Set A: {n_a} strains, Set B: {n_b} strains')
    logger.info(f'Reindex offset for set B: {reindex_offset}')

    # ── Reindex set B ──
    logger.info('Reindexing set B...')
    pg_b_reindexed = reindex_pangenome(pg_b, reindex_offset)

    # ── Create output directory ──
    os.makedirs(outdir, exist_ok=True)

    # ── Merge protein and annotation files ──
    logger.info('Merging protein and annotation files...')

    prot_a = f'{dir_a}/total.involved_prot.fa'
    annot_a = f'{dir_a}/total.involved_annot.tsv'
    prot_b = f'{dir_b}/total.involved_prot.fa'
    annot_b = f'{dir_b}/total.involved_annot.tsv'

    # Reindex B files
    prot_b_reindexed = f'{outdir}/prot_b_reindexed.fa'
    annot_b_reindexed = f'{outdir}/annot_b_reindexed.tsv'
    reindex_prot_file(prot_b, prot_b_reindexed, reindex_offset)
    reindex_annot_file(annot_b, annot_b_reindexed, reindex_offset)

    # Merge into final files
    merged_prot = f'{outdir}/total.involved_prot.fa'
    merged_annot = f'{outdir}/total.involved_annot.tsv'

    with open(merged_prot, 'w') as out_fh:
        with open(prot_a, 'r') as fh:
            out_fh.write(fh.read())
        with open(prot_b_reindexed, 'r') as fh:
            out_fh.write(fh.read())

    with open(merged_annot, 'w') as out_fh:
        with open(annot_a, 'r') as fh:
            out_fh.write(fh.read())
        with open(annot_b_reindexed, 'r') as fh:
            # Skip the header line of set B
            for i, line in enumerate(fh):
                if i == 0 and line.startswith('#'):
                    continue
                out_fh.write(line)

    # ── Merge pangenomes ──
    pg = merge_two_pangenomes(pg_a, pg_b_reindexed,
                              outdir, threads, gcode, disable)
    pg.load_annot_file(merged_annot)
    pg.load_prot_file(merged_prot)

    # ── Copy genome_index directories for retrieve mode ──
    if retrieve:
        import shutil
        genome_idx_out = f'{outdir}/genome_index'
        os.makedirs(genome_idx_out, exist_ok=True)
        # Copy A's genome_index
        genome_idx_a = f'{dir_a}/genome_index'
        if os.path.exists(genome_idx_a):
            for item in os.listdir(genome_idx_a):
                src = os.path.join(genome_idx_a, item)
                dst = os.path.join(genome_idx_out, item)
                if os.path.isdir(src) and not os.path.exists(dst):
                    shutil.copytree(src, dst)
        # Copy B's genome_index (with reindexed directory names)
        genome_idx_b = f'{dir_b}/genome_index'
        if os.path.exists(genome_idx_b):
            for dir_idx in os.listdir(genome_idx_b):
                src_dir = os.path.join(genome_idx_b, dir_idx)
                if not os.path.isdir(src_dir):
                    continue
                for strain_idx in os.listdir(src_dir):
                    try:
                        old_idx = int(strain_idx)
                        new_idx = old_idx + reindex_offset
                        new_dir_idx = new_idx // 1000
                        dst_parent = os.path.join(
                            genome_idx_out, str(new_dir_idx))
                        os.makedirs(dst_parent, exist_ok=True)
                        src = os.path.join(src_dir, strain_idx)
                        dst = os.path.join(dst_parent, str(new_idx))
                        if not os.path.exists(dst):
                            shutil.copytree(src, dst)
                    except ValueError:
                        continue

    # ── Build merged tree ──
    logger.info('Building merged tree...')
    orth_list = [dup_id, orth_id]
    tree = build_merged_tree(
        outdir=outdir, orth_list=orth_list,
        coverage=coverage, evalue=evalue, falen=falen,
        threads=threads, max_targets=max_targets,
        para_id=para_id, LD=LD, AS=AS, AL=AL,
        aligner=aligner, clust_method=clust_method,
        disable=disable
    )

    # ── Save preprocess checkpoint ──
    logger.info('Saving preprocess checkpoint...')
    file_dict_merged = OrderedDict()
    for strain in pg.strain_dict.values():
        file_dict_merged[strain.strain_name] = {}

    pickle_preprocess = PklCheck(outdir=outdir, name='preprocess')
    pickle_preprocess.load('file_dict', main_data=file_dict_merged)
    pickle_preprocess.load('pangenome', main_data=pg, parameter={
        'orth_id': orth_id, 'para_id': para_id, 'dup_id': dup_id,
        'accurate': accurate, 'id_attr_key': id_attr_key,
        'type_filter': type_filter, 'coverage': coverage,
        'AS': AS, 'AL': AL, 'LD': LD, 'retrieve': retrieve,
        'evalue': evalue, 'aligner': aligner, 'clust_method': clust_method,
        'annot': annot, 'falen': falen})
    pickle_preprocess.load('tree', main_data=tree)
    pickle_preprocess.pickle_()

    # ── Build guide dict from both previous graphs ──
    logger.info('Loading previous graphs for guide construction...')
    with open(f'{dir_a}/graph.pkl', 'rb') as fh:
        G_a: PklCheck = pickle.load(fh)
        G_a = G_a.data_dump('graph')

    with open(f'{dir_b}/graph.pkl', 'rb') as fh:
        G_b: PklCheck = pickle.load(fh)
        G_b = G_b.data_dump('graph')

    # Reindex G_b
    G_b_reindexed = reindex_graph(G_b, reindex_offset)

    # Build combined guide dict
    guide_dict = build_merged_guide_dict(G_a, G_b_reindexed)

    # Mark set B strains as "new" for the guide logic.
    # This ensures:
    # - Within A pairs: all has_new=False → guide decides instantly (fast path)
    # - Within B pairs: all has_new=True → guide returns None → normal judgment
    # - Cross-set pairs (A+B): mixed has_new → guide returns None → normal judgment
    # Cross-set pairs MUST use normal judgment because genes from different graphs
    # have different guide cluster names and would be incorrectly rejected by the
    # guide. Marking B as "new" forces _guide_says_merge to return None for any
    # pair involving B genes, letting the standard merge logic handle them.
    # B's internal structure is already well-clustered from its previous run,
    # so the normal judgment overhead for within-B pairs is minimal.
    new_strain_indices = set(range(reindex_offset, reindex_offset + n_b))
    logger.info(
        f'New strain indices (set B): {min(new_strain_indices)}-{max(new_strain_indices)}')

    # ── Set pangenome parameters ──
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

    # ─────────────────────── Partition step ──────────────────────────────────
    logger.info('Get the gene primal clust result by mcl')
    mcl(pg, tree)
    logger.info('Load the gene length information')
    pg.reload_nucl_file(tree)
    logger.info('Create synteny network')
    G, tree = generate_network(
        pg=pg, tree=tree, guide_dict=guide_dict, new_strain_indices=new_strain_indices)

    tree.load_para_id(para_id)
    tree.load_orth_id(orth_id)
    tree.load_dup_id(dup_id)
    logger.info('Build index')
    max_in_range = get_expect_identity(tree, G, pg)

    logger.info(f'Load expect identity: {max_in_range}')
    tree.load_expect_identity(max_in_range)
    logger.info('Clean up the distance graph according to paralogous genes')
    tree.update_distance_graph(disable=disable)
    logger.info('Merge by gene similarity')

    G, pg, tree = merge_by_similarity(G=G, pg=pg, tree=tree,
                                      sensitivity=sensitivity,
                                      radius=radius, fast=fast_mode,
                                      context_sim=context_similarity,
                                      flank=flank,
                                      disable=disable,
                                      guide_dict=guide_dict,
                                      new_strain_indices=new_strain_indices)

    logger.info('Double check through gene synteny')
    G = merge_by_synteny(G, pg, tree,
                         context_sim=context_similarity,
                         flank=flank,
                         sensitivity=sensitivity,
                         ins=ins,
                         guide_dict=guide_dict,
                         new_strain_indices=new_strain_indices)

    logger.info('Reload the gene annotation')
    pg.reload_annot_file(retrieve=retrieve)

    if retrieve:
        logger.info('Retrieve gene from missing')
        G, pg, tree = retrieve_gene(G, pg, tree)
        G, pg, tree = merge_by_similarity(G=G, pg=pg, tree=tree,
                                          sensitivity=sensitivity,
                                          radius=radius, fast=fast_mode,
                                          context_sim=context_similarity,
                                          flank=flank,
                                          disable=disable, step=10,
                                          guide_dict=guide_dict,
                                          new_strain_indices=new_strain_indices)
        G = merge_by_synteny(G=G, pg=pg, tree=tree,
                             context_sim=context_similarity,
                             flank=flank,
                             sensitivity=sensitivity,
                             ins=ins, step=11,
                             guide_dict=guide_dict,
                             new_strain_indices=new_strain_indices)

    # ─────────────────────── Organize results ────────────────────────────────
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

    pickle_G = PklCheck(outdir=outdir, name='graph')
    pickle_G.load('graph', main_data=G)
    pickle_G.pickle_()

    logger.info('Save basic results for downstream analysis...')
    pickle_basic = PklCheck(outdir=outdir, name='basic')
    pickle_basic.load('basic', main_data=Basic(pg=pg, params=params))
    pickle_basic.pickle_()

    return 0
