import os
import pickle
import networkx as nx

from tqdm import tqdm
from loguru import logger
from collections import OrderedDict

from pgap2.lib.tree import Tree
from pgap2.lib.basic import Basic
from pgap2.lib.pklcheck import PklCheck
from pgap2.lib.pangenome import Pangenome

from pgap2.utils.supply import tqdm_
from pgap2.utils.arrangement_detector import merge_by_synteny
from pgap2.utils.gene_retriever import retrieve_gene

from pgap2.utils.partition import (
    generate_network, get_pan_clust, merge_by_similarity,
    mcl, get_expect_identity
)

"""
Minus mode for PGAP2.

This module removes specified strains from an existing PGAP2 result and
re-partitions the remaining genomes.
"""


# ─────────────────────── Load removal list ───────────────────────────────────

def load_removal_list(filepath: str) -> set:
    """
    Load strain names to remove from a text file.
    One strain name per line.  Blank lines and lines starting with '#' are
    skipped.
    """
    names = set()
    with open(filepath, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            names.add(line)
    if not names:
        raise ValueError(f'No strain names found in {filepath}')
    logger.info(
        f'Loaded {len(names)} strain name(s) to remove from {filepath}')
    return names


def resolve_strain_indices(pg: Pangenome, names_to_remove: set) -> set:
    """
    Map strain names to their integer indices.
    Raises an error if any name is not found.
    """
    name_to_idx = {}
    for idx, strain in pg.strain_dict.items():
        name_to_idx[strain.strain_name] = idx

    indices = set()
    missing = []
    for name in names_to_remove:
        if name in name_to_idx:
            indices.add(name_to_idx[name])
        else:
            missing.append(name)

    if missing:
        logger.error(f'The following strain names were not found in the '
                     f'previous result: {sorted(missing)}')
        available = sorted(name_to_idx.keys())
        logger.info(f'Available strains ({len(available)}): '
                    f'{available[:20]}{"..." if len(available) > 20 else ""}')
        raise ValueError(f'Unknown strain names: {sorted(missing)}')

    remaining = pg.strain_num - len(indices)
    if remaining < 2:
        raise ValueError(
            f'Cannot remove {len(indices)} strain(s) — only {remaining} '
            f'would remain (minimum 2 required)')

    logger.info(f'Resolved {len(indices)} strain index(es) to remove: '
                f'{sorted(indices)}')
    return indices


# ─────────────────────── Prune helpers ───────────────────────────────────────

def _gene_belongs_to_removed(gene_id: str, removed_indices: set) -> bool:
    """Check if a gene belongs to a removed strain."""
    strain_idx = int(gene_id.split(':')[0])
    return strain_idx in removed_indices


def _extract_gene_from_node_name(node_name: str) -> str:
    """
    Extract the gene_id from an orth_tree node name.
    Root nodes: plain gene_id (e.g. '2:0:3611')
    Internal nodes: 'gene_id_identity' (e.g. '2:0:3611_0.98')
    """
    if '_' in node_name:
        return node_name.rsplit('_', 1)[0]
    return node_name


def prune_orth_tree(orth_tree: nx.DiGraph,
                    removed_indices: set) -> tuple:
    """
    Remove genes of deleted strains from the orth_tree.

    For each node:
      - Remove members belonging to deleted strains
      - Update strains set
      - Recompute has_para
      - If a node becomes empty, mark it for removal
      - If the node name's gene belongs to a removed strain but the node
        still has surviving members, relabel it to a surviving member

    Returns (pruned_tree, relabel_map) where relabel_map maps
    old_node_name → new_node_name for relabeled nodes.
    """
    pruned = orth_tree.copy()
    nodes_to_remove = []
    relabel_map = {}

    for node in list(pruned.nodes()):
        attrs = pruned.nodes[node]
        if 'members' not in attrs:
            nodes_to_remove.append(node)
            continue

        # Filter members
        new_members = {m for m in attrs['members']
                       if not _gene_belongs_to_removed(m, removed_indices)}

        if not new_members:
            nodes_to_remove.append(node)
            continue

        # Update attributes
        new_strains = {int(m.split(':')[0]) for m in new_members}
        attrs['members'] = new_members
        attrs['strains'] = new_strains
        attrs['has_para'] = len(new_members) != len(new_strains)

        # Check if the node name itself refers to a gene from a removed strain
        gene_in_name = _extract_gene_from_node_name(node)
        if _gene_belongs_to_removed(gene_in_name, removed_indices):
            # Pick a surviving member as the new representative
            new_repre = next(iter(new_members))
            if '_' in node:
                # Internal node: preserve the identity suffix
                suffix = node.rsplit('_', 1)[1]
                new_name = f'{new_repre}_{suffix}'
            else:
                # Root node: plain gene_id
                new_name = new_repre
            relabel_map[node] = new_name

    # Remove empty nodes and fix edges
    for node in nodes_to_remove:
        # If this node has a parent and children, reconnect them
        parents = list(pruned.predecessors(node))
        children = list(pruned.successors(node))
        for parent in parents:
            for child in children:
                if child not in nodes_to_remove:
                    pruned.add_edge(parent, child)
        pruned.remove_node(node)

    # Relabel nodes whose representative gene was from a removed strain
    if relabel_map:
        pruned = nx.relabel_nodes(pruned, relabel_map, copy=False)
        logger.info(
            f'Relabeled {len(relabel_map)} nodes with new representatives')

    logger.info(f'Pruned orth_tree: removed {len(nodes_to_remove)} empty nodes, '
                f'{pruned.number_of_nodes()} nodes remaining')
    return pruned, relabel_map


def prune_raw_distance_graph(raw_dg: nx.Graph,
                             surviving_roots: set,
                             relabel_map: dict) -> nx.Graph:
    """
    Remove distance graph nodes that no longer correspond to roots in
    the pruned orth_tree, and relabel nodes whose representative gene
    was from a removed strain.

    The raw_distance_graph nodes are the final clustering representatives.
    These must match the orth_tree roots.  If a root was removed because
    all its members belonged to deleted strains, its corresponding
    distance graph node must also go.
    """
    pruned = raw_dg.copy()

    # First relabel nodes that were relabeled in the orth_tree.
    # Only relabel nodes that actually exist in the graph (root nodes only).
    dg_relabel = {old: new for old, new in relabel_map.items()
                  if old in pruned}
    if dg_relabel:
        pruned = nx.relabel_nodes(pruned, dg_relabel, copy=False)
        logger.info(f'Relabeled {len(dg_relabel)} distance graph nodes')

    # Then remove nodes not in surviving roots
    nodes_to_remove = [n for n in pruned.nodes() if n not in surviving_roots]
    pruned.remove_nodes_from(nodes_to_remove)
    logger.info(f'Pruned raw_distance_graph: removed {len(nodes_to_remove)} '
                f'nodes, {pruned.number_of_nodes()} remaining, '
                f'{pruned.number_of_edges()} edges')
    return pruned


def filter_alignment_file(src_path: str, dst_path: str,
                          surviving_nodes: set,
                          relabel_map: dict):
    """
    Filter alignment result file to keep only lines where both qseqid and
    sseqid are surviving representative nodes.

    If a node was relabeled (because its representative gene belonged to a
    removed strain), we need to match using the OLD name (as it appears in
    the alignment file) and write the NEW name.
    """
    # Only root-level relabelings matter for alignment (no '_' suffix).
    reverse = {old: new for old, new in relabel_map.items()
               if '_' not in old}
    # surviving_nodes has new names; we also need to accept old names
    old_names_surviving = set(reverse.keys())

    kept = 0
    total = 0
    with open(src_path, 'r') as fin, open(dst_path, 'w') as fout:
        for line in fin:
            total += 1
            parts = line.split('\t')
            if len(parts) >= 2:
                qid = parts[0]
                sid = parts[1]
                q_ok = qid in surviving_nodes or qid in old_names_surviving
                s_ok = sid in surviving_nodes or sid in old_names_surviving
                if q_ok and s_ok:
                    # Rewrite with new names if relabeled
                    if qid in reverse:
                        parts[0] = reverse[qid]
                    if sid in reverse:
                        parts[1] = reverse[sid]
                    fout.write('\t'.join(parts))
                    kept += 1
    logger.info(f'Filtered alignment: {kept}/{total} lines kept')


def build_pruned_pangenome(pg: Pangenome, removed_indices: set,
                           outdir: str, threads: int, gcode: int,
                           disable: bool) -> Pangenome:
    """
    Create a new Pangenome without the removed strains.
    Surviving strains keep their original indices (no reindexing).
    """
    new_pg = Pangenome(
        outdir=outdir,
        threads=threads,
        gcode=gcode,
        disable=disable
    )

    removed_gene_count = 0
    for idx, strain in pg.strain_dict.items():
        if idx in removed_indices:
            removed_gene_count += sum(strain.gene_num)
            continue
        new_pg.load_strain(strain)

    new_pg.total_gene_num = pg.total_gene_num - removed_gene_count

    logger.info(f'Pruned pangenome: {new_pg.strain_num} strains, '
                f'{new_pg.total_gene_num} genes '
                f'(removed {len(removed_indices)} strains, '
                f'{removed_gene_count} genes)')
    return new_pg


def filter_fasta_file(src_path: str, dst_path: str,
                      removed_indices: set):
    """
    Filter a FASTA file, removing sequences belonging to deleted strains.
    """
    kept = 0
    removed = 0
    with open(src_path, 'r') as fin, open(dst_path, 'w') as fout:
        skip = False
        for line in fin:
            if line.startswith('>'):
                gene_id = line[1:].strip()
                if _gene_belongs_to_removed(gene_id, removed_indices):
                    skip = True
                    removed += 1
                else:
                    skip = False
                    kept += 1
                    fout.write(line)
            else:
                if not skip:
                    fout.write(line)
    logger.info(f'Filtered FASTA: kept {kept}, removed {removed} sequences')


def filter_annot_file(src_path: str, dst_path: str,
                      removed_indices: set):
    """
    Filter an annotation TSV file, removing lines for deleted strains.
    The gene ID is in the first column (format: strain_idx:contig:gene).
    """
    kept = 0
    removed = 0
    with open(src_path, 'r') as fin, open(dst_path, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue
            parts = line.split('\t', 1)
            gene_id = parts[0]
            if _gene_belongs_to_removed(gene_id, removed_indices):
                removed += 1
            else:
                fout.write(line)
                kept += 1
    logger.info(f'Filtered annotation: kept {kept}, removed {removed} lines')


def _dump_csv_sparse(pg, outdir: str, surviving_indices: list,
                     prefix: str = 'pgap2.partition'):
    """
    Write output CSV/PAV/detail files, selecting only the columns
    corresponding to ``surviving_indices`` (sorted list of strain-dict
    keys that are still present).  This avoids empty columns for removed
    strains.
    """
    from pgap2.lib.pangenome import pan_judger

    logger.info(f'Dump csv matrix to {outdir}/{prefix}.gene_content.csv')

    strain_name_list = [pg.strain_dict[i].strain_name
                        for i in surviving_indices]
    header = '#Clust,' + ','.join(strain_name_list)
    header2 = '#Clust\t' + '\t'.join(strain_name_list)
    header3 = ('#Clust\tgene_name\tproduct\tgroup\trepre_gene\t'
               'min\tmean\tvar\tuni\tinvolved_strain\tpara_strain\t'
               'involved_gene\tpara_gene\t' + ','.join(strain_name_list))

    statistic_dict = OrderedDict({
        'Strict_core': 0, 'Core': 0, 'Soft_core': 0,
        'Shell': 0, 'Cloud': 0, 'Total': 0})
    total_strain_num = len(surviving_indices)

    with (open(f'{outdir}/{prefix}.gene_content.csv', 'w') as fh,
          open(f'{outdir}/{prefix}.gene_content.pav', 'w') as fh2,
          open(f'{outdir}/{prefix}.gene_content.detail.tsv', 'w') as fh3,
          open(f'{outdir}/{prefix}.summary_statistics.txt', 'w') as fh4):

        fh.write(f'{header}\n')
        fh2.write(f'{header2}\n')
        fh3.write(f'{header3}\n')

        for i, full_pan in enumerate(pg.pan_array):
            full_pav = pg.pav_array[i]
            full_sym = pg.pan_array_symbol[i]

            # Extract only surviving columns
            one_pan = [full_pan[j] for j in surviving_indices]
            one_pav = [full_pav[j] for j in surviving_indices]
            one_sym = [full_sym[j] for j in surviving_indices]

            para_strain = 0
            para_gene_num = 0
            involved_strain = 0
            involved_gene = 0
            for eg in one_pav:
                if eg > 0:
                    involved_strain += 1
                    involved_gene += eg
                if eg > 1:
                    para_strain += 1
                    para_gene_num += eg - 1

            group = pan_judger(
                ortho_num=involved_strain, total_num=total_strain_num)
            statistic_dict[group] += 1
            statistic_dict['Total'] += 1

            mn = pg.pan_attr[i]['min']
            uni = pg.pan_attr[i]['uni']
            mean = pg.pan_attr[i]['mean']
            var = pg.pan_attr[i]['var']
            repre_node = pg.pan_attr[i]['repre_node']
            gene_name = pg.pan_attr[i]['gene_name'] or '[]'
            gene_product = pg.pan_attr[i]['gene_product'] or '[]'

            fh.write('clust_{},{}\n'.format(i, ','.join(one_pan)))
            fh2.write('clust_{}\t{}\n'.format(
                i, '\t'.join(str(x) for x in one_pav)))
            fh3.write(
                'clust_{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'
                .format(i, gene_name, gene_product, group, repre_node,
                        mn, mean, var, uni, involved_strain, para_strain,
                        involved_gene, para_gene_num,
                        ','.join(one_sym)))

        for k, v in statistic_dict.items():
            freq = round(int(v) * 100 / int(statistic_dict['Total']), 2)
            limit = ''
            if k == 'Strict_core':
                limit = '(strains = 100%)'
                k = 'Strict core genes'
            elif k == 'Core':
                limit = '(99% <= strains < 100%)'
                k = 'Core genes'
            elif k == 'Soft_core':
                limit = '(95% <= strains < 99%)'
                k = 'Soft core genes'
            elif k == 'Shell':
                limit = '(15% <= strains < 95%)'
                k = 'Shell genes'
            elif k == 'Cloud':
                limit = '(0% <= strains < 15%)'
                k = 'Cloud genes'
            elif k == 'Total':
                limit = '(0% <= strains <= 100%)'
                k = 'Total genes'
            else:
                logger.error('ERROR: tell me at github')
            fh4.write('{}\t{}\t{}\t{}%\n'.format(k, limit, v, freq))


def build_guide_dict_minus(G: nx.Graph, removed_indices: set) -> dict:
    """
    Build a guide dictionary from the previous graph, excluding genes
    from removed strains.
    """
    guide_dict = {}
    for node in G.nodes():
        members = G.nodes[node].get('members', set())
        for member in members:
            if not _gene_belongs_to_removed(member, removed_indices):
                guide_dict[member] = node
    logger.info(f'Built guide dict: {len(guide_dict)} surviving genes')
    return guide_dict


# ─────────────────────── Main ─────────────────────────────────────────────────

def main(prev_dir: str, removal_file: str, outdir: str,
         threads: int, disable: bool):

    # ── Load parameters from previous result ──
    logger.info(f'Loading parameters from {prev_dir}/basic.pkl')
    with open(f'{prev_dir}/basic.pkl', 'rb') as fh:
        pkl: PklCheck = pickle.load(fh)
        basic = pkl.data_dump('basic')
        params = basic.params

    # Inherit all parameters
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
    retrieve = params.get('retrieve', False)
    aligner = params.get('aligner', 'diamond')
    gcode = params.get('gcode', 11)
    falen = params.get('falen', 20)
    clust_method = params.get('clust_method', 'mmseqs2')
    id_attr_key = params.get('id_attr_key', 'ID')
    type_filter = params.get('type_filter', 'CDS')
    annot = params.get('annot', False)

    logger.info('Inherited parameters from previous result:')
    for k, v in params.items():
        logger.debug(f'  {k}: {v}')

    # ── Load pangenome and tree ──
    logger.info(f'Loading pangenome and tree from {prev_dir}/preprocess.pkl')
    with open(f'{prev_dir}/preprocess.pkl', 'rb') as fh:
        pkl_pre: PklCheck = pickle.load(fh)
        pg_prev = pkl_pre.data_dump('pangenome')
        tree_prev = pkl_pre.data_dump('tree')

    # ── Load removal list and resolve indices ──
    names_to_remove = load_removal_list(removal_file)
    removed_indices = resolve_strain_indices(pg_prev, names_to_remove)

    removed_names = sorted(names_to_remove)
    logger.info(f'Strains to remove: {removed_names}')

    # ── Create output directory ──
    os.makedirs(outdir, exist_ok=True)

    # ── Prune orth_tree ──
    logger.info('Pruning orth_tree...')
    orth_tree_prev = tree_prev.orth_identity_tree
    pruned_orth_tree, relabel_map = prune_orth_tree(
        orth_tree_prev, removed_indices)

    # Get the set of surviving roots (these match raw_distance_graph nodes)
    surviving_roots = {n for n in pruned_orth_tree.nodes()
                       if pruned_orth_tree.in_degree(n) == 0}
    logger.info(f'Surviving orth_tree roots: {len(surviving_roots)}')

    # ── Prune raw_distance_graph ──
    logger.info('Pruning raw distance graph...')
    raw_dg_prev = tree_prev.raw_distance_graph
    pruned_raw_dg = prune_raw_distance_graph(raw_dg_prev, surviving_roots,
                                             relabel_map)

    # ── Filter alignment file ──
    logger.info('Filtering alignment result...')
    prev_alignment = tree_prev.alignment_result
    new_alignment = f'{outdir}/{aligner}.filtered.result'
    filter_alignment_file(prev_alignment, new_alignment, surviving_roots,
                          relabel_map)

    # ── Build pruned pangenome ──
    logger.info('Building pruned pangenome...')
    pg = build_pruned_pangenome(pg_prev, removed_indices, outdir,
                                threads, gcode, disable)

    # ── Filter protein and annotation files ──
    logger.info('Filtering protein and annotation files...')
    prev_prot = f'{prev_dir}/total.involved_prot.fa'
    prev_annot = f'{prev_dir}/total.involved_annot.tsv'
    new_prot = f'{outdir}/total.involved_prot.fa'
    new_annot = f'{outdir}/total.involved_annot.tsv'

    filter_fasta_file(prev_prot, new_prot, removed_indices)
    filter_annot_file(prev_annot, new_annot, removed_indices)

    pg.load_annot_file(new_annot)
    pg.load_prot_file(new_prot)

    # ── Copy genome_index for surviving strains (if retrieve mode) ──
    if retrieve:
        import shutil
        genome_idx_prev = f'{prev_dir}/genome_index'
        genome_idx_out = f'{outdir}/genome_index'
        if os.path.exists(genome_idx_prev):
            os.makedirs(genome_idx_out, exist_ok=True)
            for contig_dir in os.listdir(genome_idx_prev):
                src_dir = os.path.join(genome_idx_prev, contig_dir)
                if not os.path.isdir(src_dir):
                    continue
                for entry in os.listdir(src_dir):
                    try:
                        if int(entry) in removed_indices:
                            continue
                    except ValueError:
                        continue
                    src = os.path.join(src_dir, entry)
                    dst = os.path.join(genome_idx_out, contig_dir, entry)
                    if not os.path.exists(dst):
                        os.makedirs(os.path.dirname(dst), exist_ok=True)
                        shutil.copytree(src, dst)

    # ── Build pruned Tree object ──
    logger.info('Building pruned Tree object...')
    tree = Tree()
    tree.load_alignment_result(new_alignment)
    tree.load_ortho_identity_tree(pruned_orth_tree)
    tree.load_distance_graph(pruned_raw_dg, raw=True)

    # ── Save preprocess checkpoint ──
    logger.info('Saving preprocess checkpoint...')
    file_dict = OrderedDict()
    for strain in pg.strain_dict.values():
        file_dict[strain.strain_name] = {}

    pickle_preprocess = PklCheck(outdir=outdir, name='preprocess')
    pickle_preprocess.load('file_dict', main_data=file_dict)
    pickle_preprocess.load('pangenome', main_data=pg, parameter={
        'orth_id': orth_id, 'para_id': para_id, 'dup_id': dup_id,
        'accurate': accurate, 'id_attr_key': id_attr_key,
        'type_filter': type_filter, 'coverage': coverage,
        'AS': AS, 'AL': AL, 'LD': LD, 'retrieve': retrieve,
        'evalue': evalue, 'aligner': aligner, 'clust_method': clust_method,
        'annot': annot, 'falen': falen})
    pickle_preprocess.load('tree', main_data=tree)
    pickle_preprocess.pickle_()

    # ── Build guide dict from previous graph ──
    logger.info('Loading previous graph for guide construction...')
    with open(f'{prev_dir}/graph.pkl', 'rb') as fh:
        G_prev: PklCheck = pickle.load(fh)
        G_prev = G_prev.data_dump('graph')

    guide_dict = build_guide_dict_minus(G_prev, removed_indices)

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
        pg=pg, tree=tree, guide_dict=guide_dict, new_strain_indices=set())

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
                                      new_strain_indices=set())

    logger.info('Double check through gene synteny')
    G = merge_by_synteny(G, pg, tree,
                         context_sim=context_similarity,
                         flank=flank,
                         sensitivity=sensitivity,
                         ins=ins,
                         guide_dict=guide_dict,
                         new_strain_indices=set())

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
                                          new_strain_indices=set())
        G = merge_by_synteny(G=G, pg=pg, tree=tree,
                             context_sim=context_similarity,
                             flank=flank,
                             sensitivity=sensitivity,
                             ins=ins, step=11,
                             guide_dict=guide_dict,
                             new_strain_indices=set())

    # ─────────────────────── Organize results ────────────────────────────────
    logger.info('Organize the results')
    # Because surviving strains keep original indices (sparse), we need
    # the pan arrays to be large enough for the max strain index, not just
    # the count of surviving strains.
    max_idx = max(pg.strain_dict.keys())
    n_slots = max_idx + 1
    pg.one_pan = [""] * n_slots
    pg.one_pan_symbol = [""] * n_slots
    pg.one_pav = [0] * n_slots
    pg.pan_array = []
    pg.pan_array_symbol = []
    pg.pav_array = []
    pg.pan_attr = []

    # Collect the sorted list of surviving strain indices for output
    surviving_indices = sorted(pg.strain_dict.keys())

    bar = tqdm(range(G.number_of_nodes()), unit=f" Organize",
               disable=disable, desc=tqdm_.step(-1))
    for node in G.nodes():
        bar.update()
        my_pan_clust = get_pan_clust(G, pg, tree, node)
        pg.load_one_pan(pan_clust=my_pan_clust)
    bar.close()

    logger.info('Dump the gene content matrix')
    _dump_csv_sparse(pg, outdir, surviving_indices,
                     prefix='pgap2.partition')

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
