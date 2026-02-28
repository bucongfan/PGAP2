import os
import pickle
import argparse
import networkx as nx

from tqdm import tqdm
from loguru import logger
from collections import defaultdict

from pgap2.lib.tree import Tree
from pgap2.lib.basic import Basic
from pgap2.lib.pklcheck import PklCheck
from pgap2.lib.pangenome import Pangenome

from pgap2.utils.supply import tqdm_
from pgap2.utils.generate_tree import generate_tree
from pgap2.utils.gene_retriever import retrieve_gene
from pgap2.utils.arrangement_detector import merge_by_synteny
from pgap2.utils.data_loader import file_parser, get_file_dict

from pgap2.utils.partition import (
    generate_network, get_pan_clust, merge_by_similarity,
    mcl, get_expect_identity
)

"""
Add-genome mode for PGAP2.

This module adds one or more genomes to an existing pangenome analysis by reusing
the same pipeline as partition (mcl → generate_network → merge_by_similarity →
merge_by_synteny), but with a guide dictionary built from the previous run's
results. The guide dictionary enables skipping complex merge judgments for
node pairs that contain only old genes.

Guide logic:
- If both nodes in a merge candidate contain ONLY old genes:
  → Use the guide dict to decide (were they in the same cluster before?
    yes → merge, no → skip)
- If either node contains new genes:
  → Proceed with normal merge judgment

"""


def merge_pangenome(previous_pg: Pangenome, new_pg: Pangenome) -> Pangenome:
    """
    Merge a new pangenome with the previous pangenome.

    This function combines the strain information from both pangenomes,
    ensuring that the new genome's strain indices continue from where 
    the previous pangenome ended.

    Args:
        previous_pg: The existing pangenome object with previous strains
        new_pg: The new pangenome object containing newly added strain(s)

    Returns:
        Pangenome: A merged pangenome object containing all strains
    """
    logger.info(
        f'Merging pangenomes: {previous_pg.strain_num} previous strains + {new_pg.strain_num} new strain(s)')

    # Create a new merged pangenome with the same configuration as previous
    merged_pg = Pangenome(
        outdir=new_pg.outdir,  # Use new output directory
        threads=new_pg.threads,
        gcode=new_pg.gcode,
        disable=new_pg.disable_tqdm
    )

    # Copy all strains from previous pangenome
    for strain_index, strain in previous_pg.strain_dict.items():
        merged_pg.load_strain(strain)

    # Add all strains from new pangenome (they should already have correct indices)
    for strain_index, strain in new_pg.strain_dict.items():
        merged_pg.load_strain(strain)
        logger.debug(
            f'Added new strain: {strain.strain_name} (index: {strain_index})')

    # Merge total gene numbers
    merged_pg.total_gene_num = previous_pg.total_gene_num + new_pg.total_gene_num

    # Copy other important attributes from previous pangenome
    merged_pg.orth_id = previous_pg.orth_id
    merged_pg.para_id = previous_pg.para_id
    merged_pg.dup_id = previous_pg.dup_id
    merged_pg.accurate = previous_pg.accurate
    merged_pg.exhaust_orth = previous_pg.exhaust_orth
    merged_pg.retrieve = previous_pg.retrieve
    merged_pg.evalue = previous_pg.evalue
    merged_pg.aligner = previous_pg.aligner
    merged_pg.LD = previous_pg.LD
    merged_pg.AL = previous_pg.AL
    merged_pg.AS = previous_pg.AS

    logger.info(
        f'Merged pangenome: {merged_pg.strain_num} total strains, {merged_pg.total_gene_num} total genes')

    return merged_pg


def build_guide_dict(previous_G: nx.Graph) -> dict:
    """
    Build a guide dictionary from the previous run's final graph.

    Maps each gene member to its node name (cluster ID) in the previous graph.
    This allows merge functions to quickly check whether two genes were in the
    same cluster in the previous run.

    Args:
        previous_G: The final network graph from the previous partition run.

    Returns:
        dict: {gene_id: cluster_node_name} for every gene in the previous graph.
    """
    guide_dict = {}
    for node in previous_G.nodes():
        members = previous_G.nodes[node].get('members', set())
        for member in members:
            guide_dict[member] = node
    logger.info(
        f'Built guide dict: {len(guide_dict)} genes → {previous_G.number_of_nodes()} clusters')
    return guide_dict


def merge_files(previous_file: str, current_file: str, output_file: str, skip_header: bool = False) -> str:
    """
    Merge the previous file with the current file.

    Args:
        previous_file: Path to the previous file
        current_file: Path to the current (new) file
        output_file: Path to the output merged file
        skip_header: If True, skip the first line (header) of the current file when appending

    Returns:
        str: Path to the merged file
    """
    logger.info(f'Merging files:')
    logger.info(f'  Previous: {previous_file}')
    logger.info(f'  Current:  {current_file}')

    with open(output_file, 'w') as out_fh:
        # Write previous file content first
        with open(previous_file, 'r') as fh:
            out_fh.write(fh.read())

        # Append current file content
        with open(current_file, 'r') as fh:
            if skip_header:
                # Skip header line
                next(fh, None)
            out_fh.write(fh.read())

    logger.info(f'  Merged:   {output_file}')

    return output_file


def main(indir: str, outdir: str, previous_dir: str, aligner: str, clust_method: str, falen: int, threads: int, id_attr_key: str, type_filter: str, annot: bool, gcode: int, retrieve: bool, disable: bool):

    # Load previous parameters
    logger.info(f'Loading parameters from {previous_dir}/basic.pkl')
    with open(f'{previous_dir}/basic.pkl', 'rb') as fh:
        previous: PklCheck = pickle.load(fh)
        basic = previous.data_dump('basic')
        params = basic.params

    # Log current arguments
    logger.info('Current arguments:')
    current_args = {
        'aligner': aligner,
        'clust_method': clust_method,
        'falen': falen,
        'threads': threads,
        'id_attr_key': id_attr_key,
        'type_filter': type_filter,
        'annot': annot,
        'gcode': gcode,
        'retrieve': retrieve,
        'disable': disable
    }
    for k, v in current_args.items():
        logger.info(f'{k}: {v}')

    # Check for consistency
    check_keys = ['aligner', 'clust_method', 'falen',
                  'id_attr_key', 'type_filter', 'annot', 'gcode', 'retrieve']
    for key in check_keys:
        if key in params and params[key] != current_args[key]:
            logger.warning(
                f"Parameter mismatch: {key} (Previous: {params[key]}, Current: {current_args[key]})")

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

    logger.info('Inherited parameters:')
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
        logger.info(f'{k}: {v}')

    decode_status = False
    file_dict = get_file_dict(indir)

    # Validate that at least one genome is being added
    if len(file_dict) < 1:
        logger.error(f'No genome file found in {indir}')
        raise ValueError(
            f'Invalid input: no genome found in {indir}')

    logger.info(
        f'Found {len(file_dict)} genome(s) in {indir}')

    # Check for duplicate strain names between new and previous genomes
    existing_strain_names = {
        strain.strain_name for strain in basic.strain_dict.values()}
    new_strain_names = set(file_dict.keys())
    duplicates = existing_strain_names & new_strain_names
    if duplicates:
        logger.warning(
            f'Found {len(duplicates)} duplicate strain(s) already in previous results, skipping:')
        for dup in sorted(duplicates):
            logger.warning(f'  - {dup}')
        # Remove duplicates from file_dict before index assignment
        for dup in duplicates:
            del file_dict[dup]

    if len(file_dict) < 1:
        logger.error(
            f'No new genome to add after removing {len(duplicates)} duplicate(s).')
        raise ValueError(
            f'All input genomes already exist in the previous results.')

    num_new_genomes = len(file_dict)
    logger.info(
        f'{num_new_genomes} new genome(s) to be added from {indir}')
    for strain_name in file_dict.keys():
        logger.info(f'  - {strain_name}')

    if os.path.exists(f'{outdir}/preprocess.pkl'):
        '''
        Found a previous preprocess.pkl file, loading...
        This file contains the parameters and file structure of the previous run.
        If the parameters are not match, it will raise a ValueError.
        If the file structure has changed, it will warn the user and reload the file structure from the current input directory.
        If the file structure is the same, it will load the pangenome and tree from the previous run.
        If the pangenome has invalid genes, it will warn the user and continue to the next step.
        '''
        logger.info(f'Found {outdir}/preprocess.pkl')
        logger.info(f'Loading...')
        with open(f'{outdir}/preprocess.pkl', 'rb') as fh:
            previous: PklCheck = pickle.load(fh)
            logger.info(f'Check the previous file parameters...')
            decode_status = previous.decode(
                orth_id=orth_id, para_id=para_id, dup_id=dup_id, accurate=accurate, coverage=coverage, id_attr_key=id_attr_key, type_filter=type_filter, LD=LD, AS=AS, AL=AL, evalue=evalue, aligner=aligner, clust_method=clust_method, falen=falen, annot=annot, retrieve=retrieve,)
            if decode_status:
                # success
                pg = previous.data_dump('pangenome')
                tree = previous.data_dump('tree')
                previous_file_dict = previous.data_dump('file_dict')

                if previous_file_dict != file_dict:
                    logger.warning(
                        f'File structure has changed')
                    total_name = list(previous_file_dict.keys(
                    )) + [k for k in file_dict.keys() if k not in previous_file_dict]

                    max_width = max(
                        [len(name) for name in total_name+['Previous', 'Current']])+2
                    logger.warning(
                        f'{"Previous":<{max_width}}\t{"Current":<{max_width}}')
                    logger.warning(f'{"-"*(max_width*2)}')
                    new_add = 0
                    loaded_count = 0
                    for strain in total_name:
                        cur_name = strain if strain in file_dict else ''
                        pre_name = strain if strain in previous_file_dict else ''
                        if cur_name != pre_name:  # empty
                            logger.warning(
                                f'{pre_name:<{max_width}}\t{cur_name:<{max_width}}')
                        if cur_name and not pre_name:
                            new_add += 1
                        else:
                            loaded_count += 1
                    logger.warning(f'{"-"*(max_width*2)}')
                    len_prev = len(previous_file_dict)
                    len_cur = len(file_dict)
                    logger.warning(
                        f'{len_prev:<{max_width}}\t{len_cur:<{max_width}}')

                    if new_add:
                        logger.warning(
                            f'Total {new_add} new strain added. Make sure the preprocess.pkl I loaded is the right one!!!')
                        logger.warning(
                            f'I will reload the file structure from the current input: {indir}')
                        decode_status = False
                    if loaded_count < 2:
                        logger.error(
                            f'Loaded file has less than 2 strains, it is not a valid file that may cause the --exclude_outlier parameter in preprocess step filtered much strains')
                        logger.error(
                            'Please check the input file quality and rerun the preprocess step or just begin from the partition step')
                        raise ValueError('Invalid preprocess.pkl file')

                    if decode_status:
                        file_dict = previous_file_dict
                else:
                    logger.info(
                        f'Load previous file structure from {outdir}/pgap2.pkl')

                if decode_status:
                    total_bad_gene_num = 0
                    for strain in tqdm(pg.strain_dict, unit=' strain', disable=disable, desc=tqdm_.step(1)):
                        bad_num = pg.strain_dict[strain].bad_gene_num
                        if bad_num > 0:
                            total_bad_gene_num += bad_num
                            logger.warning(
                                f'{strain} invalid gene count: {bad_num}')
                    if total_bad_gene_num > 0:
                        logger.info(
                            f'Total invalid gene count: {total_bad_gene_num}. Check it in log file: {outdir}/preprocess.log')

                    for _ in tqdm([dup_id, orth_id], unit=f" clust iteration", disable=disable, desc=tqdm_.step(2)):
                        ...

                    # Get the new strain index (it's the last one added)
                    # Load previous_dir pangenome to get the original strain count
                    with open(f'{previous_dir}/preprocess.pkl', 'rb') as prev_fh:
                        prev_pkl: PklCheck = pickle.load(prev_fh)
                        prev_pg = prev_pkl.data_dump('pangenome')
                        start_strain_index = prev_pg.strain_num
                        logger.info(f'New strain index: {start_strain_index}')
            else:
                logger.warning(
                    f'Previous file parameters is not match, start partition from the begining')

    if decode_status is False:
        '''
        Load strain from input directory
        If the previous preprocess.pkl file is not found or the parameters are not match,
        it will load the strain from the input directory and create a new pangenome object.
        '''
        # load new genome(s)
        logger.info(f'Load {num_new_genomes} strain(s) from {indir}')

        # Load previous pangenome to get the starting strain index
        with open(f'{previous_dir}/preprocess.pkl', 'rb') as fh:
            previous_pkl: PklCheck = pickle.load(fh)
            previous_pg = previous_pkl.data_dump('pangenome')
            previous_tree = previous_pkl.data_dump('tree')
            start_strain_index = previous_pg.strain_num
            logger.info(f'Starting strain index: {start_strain_index}')

        # Parse new genome(s) with correct starting index
        new_pg = file_parser(
            indir=indir, outdir=outdir, annot=annot, threads=threads, disable=disable,
            retrieve=retrieve, falen=falen, gcode=gcode, id_attr_key=id_attr_key,
            type_filter=type_filter, prefix='partition', start_index=start_strain_index, run_type='add')
        # Merge previous pangenome with new pangenome
        pg = merge_pangenome(previous_pg, new_pg)

        file_prot = f'{outdir}/add.involved_prot.fa'
        file_annot = f'{outdir}/add.involved_annot.tsv'

        # Merge protein and annotation files
        merged_prot = merge_files(
            previous_file=f'{previous_dir}/total.involved_prot.fa',
            current_file=file_prot,
            output_file=f'{outdir}/total.involved_prot.fa',
            skip_header=False
        )
        merged_annot = merge_files(
            previous_file=f'{previous_dir}/total.involved_annot.tsv',
            current_file=file_annot,
            output_file=f'{outdir}/total.involved_annot.tsv',
            skip_header=True  # Skip header in new annot file
        )

        pg.load_annot_file(merged_annot)
        pg.load_prot_file(merged_prot)
        logger.info(
            f'Create distane tree with {pg.strain_num} strains')
        logger.info(
            f'Clustering with orth_id: {orth_id}, para_id: {para_id}, dup_id: {dup_id}')
        tree = generate_tree(
            input_file=merged_prot, orth_list=[dup_id, orth_id], outdir=pg.outdir, evalue=evalue, aligner=aligner, falen=falen, disable=disable, threads=threads, max_targets=max_targets, coverage=coverage, ID=para_id, LD=LD, AS=AS, AL=AL, clust_method=clust_method)

        logger.info(f'Pangenome and tree loaded successfully.')
        logger.info(
            f'To save the complete information of this project for breakpoint resume...')
        pickle_preprocess = PklCheck(outdir=outdir, name='preprocess')
        pickle_preprocess.load('file_dict', main_data=file_dict)
        pickle_preprocess.load('pangenome', main_data=pg, parameter={'orth_id': orth_id, 'para_id': para_id, 'dup_id': dup_id, 'accurate': accurate,
                                                                     'id_attr_key': id_attr_key, 'type_filter': type_filter,
                                                                     'coverage': coverage, 'AS': AS, 'AL': AL, 'LD': LD, 'retrieve': retrieve,
                                                                     'evalue': evalue, 'aligner': aligner, 'clust_method': clust_method,
                                                                     'annot': annot, 'falen': falen})
        pickle_preprocess.load('tree', main_data=tree)
        pickle_preprocess.pickle_()

    with open(f'{previous_dir}/graph.pkl', 'rb') as fh:
        previous_G: PklCheck = pickle.load(fh)
        previous_G = previous_G.data_dump('graph')

    # Build guide dict from previous results
    guide_dict = build_guide_dict(previous_G)
    new_strain_indices = set(
        range(start_strain_index, start_strain_index + num_new_genomes))
    logger.info(f'New strain indices: {new_strain_indices}')

    '''
    load necessary parameters and file paths used to quick downstream analysis
    '''
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

    with open(f'{previous_dir}/preprocess.pkl', 'rb') as fh:
        previous_pkl: PklCheck = pickle.load(fh)
        previous_tree = previous_pkl.data_dump('tree')
    # -----------------------------------partition step-----------------------------------#
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
    logger.info(f'Clean up the distance graph according to paralogous genes')
    tree.update_distance_graph(disable=disable)
    logger.info(f'Merge by gene similarity')

    G, pg, tree = merge_by_similarity(G=G, pg=pg, tree=tree,
                                      sensitivity=sensitivity,
                                      radius=radius, fast=fast_mode,
                                      context_sim=context_similarity,
                                      flank=flank,
                                      disable=disable,
                                      guide_dict=guide_dict,
                                      new_strain_indices=new_strain_indices,)

    logger.info(f'Double check through gene synteny')
    G = merge_by_synteny(G, pg, tree,
                         context_sim=context_similarity,
                         flank=flank,
                         sensitivity=sensitivity,
                         ins=ins,
                         guide_dict=guide_dict,
                         new_strain_indices=new_strain_indices,
                         )

    logger.info(f'Reload the gene annotation')
    pg.reload_annot_file(retrieve=retrieve)

    if retrieve:
        logger.info(f'Retrieve gene from missing')
        G, pg, tree = retrieve_gene(G, pg, tree)
        G, pg, tree = merge_by_similarity(G=G, pg=pg, tree=tree,
                                          sensitivity=sensitivity,
                                          radius=radius, fast=fast_mode,
                                          context_sim=context_similarity,
                                          flank=flank,
                                          disable=disable, step=10,
                                          guide_dict=guide_dict,
                                          new_strain_indices=new_strain_indices,)
        G = merge_by_synteny(G=G, pg=pg, tree=tree,
                             context_sim=context_similarity,
                             flank=flank,
                             sensitivity=sensitivity,
                             ins=ins, step=11,
                             guide_dict=guide_dict,
                             new_strain_indices=new_strain_indices,
                             )

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

    logger.info(
        f'To save the basic results of this project for downstream visulization...')
    pickle_basic = PklCheck(outdir=outdir, name='basic')
    pickle_basic.load('basic', main_data=Basic(pg=pg))
    pickle_basic.pickle_()
    return 0
