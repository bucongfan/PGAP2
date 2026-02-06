import os
import shutil
import networkx as nx

from loguru import logger
from tqdm import tqdm
from functools import partial
from concurrent.futures import ThreadPoolExecutor

from pgap2.lib.tree import Tree
from pgap2.utils.supply import run_command
from pgap2.utils.supply import sfw, tqdm_

"""
clustering using cdhit or mmseqs2, and then build a tree based on the clustering result.

input:
- input_file: Input file containing gene sequences.
- orth_list: List of orthologous gene identities for clustering.
- outdir: Output directory for storing results.

params:
- coverage: Coverage threshold for clustering.
- evalue: E-value threshold for filtering alignments.
- falen: Flank length for context comparison.
- disable: Boolean to disable progress bar.
- threads: Number of threads to use for parallel processing.
- ID: Identity threshold for filtering alignments.
- LD: Length difference threshold for filtering alignments.
- AL: Alignment coverage threshold for filtering alignments.
- AS: Alignment score threshold for filtering alignments.
- aligner: Aligner to use for sequence alignment ('diamond' or 'blastp').
- clust_method: Clustering method to use ('cdhit' or 'mmseqs2').

output:
- Tree object containing the orthologous gene relationships and distances.
"""


def clust_recorder(subject: dict, query: dict, tag: str) -> dict:
    '''
    subject and query:
    dict={clust:[fa1,fa2,fa3]}
    subject: tmpdict is the previous iterative dict
    query: mydict is the current iterative dict
    '''
    previous_group = len(subject)
    if not subject:  # 初次进入
        changed_group = len(query)
        previous_group = sum([len(query[_]) for _ in query])
        logger.debug(
            f'clust_recorder: In {tag}, group were updated from {previous_group} to {changed_group}')
        return query
    else:
        for group in query:
            for header in query[group]:
                if group == header:
                    continue
                if group in subject and header in subject:
                    subject[group].extend(subject[header])
                    del subject[header]
                else:
                    raise Exception(
                        f'{group} or {header} not in previous clusters')
    changed_group = len(subject)
    logger.debug(
        f'clust_recorder: In {tag}, group were updated from {previous_group} to {changed_group}')
    return subject


def process_lines(lines, ID: int = 0, LD: int = 0.7, AL: float = 0, AS: float = 0):
    edges = []
    result = []
    for line in lines:
        lines = line.rstrip().split('\t')
        qseqid = lines[0]
        sseqid = lines[1]
        pident = lines[2]
        hsp = int(lines[3])
        qlen = int(lines[12])
        slen = int(lines[13])
        pident = round(float(pident) / 100, 3)
        len_diff = min(int(qlen), int(slen)) / max(int(qlen), int(slen))
        if qlen >= slen:
            al_cov = hsp / qlen
            as_cov = hsp / slen
        else:
            al_cov = hsp / slen
            as_cov = hsp / qlen

        if pident < ID or len_diff < LD or al_cov < AL or as_cov < AS:
            continue
        result.append(line)
        if qseqid != sseqid:
            edges.append((qseqid, sseqid, pident))
    return result, edges


def load_cdhit_result(output_file_clstr: str) -> dict:
    mydict = {}
    with open(output_file_clstr) as fh:
        cluster = []
        repre_node = None
        for line in fh:
            line = line.strip()
            if line[0] == ">":
                if cluster and repre_node:
                    mydict.update({repre_node: cluster})
                cluster = []
                repre_node = None
            else:
                seq_name = line.split(">")[1].split("...")[0]
                cluster.append(seq_name)
                repre_node = seq_name if line.endswith('*') else repre_node
        if cluster and repre_node:
            mydict.update({repre_node: cluster})
    return mydict


def run_cdhit(input_file: str,
              outdir: str,
              id: float,  # identity -c
              l: int = 10,  # length of throw_away_sequences, default 10
              s: float = 0.0,  # length difference cutoff, default 0.0
              b: int = 20,  # band_width of alignment, default 20
              threads: int = 1):
    s = s if id > 0.95 else 0
    output_file = f'{outdir}/repre_node'
    output_file_clstr = f'{output_file}.clstr'
    # if not os.path.exists(output_file_clstr):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    # recommanded word_size by cdhit manual
    if 0.7 < id <= 1:
        word_size = 5
    else:
        word_size = 4
    run_command(
        '{} -c {} -i {} -o {} -T {} -s {} -M 0 -d 256 -n {} -b {} -g 1 -l {}'.format(sfw.cdhit, id, input_file, output_file, threads, s, word_size, b, l))

    mydict = load_cdhit_result(output_file_clstr)
    return mydict, output_file


def run_mmseq2(data, data_type: str, id: float, coverage: float, outdir: str, threads: int = 8):
    '''
    Deprecated method
    '''
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if data_type == "fasta":
        fa = data
        run_command(
            f'mmseqs createdb --shuffle 1 {fa} {outdir}/seq.db -v 0')
        data_index = f'{outdir}/seq.db'
    elif data_type == "index":
        data_index = data
    else:
        raise (f'ERROR. Tell me on github')

    run_command(
        f'{sfw.mmseqs2} linclust {data_index} {outdir}/seq.clst {outdir}/tmp --min-seq-id {id} -c {coverage}  --threads {threads} -v 0')
    run_command(
        f'{sfw.mmseqs2} createtsv --first-seq-as-repr 1 {data_index} {data_index} {outdir}/seq.clst {outdir}/this_clust.tab --threads {threads}')
    run_command(
        f'{sfw.mmseqs2} createsubdb {outdir}/seq.clst {data_index} {outdir}/seq.clst.rep')

    mydict = {}
    with open(f'{outdir}/this_clust.tab') as fh:
        for line in fh:
            group, header = line.strip().split('\t')
            if group not in mydict:
                mydict[group] = []
            mydict[group].append(header)
    return mydict, f'{outdir}/seq.clst.rep'


def run_alignment(input_file: str, outdir: str, threads: int = 1, max_targets: int = 2000, evalue=1E-6, aligner: str = 'diamond') -> str:
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if aligner == 'diamond':
        # make diamond database
        diamond_db = os.path.join(outdir, 'diamond_db')
        run_command(
            f'{sfw.diamond} makedb --in {input_file} -d {diamond_db} -p {threads} --quiet')

        # run diamond blastp
        diamond_result = os.path.join(outdir, 'diamond.tsv')
        run_command(f'{sfw.diamond} blastp -q {input_file} -d {diamond_db} -p {threads} -e {evalue} -k {max_targets} -o {diamond_result} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen  &>/dev/null')
        return diamond_result
    elif aligner == 'blastp':
        # make blast database
        blast_db = os.path.join(outdir, 'blast_db')
        run_command(
            f'{sfw.makeblastdb} -in {input_file} -dbtype prot -out {blast_db} -parse_seqids')

        # run blastp
        blast_result = os.path.join(outdir, 'blast.tsv')
        run_command(
            f'{sfw.blastp} -query {input_file} -db {blast_db} -num_threads {threads} -out {blast_result} -evalue {evalue} -max_target_seqs {max_targets} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"')
        return blast_result


def run_incremental_alignment(new_file: str, prev_file: str, outdir: str, threads: int = 1, max_targets: int = 2000, evalue=1E-6, aligner: str = 'diamond') -> str:
    """
    Run incremental alignment for new sequences:
    1. Align new sequences against themselves (new vs new)
    2. Align new sequences against previous representatives (new vs prev)

    This avoids re-aligning all previous sequences, significantly improving performance.

    Args:
        new_file: FASTA file with new representative sequences
        prev_file: FASTA file with previous representative sequences
        outdir: Output directory
        threads: Number of threads
        max_targets: Maximum number of alignment targets
        evalue: E-value threshold
        aligner: Aligner to use ('diamond' or 'blastp')

    Returns:
        Path to merged alignment result file
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    logger.info(f'- Running incremental alignment: new vs new + new vs previous')

    if aligner == 'diamond':
        # Step 1: Align new sequences against themselves (new vs new)
        logger.debug(f'  Step 1/2: Aligning new sequences against themselves')
        new_db = os.path.join(outdir, 'new_diamond_db')
        run_command(
            f'{sfw.diamond} makedb --in {new_file} -d {new_db} -p {threads} --quiet')

        new_vs_new_result = os.path.join(outdir, 'diamond_new_vs_new.tsv')
        run_command(
            f'{sfw.diamond} blastp -q {new_file} -d {new_db} -p {threads} -e {evalue} '
            f'-k {max_targets} -o {new_vs_new_result} '
            f'--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen &>/dev/null')

        # Step 2: Align new sequences against previous representatives (new vs prev)
        logger.debug(
            f'  Step 2/2: Aligning new sequences against previous representatives')
        prev_db = os.path.join(outdir, 'prev_diamond_db')
        run_command(
            f'{sfw.diamond} makedb --in {prev_file} -d {prev_db} -p {threads} --quiet')

        new_vs_prev_result = os.path.join(outdir, 'diamond_new_vs_prev.tsv')
        run_command(
            f'{sfw.diamond} blastp -q {new_file} -d {prev_db} -p {threads} -e {evalue} '
            f'-k {max_targets} -o {new_vs_prev_result} '
            f'--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen &>/dev/null')

        # Merge the two results
        merged_result = os.path.join(outdir, 'diamond.tsv')
        with open(merged_result, 'w') as out_fh:
            # Read new vs new
            with open(new_vs_new_result, 'r') as fh:
                out_fh.write(fh.read())
            # Read new vs prev
            with open(new_vs_prev_result, 'r') as fh:
                out_fh.write(fh.read())

        return merged_result

    elif aligner == 'blastp':
        # Step 1: Align new sequences against themselves (new vs new)
        logger.debug(f'  Step 1/2: Aligning new sequences against themselves')
        new_db = os.path.join(outdir, 'new_blast_db')
        run_command(
            f'{sfw.makeblastdb} -in {new_file} -dbtype prot -out {new_db} -parse_seqids')

        new_vs_new_result = os.path.join(outdir, 'blast_new_vs_new.tsv')
        run_command(
            f'{sfw.blastp} -query {new_file} -db {new_db} -num_threads {threads} '
            f'-out {new_vs_new_result} -evalue {evalue} -max_target_seqs {max_targets} '
            f'-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"')

        # Step 2: Align new sequences against previous representatives (new vs prev)
        logger.debug(
            f'  Step 2/2: Aligning new sequences against previous representatives')
        prev_db = os.path.join(outdir, 'prev_blast_db')
        run_command(
            f'{sfw.makeblastdb} -in {prev_file} -dbtype prot -out {prev_db} -parse_seqids')

        new_vs_prev_result = os.path.join(outdir, 'blast_new_vs_prev.tsv')
        run_command(
            f'{sfw.blastp} -query {new_file} -db {prev_db} -num_threads {threads} '
            f'-out {new_vs_prev_result} -evalue {evalue} -max_target_seqs {max_targets} '
            f'-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"')

        # Merge the two results
        merged_result = os.path.join(outdir, 'blast.tsv')
        with open(merged_result, 'w') as out_fh:
            # Read new vs new
            with open(new_vs_new_result, 'r') as fh:
                out_fh.write(fh.read())
            # Read new vs prev
            with open(new_vs_prev_result, 'r') as fh:
                out_fh.write(fh.read())

        return merged_result


def generate_tree(input_file, orth_list: list, outdir: str, coverage: float, evalue: float, falen: int, disable: bool = False, threads: int = 1, max_targets: int = 2000, ID: int = 0, LD: int = 0.7, AL: float = 0, AS: float = 0, aligner='diamond', clust_method='cdhit') -> Tree:
    falen -= 1  # l is the length of the word, so the actual length of the word is l+1
    orth_tree = nx.DiGraph()
    bar = tqdm(range(len(orth_list)),
               unit=f" clust iteration", disable=disable, desc=tqdm_.step(2))

    hier_dict = {}
    data_type = 'fasta'
    logger.info(f'- Hierarchical clustering...')
    for i, identity in enumerate(orth_list):
        bar.update()
        identity = round(identity, 3)
        logger.debug(
            f'- Running iterative clust: identity={identity}')

        this_outdir = f'{outdir}/clust_{identity}'

        '''
        mydict: current_repre_clust:[sub_clust1,sub_clust2,....]
        tmpdict: repre_clust:[gene1,gene2,....]
        '''
        if clust_method == 'cdhit':
            mydict, input_file = run_cdhit(
                input_file=input_file, outdir=this_outdir, id=identity, s=coverage, l=falen, threads=threads, b=20)
        elif clust_method == 'mmseqs2':
            # delete this_oudir if it exists
            if os.path.exists(this_outdir):
                shutil.rmtree(this_outdir)
            mydict, input_file = run_mmseq2(
                input_file, data_type, id=identity, coverage=coverage, outdir=this_outdir, threads=threads)
            data_type = 'index'

        # only for debug
        # tmpdict = clust_recorder(
        #     subject=tmpdict, query=mydict, tag=f'clust_{identity}')
        # with open(f'{outdir}/clust_{identity}.list', 'w') as fh:
        #     for each in tmpdict:
        #         fh.write('{}\t{}\n'.format(each, tmpdict[each]))

        need_added_node = []
        need_relabeled = {}
        for repre, sub_clusters in mydict.items():
            if i == 0:
                strains = {int(_.split(':')[0]) for _ in sub_clusters}
                members = set(sub_clusters)
                has_para = len(members) != len(strains)

                pse_repre = f'{repre}_{identity}'
                hier_dict[repre] = pse_repre
                need_added_node.append(
                    (pse_repre, {
                        'mci': identity, 'uni': identity,
                        'members': members, 'strains': strains,
                        'has_para': has_para
                    })
                )
            else:
                if len(mydict[repre]) == 1:
                    tmp_repre = hier_dict[repre]
                    orth_tree.nodes[tmp_repre]['uni'] = identity
                    need_relabeled[tmp_repre] = repre
                else:
                    pse_repre = f'{repre}_{identity}'
                    strains = set()
                    hier_node = []
                    members = set()
                    for sub_repre in sub_clusters:
                        pse_sub_repre = hier_dict[sub_repre]
                        hier_node.append((pse_repre, pse_sub_repre))
                        strains |= orth_tree.nodes[pse_sub_repre]['strains']
                        members.update(
                            orth_tree.nodes[pse_sub_repre]['members'])

                    has_para = len(members) != len(strains)
                    orth_tree.add_node(pse_repre, mci=identity, uni=identity,
                                       has_para=has_para, strains=strains, members=members)
                    orth_tree.add_edges_from(hier_node)
                    hier_dict[repre] = pse_repre
                    need_relabeled[pse_repre] = repre
        if need_added_node:
            orth_tree.add_nodes_from(need_added_node)
    nx.relabel_nodes(orth_tree, need_relabeled, copy=False)
    bar.close()

    logger.info(
        f'- Running diamond to get the ortholog node distance graph...')
    if clust_method == 'mmseqs2':
        # convert the index to fasta
        run_command(
            f'{sfw.mmseqs2} convert2fasta {input_file} {input_file}.fasta')
        input_file = f'{input_file}.fasta'
    alignment_result = run_alignment(
        input_file=input_file, outdir=outdir, threads=threads, evalue=evalue, aligner=aligner, max_targets=max_targets)
    logger.info(f'- Loading the ortholog node distance graph...')
    edges = []
    filtered_alignment_result = f'{outdir}/{aligner}.filtered.result'

    with open(alignment_result, 'r') as f:
        lines = f.readlines()

    # split the lines into chunks for parallel processing
    num_threads = threads
    chunk_size = len(lines) // num_threads
    chunks = [lines[i:i + chunk_size]
              for i in range(0, len(lines), chunk_size)]

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        results = executor.map(
            partial(process_lines, ID=ID, LD=LD, AL=AL, AS=AS), chunks)

    # collect all results
    final_lines = []
    edges = []
    for result, edge in results:
        final_lines.extend(result)
        edges.extend(edge)

    with open(filtered_alignment_result, 'w') as fh:
        fh.writelines(final_lines)

    G = nx.Graph()
    G.add_nodes_from(mydict.keys())
    G.add_weighted_edges_from(edges)
    tree = Tree()
    tree.load_alignment_result(filtered_alignment_result)  # only for mcl
    logger.info(f'- Recording the paralog node in the distance graph...')
    tree.load_ortho_identity_tree(orth_tree)
    logger.info(f'- Extracting the node relationship...')
    tree.load_distance_graph(G, raw=True)

    return tree


def run_cdhit_2d(new_file: str, ref_file: str, outdir: str, id: float, l: int = 10, s: float = 0.0, b: int = 20, threads: int = 1):
    """
    Use cd-hit-2d to compare new sequences against reference clusters.

    Args:
        new_file: Input file with new sequences to be clustered
        ref_file: Reference file with existing representative sequences
        outdir: Output directory
        id: Identity threshold
        l: Length of throw_away_sequences
        s: Length difference cutoff
        b: Band width of alignment
        threads: Number of threads

    Returns:
        mydict: Dictionary mapping representative to cluster members
        new_singletons_file: Path to FASTA file containing only new unmatched sequences (singletons)
    """
    s = s if id > 0.95 else 0
    output_file = f'{outdir}/repre_node'
    output_file_clstr = f'{output_file}.clstr'

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # Recommended word_size by cdhit manual
    if 0.7 < id <= 1:
        word_size = 5
    else:
        word_size = 4

    # Run cd-hit-2d to compare new sequences against reference
    run_command(
        f'{sfw.cdhit_2d} -c {id} -i {ref_file} -i2 {new_file} -o {output_file} -T {threads} -s {s} -M 0 -d 256 -n {word_size} -b {b} -g 1 -l {l}')

    mydict = load_cdhit_result(output_file_clstr)

    # IMPORTANT: cd-hit-2d only outputs sequences that matched to reference
    # We need to find sequences from new_file that are NOT in mydict and add them as singletons
    from Bio import SeqIO

    # Get all sequences that were assigned to clusters
    assigned = set()
    for members in mydict.values():
        assigned.update(members)

    # Get all sequences from new_file
    all_new_seqs = {}
    with open(new_file) as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            all_new_seqs[record.id] = record

    # Find unassigned sequences (didn't match any reference)
    unassigned = set(all_new_seqs.keys()) - assigned

    # Write only the new singletons (unassigned) to a new FASTA file for next iteration
    new_singletons_file = f'{outdir}/new_singletons.fa'
    with open(new_singletons_file, 'w') as fh:
        for seq_id in unassigned:
            SeqIO.write(all_new_seqs[seq_id], fh, 'fasta')

    logger.debug(f'cd-hit-2d: {len(assigned)} sequences assigned to existing clusters, '
                 f'{len(unassigned)} new singleton clusters')

    return mydict, unassigned, new_singletons_file


def run_mmseqs2_search(new_data, new_data_type: str, ref_db, outdir: str, id: float, coverage: float, threads: int = 1):
    """
    Use MMseqs2 search/map to compare new sequences against reference database.

    This is the MMseqs2 equivalent of cd-hit-2d for incremental clustering.

    Args:
        new_data: New sequences (fasta file or mmseqs2 db)
        new_data_type: Type of new_data ('fasta' or 'index')
        ref_db: Reference MMseqs2 database from previous clustering
        outdir: Output directory
        id: Identity threshold
        coverage: Coverage threshold
        threads: Number of threads

    Returns:
        mydict: Dictionary mapping representative to cluster members
        new_singletons_db: Path to database containing only new unmatched sequences (singletons)
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Create database for new sequences if needed
    if new_data_type == 'fasta':
        new_db = f'{outdir}/new_seq.db'
        run_command(
            f'{sfw.mmseqs2} createdb --shuffle 1 {new_data} {new_db} -v 0')
    else:
        new_db = new_data

    # Use mmseqs2 search to map new sequences to reference clusters
    # search is more sensitive than map and suitable for clustering
    result_db = f'{outdir}/search_result'
    tmp_dir = f'{outdir}/tmp'

    run_command(
        f'{sfw.mmseqs2} search {new_db} {ref_db} {result_db} {tmp_dir} '
        f'--min-seq-id {id} -c {coverage} --threads {threads} -v 0 --cov-mode 0')

    # Create a clustering result by assigning new sequences to their best hit
    clust_db = f'{outdir}/seq.clst'

    # Use result2profile or result2repseq to get cluster assignments
    # We'll use createtsv to export the results
    result_tsv = f'{outdir}/search_result.tsv'
    run_command(
        f'{sfw.mmseqs2} convertalis {new_db} {ref_db} {result_db} {result_tsv} '
        f'--format-output query,target -v 0')

    # Parse the results to create mydict
    mydict = {}
    assigned = set()

    with open(result_tsv) as fh:
        for line in fh:
            query, target = line.strip().split('\t')
            if query not in assigned:
                if target not in mydict:
                    mydict[target] = []
                mydict[target].append(query)
                assigned.add(query)

    # Also need to handle sequences that didn't match (create singleton clusters)
    # Get all sequence IDs from new data. Prefer parsing the original FASTA when
    # available; if we only have an mmseqs2 DB, convert it to FASTA first.
    all_seqs = set()
    try:
        if new_data_type == 'fasta':
            fasta_path = new_data
        else:
            # convert mmseqs2 db to fasta
            fasta_path = f'{outdir}/new_db.fasta'
            run_command(f'{sfw.mmseqs2} convert2fasta {new_db} {fasta_path}')

        with open(fasta_path) as fh:
            for line in fh:
                if line.startswith('>'):
                    seq_id = line[1:].split()[0]
                    all_seqs.add(seq_id)
    except Exception as e:
        logger.warning(f'Failed to read/convert new sequence headers: {e}')
        # As a fallback, try to use the assigned set only
        all_seqs = set(assigned)

    # Add unassigned sequences as singletons
    unassigned = all_seqs - assigned

    # Create output database containing ONLY new singletons (unassigned sequences)
    # This will be used as input for the next iteration
    new_singletons_db = f'{outdir}/new_singletons.db'

    # Write unassigned sequence IDs to a file
    unassigned_file = f'{outdir}/unassigned_list.txt'
    with open(unassigned_file, 'w') as fh:
        for seq_id in unassigned:
            fh.write(f'{seq_id}\n')

    # Create subdb with only unassigned sequences from new_db
    if unassigned:
        run_command(
            f'{sfw.mmseqs2} createsubdb {unassigned_file} {new_db} {new_singletons_db} -v 0')
    else:
        # If no unassigned, create an empty database
        logger.warning(
            'All sequences were assigned to existing clusters; no new singletons created.')
        empty_fasta = f'{outdir}/empty.fa'
        with open(empty_fasta, 'w') as fh:
            pass  # Empty file
        run_command(
            f'{sfw.mmseqs2} createdb {empty_fasta} {new_singletons_db} -v 0')

    logger.debug(f'mmseqs2 search: {len(assigned)} sequences assigned to existing clusters, '
                 f'{len(unassigned)} new singleton clusters')

    return mydict, unassigned, new_singletons_db


def merge_tree(input_file, previous_dir, orth_list: list, outdir: str, coverage: float, evalue: float, falen: int, disable: bool = False, threads: int = 1, max_targets: int = 2000, ID: int = 0, LD: int = 0.7, AL: float = 0, AS: float = 0, aligner='diamond', clust_method='cdhit') -> Tree:
    """
    Merge new sequences into existing tree using incremental clustering.

    This function uses the clustering results from previous_dir and applies:
    - cd-hit-2d (for cdhit method) or
    - mmseqs2 search (for mmseqs2 method)
    to quickly assign new sequences to existing clusters, avoiding full re-clustering.
    """
    falen -= 1  # l is the length of the word, so the actual length of the word is l+1

    logger.info(
        f'- Merging new sequences with previous clustering results from {previous_dir}')

    # Load previous tree structure from pickle
    import pickle
    from pgap2.lib.pklcheck import PklCheck

    with open(f'{previous_dir}/preprocess.pkl', 'rb') as fh:
        previous_pkl: PklCheck = pickle.load(fh)
        previous_tree = previous_pkl.data_dump('tree')

    # Get the orthology tree structure from previous run
    orth_tree = previous_tree.orth_identity_tree.copy()

    bar = tqdm(range(len(orth_list)),
               unit=f" clust iteration", disable=disable, desc=tqdm_.step(2))

    hier_dict = {}
    data_type = 'fasta'  # For mmseqs2

    if clust_method == 'cdhit':
        logger.info(f'- Incremental hierarchical clustering using cd-hit-2d...')
    elif clust_method == 'mmseqs2':
        logger.info(
            f'- Incremental hierarchical clustering using mmseqs2 search...')
    else:
        logger.error(f'Unsupported clustering method: {clust_method}')
        raise ValueError(
            f'clust_method must be "cdhit" or "mmseqs2", got "{clust_method}"')

    for i, identity in enumerate(orth_list):
        bar.update()
        identity = round(identity, 3)

        if clust_method == 'cdhit':
            logger.debug(
                f'- Running iterative clust with cd-hit-2d: identity={identity}')
        else:
            logger.debug(
                f'- Running iterative clust with mmseqs2 search: identity={identity}')

        this_outdir = f'{outdir}/clust_{identity}'
        prev_outdir = f'{previous_dir}/clust_{identity}'

        # Check if previous clustering result exists and run appropriate method
        if clust_method == 'cdhit':
            prev_repre_file = f'{prev_outdir}/repre_node'
            if not os.path.exists(prev_repre_file):
                logger.error(
                    f'Previous clustering result not found: {prev_repre_file}')
                logger.warning(f'Cannot proceed with incremental clustering.')
                assert False, "Previous clustering result not found, cannot proceed with incremental clustering."

            if not os.path.exists(this_outdir):
                os.makedirs(this_outdir)

            # Use cd-hit-2d to compare new sequences against previous representatives
            logger.debug(
                f'  Comparing new sequences against {prev_repre_file}')
            mydict, unassigned, output_file = run_cdhit_2d(
                new_file=input_file,
                ref_file=prev_repre_file,
                outdir=this_outdir,
                id=identity,
                s=coverage,
                l=falen,
                threads=threads,
                b=20
            )
            input_file = output_file  # Update for next iteration

        elif clust_method == 'mmseqs2':
            prev_repre_db = f'{prev_outdir}/seq.clst.rep'
            if not os.path.exists(prev_repre_db):
                logger.error(
                    f'Previous clustering database not found: {prev_repre_db}')
                logger.warning(f'Cannot proceed with incremental clustering.')
                assert False, "Previous clustering database not found, cannot proceed with incremental clustering."

            if os.path.exists(this_outdir):
                shutil.rmtree(this_outdir)
            os.makedirs(this_outdir)

            # Use mmseqs2 search to compare new sequences against previous database
            logger.debug(f'  Searching new sequences against {prev_repre_db}')
            mydict, unassigned, output_db = run_mmseqs2_search(
                new_data=input_file,
                new_data_type=data_type,
                ref_db=prev_repre_db,
                outdir=this_outdir,
                id=identity,
                coverage=coverage,
                threads=threads
            )
            input_file = output_db  # Update for next iteration
            data_type = 'index'  # Now it's a mmseqs2 database

        # Update the orthology tree with new sequences
        # Similar logic to generate_tree:
        # - First iteration (i==0): new sequences map to old representatives -> update old node's members
        #                           unassigned -> create new leaf nodes
        # - Later iterations (i>0): if multiple nodes merge -> create new parent node with edges
        #                           if single node -> just update uni and prepare for relabel

        need_added_node = []
        need_relabeled = {}

        # Helper function to propagate members up to all ancestors
        def propagate_to_ancestors(node, members, strains):
            """Propagate new members and strains to all ancestor nodes"""
            for ancestor in nx.ancestors(orth_tree, node):
                orth_tree.nodes[ancestor]['members'].update(members)
                orth_tree.nodes[ancestor]['strains'].update(strains)
                orth_tree.nodes[ancestor]['has_para'] = len(
                    orth_tree.nodes[ancestor]['members']) != len(orth_tree.nodes[ancestor]['strains'])

        if i == 0:
            # First iteration: map new sequences to old representatives, create leaf nodes for unassigned

            # Process mapped sequences - update existing nodes
            for repre, new_members in mydict.items():
                strains = {int(_.split(':')[0]) for _ in new_members}
                members = set(new_members)
                pse_repre = f'{repre}_{identity}'

                # Find existing node (could be 'repre_identity' or just 'repre' if relabeled)
                if pse_repre in orth_tree.nodes:
                    existing_node = pse_repre
                elif repre in orth_tree.nodes:
                    existing_node = repre
                else:
                    raise KeyError(
                        f"Representative {repre} not found in orth_tree")

                logger.debug(
                    f'  Updating existing node {existing_node} with {len(members)} new members')
                orth_tree.nodes[existing_node]['members'].update(members)
                orth_tree.nodes[existing_node]['strains'].update(strains)
                orth_tree.nodes[existing_node]['has_para'] = len(
                    orth_tree.nodes[existing_node]['members']) != len(orth_tree.nodes[existing_node]['strains'])

                # Propagate to ancestors
                propagate_to_ancestors(existing_node, members, strains)
                hier_dict[repre] = existing_node

            # Process unassigned - create new leaf nodes
            for seq_id in unassigned:
                strains = {int(seq_id.split(':')[0])}
                members = {seq_id}
                pse_repre = f'{seq_id}_{identity}'

                logger.debug(f'  Creating new leaf node {pse_repre}')
                hier_dict[seq_id] = pse_repre
                need_added_node.append(
                    (pse_repre, {
                        'mci': identity, 'uni': identity,
                        'members': members, 'strains': strains,
                        'has_para': False
                    })
                )

        else:
            # Later iterations: handle merging like generate_tree
            # mydict now contains: old_repre -> [new singletons from previous iteration that merged with it]

            for repre, sub_clusters in mydict.items():
                # sub_clusters are the new singleton IDs that merged with this old representative
                # The old representative (repre) already exists as a node in orth_tree

                pse_repre = f'{repre}_{identity}'

                if pse_repre in orth_tree.nodes:
                    # Case 1: The old repre already has a parent node at this identity level
                    # Just add edges from parent to new children and update parent
                    parent_node = pse_repre

                    new_strains = set()
                    new_members = set()
                    hier_node = []

                    for sub_repre in sub_clusters:
                        if sub_repre in hier_dict:
                            pse_sub_repre = hier_dict[sub_repre]
                            hier_node.append((parent_node, pse_sub_repre))
                            new_strains |= orth_tree.nodes[pse_sub_repre]['strains']
                            new_members.update(
                                orth_tree.nodes[pse_sub_repre]['members'])

                    if hier_node:
                        orth_tree.add_edges_from(hier_node)

                    orth_tree.nodes[parent_node]['members'].update(new_members)
                    orth_tree.nodes[parent_node]['strains'].update(new_strains)
                    orth_tree.nodes[parent_node]['has_para'] = len(
                        orth_tree.nodes[parent_node]['members']) != len(orth_tree.nodes[parent_node]['strains'])

                    propagate_to_ancestors(
                        parent_node, new_members, new_strains)
                    hier_dict[repre] = parent_node

                elif repre in orth_tree.nodes:
                    # Case 2: The old repre was a singleton (no parent at this identity level)
                    # Now it clusters with new sequences -> create a new parent node
                    old_node = repre

                    new_strains = orth_tree.nodes[old_node]['strains'].copy()
                    new_members = orth_tree.nodes[old_node]['members'].copy()
                    # parent -> old singleton
                    hier_node = [(pse_repre, old_node)]

                    for sub_repre in sub_clusters:
                        if sub_repre in hier_dict:
                            pse_sub_repre = hier_dict[sub_repre]
                            hier_node.append((pse_repre, pse_sub_repre))
                            new_strains |= orth_tree.nodes[pse_sub_repre]['strains']
                            new_members.update(
                                orth_tree.nodes[pse_sub_repre]['members'])

                    # Create the new parent node
                    need_added_node.append((pse_repre, {
                        'members': new_members,
                        'strains': new_strains,
                        'has_para': len(new_members) != len(new_strains),
                        'uni': identity
                    }))
                    orth_tree.add_edges_from(hier_node)
                    hier_dict[repre] = pse_repre

                else:
                    raise KeyError(
                        f"Old representative {repre} not found in orth_tree")

            # Process unassigned - these are still singletons, update their uni
            for seq_id in unassigned:
                if seq_id in hier_dict:
                    tmp_repre = hier_dict[seq_id]
                    orth_tree.nodes[tmp_repre]['uni'] = identity
                    need_relabeled[tmp_repre] = seq_id

        if need_added_node:
            orth_tree.add_nodes_from(need_added_node)

    # Relabel nodes that remained singletons throughout
    nx.relabel_nodes(orth_tree, need_relabeled, copy=False)
    bar.close()

    # The sequences that need alignment are those that remained singletons
    # throughout all clustering iterations (i.e., need_relabeled.values())
    final_new_singletons = set(need_relabeled.values())

    logger.info(
        f'- Running incremental {aligner} alignment to get the ortholog node distance graph...')
    logger.info(
        f'- Total new singleton sequences for alignment: {len(final_new_singletons)}')

    # Create a FASTA file with only the new singleton sequences for alignment
    # Read from the original input file to get sequences
    from Bio import SeqIO

    # Determine the original input file path
    # This is created by partition_add
    original_input_fa = f'{outdir}/total.involved_prot.fa'

    new_singletons_fa = f'{outdir}/new_singletons_final.fa'
    with open(new_singletons_fa, 'w') as out_fh:
        with open(original_input_fa) as in_fh:
            for record in SeqIO.parse(in_fh, 'fasta'):
                if record.id in final_new_singletons:
                    SeqIO.write(record, out_fh, 'fasta')

    logger.debug(
        f'Wrote {len(final_new_singletons)} new singleton sequences to {new_singletons_fa}')

    # Get previous representative file
    # Use the final clustering result from previous run
    prev_final_identity = orth_list[-1]
    prev_final_outdir = f'{previous_dir}/clust_{prev_final_identity}'

    if clust_method == 'cdhit':
        prev_repre_file = f'{prev_final_outdir}/repre_node'
    else:  # mmseqs2
        prev_repre_db = f'{prev_final_outdir}/seq.clst.rep'
        # Convert previous mmseqs2 db to fasta if needed
        prev_repre_file = f'{prev_final_outdir}/repre_node.fasta'
        if not os.path.exists(prev_repre_file):
            logger.debug(f'Converting previous mmseqs2 database to fasta...')
            run_command(
                f'{sfw.mmseqs2} convert2fasta {prev_repre_db} {prev_repre_file}')

    # Run incremental alignment: new singletons vs themselves + new singletons vs previous representatives
    new_alignment_result = run_incremental_alignment(
        new_file=new_singletons_fa,
        prev_file=prev_repre_file,
        outdir=outdir,
        threads=threads,
        evalue=evalue,
        aligner=aligner,
        max_targets=max_targets
    )

    # Merge with previous alignment result to get complete alignment data
    prev_alignment_result = f'{previous_dir}/{aligner}.tsv'
    merged_alignment_result = f'{outdir}/{aligner}.tsv'

    logger.info(f'- Merging new alignment with previous alignment results...')
    with open(merged_alignment_result, 'w') as out_fh:
        # Write new alignment results first
        with open(new_alignment_result, 'r') as fh:
            out_fh.write(fh.read())

        # Append previous alignment results if they exist
        if os.path.exists(prev_alignment_result):
            with open(prev_alignment_result, 'r') as fh:
                out_fh.write(fh.read())
        else:
            logger.warning(
                f'Previous alignment result not found: {prev_alignment_result}')

    logger.info(f'- Loading the ortholog node distance graph...')
    edges = []
    filtered_alignment_result = f'{outdir}/{aligner}.filtered.result'

    with open(merged_alignment_result, 'r') as f:
        lines = f.readlines()

    # Split the lines into chunks for parallel processing
    num_threads = threads
    chunk_size = max(1, len(lines) // num_threads)
    chunks = [lines[i:i + chunk_size]
              for i in range(0, len(lines), chunk_size)]

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        results = executor.map(
            partial(process_lines, ID=ID, LD=LD, AL=AL, AS=AS), chunks)

    # Collect all results
    final_lines = []
    edges = []
    for result, edge in results:
        final_lines.extend(result)
        edges.extend(edge)

    with open(filtered_alignment_result, 'w') as fh:
        fh.writelines(final_lines)

    # Load previous distance graph and add new edges
    G = previous_tree.raw_distance_graph.copy() if hasattr(previous_tree,
                                                           'raw_distance_graph') and previous_tree.raw_distance_graph else nx.Graph()

    # Add new nodes from mydict
    G.add_nodes_from(mydict.keys())

    # Add new edges
    G.add_weighted_edges_from(edges)

    # Create new tree object
    tree = Tree()
    tree.load_alignment_result(filtered_alignment_result)
    logger.info(f'- Recording the paralog node in the distance graph...')
    tree.load_ortho_identity_tree(orth_tree)
    logger.info(f'- Extracting the node relationship...')
    tree.load_distance_graph(G, raw=True)

    logger.info(f'- Successfully merged new sequences into existing tree')
    return tree
