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


def generate_tree(input_file, orth_list: list, outdir: str, coverage: float, evalue: float, falen: int, disable: bool = False, threads: int = 1, max_targets: int = 2000, ID: int = 0, LD: int = 0.7, AL: float = 0, AS: float = 0, aligner='diamond', clust_method='cdhit') -> Tree:
    falen -= 1  # l is the length of the word, so the actual length of the word is l+1
    orth_tree = nx.DiGraph()
    bar = tqdm(range(len(orth_list)),
               unit=f" clust iteration", disable=disable, desc=tqdm_.step(2))

    hier_dict = {}
    tmpdict = {}
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
    tree.load_alignment_result(filtered_alignment_result)
    logger.info(f'- Recording the paralog node in the distance graph...')
    tree.load_ortho_identity_tree(orth_tree)
    logger.info(f'- Extracting the node relationship...')
    tree.load_distance_graph(G, raw=True)

    return tree
