import edlib
from Bio import SeqIO
from Bio.Seq import Seq
import sys


def run_pw(seqA, seqB, dna=False):
    def calculate_pwid(sA, sB, dna):
        additional_equalities = [('A', 'N'), ('C', 'N'), ('G', 'N'), ('T', 'N')] if dna else \
            [('*', 'X'), ('A', 'X'), ('C', 'X'), ('B', 'X'),
             ('E', 'X'), ('D', 'X'), ('G', 'X'), ('F', 'X'),
             ('I', 'X'), ('H', 'X'), ('K', 'X'), ('M', 'X'),
             ('L', 'X'), ('N', 'X'), ('Q', 'X'), ('P', 'X'),
             ('S', 'X'), ('R', 'X'), ('T', 'X'), ('W', 'X'),
             ('V', 'X'), ('Y', 'X'), ('X', 'X'), ('Z', 'X'),
             ('D', 'B'), ('N', 'B'), ('E', 'Z'), ('Q', 'Z')]
        aln = edlib.align(sA, sB, mode="HW", task='distance', k=0.5 *
                          len(sA), additionalEqualities=additional_equalities)
        if aln['editDistance'] == -1:
            return 0.0
        return 1.0 - aln['editDistance'] / float(len(sA))
    if len(seqA) > len(seqB):
        seqA, seqB = seqB, seqA
    if dna:
        pwid = max(calculate_pwid(sA, seqB, dna)
                   for sA in [seqA, str(Seq(seqA).reverse_complement())])
    else:
        pwid = calculate_pwid(seqA, seqB, dna)
    return pwid


fasta_file = sys.argv[1]

query_id = sys.argv[2]
target_id = sys.argv[3]

# read the fasta file and find the sequences for the query and target
target_sequence = []
for sequence in SeqIO.parse(fasta_file, "fasta"):
    if sequence.id in (target_id, query_id):
        target_sequence.append((sequence.seq, sequence.id))
    if len(target_sequence) == 2:
        pwid = run_pw(target_sequence[0][0], target_sequence[1][0])
        print('{}\t{}\t{}'.format(
            target_sequence[0][1], target_sequence[1][1], pwid))
        break
