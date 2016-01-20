import sys
import random
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

alphabet = "ACGT"
def random_string(length):
	return Seq(''.join([alphabet[random.randint(0, len(alphabet) - 1)] for _ in xrange(length)]), generic_dna)

def revseq(seq):
	return seq.reverse_complement()

def mutate(seq, similarity):
	mutatedSeq = [' '] * len(seq)
	for i in xrange(len(seq)):
		mutatedSeq[i] = seq[i] if random.random() <= float(similarity) else alphabet[random.randint(0, len(alphabet) - 1)]
	return Seq(''.join(mutatedSeq), generic_dna)

N_MODE = 1
ORDINARY_MODE = 0

def mutateN(seq, prob_ordinary_to_n, prob_n_to_ordinary):
	mode = ORDINARY_MODE if random.random() <= 0.8 else N_MODE
	mutatedSeq = [' '] * len(seq)
	for i in xrange(len(seq)):
		mutatedSeq[i] = seq[i] if mode == ORDINARY_MODE else 'N'
		if mode == ORDINARY_MODE:
			mode = N_MODE if random.random() <= prob_ordinary_to_n else ORDINARY_MODE
		else:
			mode = ORDINARY_MODE if random.random() <= prob_n_to_ordinary else N_MODE
	return Seq(''.join(mutatedSeq), generic_dna)

base_length = 5000
repeat_length = 2000

block = [random_string(repeat_length), random_string(repeat_length), random_string(repeat_length)]
template = [
[block[0], block[1], random_string(base_length), block[2], block[1], block[1], revseq(block[2]), random_string(base_length)],
[block[1], block[2], random_string(base_length), revseq(block[0]), random_string(base_length), block[2], revseq(block[1])],
[block[2], random_string(base_length), block[2], revseq(block[1]), block[0], random_string(base_length), block[1]]
]

glist = []
for g in template:
	glist.append([mutate(b, 0.75) for b in g])

for gid, g in enumerate(glist):	
	seq = mutateN(Seq(''.join([str(x) for x in g]), generic_dna), 0.01, 0.1)
	SeqIO.write(SeqRecord(seq, id = 'genome' + str(gid)), sys.stdout, 'fasta')	