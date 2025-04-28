#!/usr/bin/env python3
# Script by Jason Kwong & Torsten Seemann
# Script to mask recombination from CFML output and draw SVG of recombinant regions

# Usage
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(
	formatter_class=RawTextHelpFormatter,
	description='Mask recombination from ClonalFrameML/Gubbins output and draw SVG of recombinant regions',
	usage='\n  %(prog)s --aln FASTA --out OUTFILE [--gubbins] <PREFIX>')
parser.add_argument('prefix', metavar='PREFIX', help='prefix used for CFML/Gubbins input files (required)')
parser.add_argument('--gubbins', action='store_true', help='parse as Gubbins instead of ClonalFrameML')
parser.add_argument('--aln', metavar='FASTA', required=True, help='multiFASTA alignment used as input for CFML (required)')
parser.add_argument('--out', metavar='OUTFILE', default='maskrc.aln', help='output file for masked alignment (default="maskrc.aln")')
parser.add_argument('--symbol', metavar='CHAR', default='?', help='symbol to use for masking (default="?")')
parser.add_argument('--regions', metavar='FILE', help='output recombinant regions to file')
parser.add_argument('--svg', metavar='FILE', help='draw SVG output of recombinant regions and save as specified file')
parser.add_argument('--svgsize', metavar='WIDExHIGH', default='800x600', help='specify width and height of SVG in pixels (default="800x600")')
parser.add_argument('--svgorder', metavar='FILE', help='specify file containing list of taxa (1 per line) in desired order')
parser.add_argument('--svgcolour', metavar='COLOUR', default='black', help='specify colour of recombination regions in HEX format (default=black)')
parser.add_argument('--consensus', action='store_true', help='add consensus row of recombination hotspots')
parser.add_argument('--version', action='version', version='%(prog)s 0.5')

args = parser.parse_args()

prefix = str(args.prefix)
treefn = prefix + ('.labelled_tree.newick', '.final_tree.tre')[args.gubbins]
coordfn = prefix + ('.importation_status.txt', '.recombination_predictions.gff')[args.gubbins]

import os
import sys
import csv
from collections import defaultdict
from ete3 import Tree
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
import svgwrite
from BCBio import GFF
import pprint

# Functions
def msg(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

def err(*args, **kwargs):
	msg(*args, **kwargs)
	sys.exit(1);

def check_file(f):
	msg('Checking file: {}'.format(f))
	if os.path.isfile(f) == False:
		err('ERROR: Cannot find "{}". Check output files exist in this directory.'.format(f))

# Check input files (output from CFML)
check_file(treefn)
check_file(coordfn)

# Setup lists and dictionary
d = defaultdict(list)
seqLIST = []
leafLIST = []
t = Tree(treefn, format=1)

# Set taxa order for SVG
if args.svgorder:
	msg('Setting taxa order using --svgorder {}'.format(args.svgorder));
	with open(args.svgorder) as file:
		for row in file:
			leafLIST.append(row.rstrip())
else:
	msg('Setting taxa order using leaves from: {}'.format(treefn));
	for leaf in t:
		leafLIST.append(leaf.name)

# gubbins GFF
# node="internal_13->2016-04538";neg_log_likelihood="3764.359443";taxa="2016-04538";snp_count="108"; 
# node="internal_ROOT->internal_10";neg_log_likelihood="47939.649932";taxa=" 2016-23066 2015-19014";snp_count="8";
# node="internal_15->internal_16";neg_log_likelihood="574.375212";taxa=" 2008-16004 2010-22387";snp_count="39";
# node="internal_15->internal_16";neg_log_likelihood="418.840895";taxa=" 2008-16004 2010-22387";snp_count="38";
# node="internal_12->internal_14";neg_log_likelihood="17258.115166";taxa=" 2012-21203  2008-20031  2008-16004 2010-22387";

# clonalframeml "BED"
#Node    		Beg     End
#2016-11776      5624    6619
#2016-11776      10778   11326
#Reference       2489729 2490893
#NODE_10 		1165936 1166465

# internal
# ['Reference', ['2434777', '2434933', 'ancestral']],
# ['2016-17702', ['2483546', '2485711', 'ancestral']],
# ['2016-17705', ['2483546', '2485711', 'ancestral']],
# ['2016-11776', ['2483546', '2485711', 'ancestral']],
# ['2016-11778', ['2483546', '2485711', 'ancestral']],
# ['2016-11777', ['2483546', '2485711', 'ancestral']],
# ['2016-17701', ['2483546', '2485711', 'ancestral']],
# ['Reference', ['2483546', '2485711', 'ancestral']]]
 
if args.gubbins:
	msg('Parsing Gubbins file: {}'.format(coordfn))
	with open(coordfn) as gff:
		for rec in GFF.parse(gff):
			for f in rec.features:
				taxa = f.qualifiers['taxa'][0].strip().split()
				state = ('ancestral', 'extant')[ len(taxa)==1 ]
				for taxon in taxa:
					# BCBio-GFF module parser seems to extract coordinate-1 for f.location.start
					seqLIST.append( [ taxon, [ int(f.location.start)+1, int(f.location.end), state ] ] )
else:
	msg('Parsing ClonalFrameML file: {}'.format(coordfn))
	with open(coordfn) as csvfile:
		RCseqs = csv.reader(csvfile, delimiter='\t')
		next(csvfile)
		for row in RCseqs:
			seq = row[0]
			RCstart = row[1]
			RCstop = row[2]
			node = t.search_nodes(name=seq)[0]
			if node.is_leaf():
				leaf_names = [node.name]
				z = 'extant'
			else:
				leaves = node.get_leaves()
				leaf_names = [leaf.name for leaf in leaves]
				z = 'ancestral'
			for l in leaf_names:
				seqLIST.append([l, [RCstart, RCstop, z]])

# Build dictionary of recombinant regions
for k,v in seqLIST:
	d[k].append(v)

# Mask recombination in sequences
seqALN = []
for record in SeqIO.parse(args.aln, 'fasta'):
	msg('Reading {} ... '.format(record.id))
	seqlen = len(record.seq)
	regions = d.get(record.id, None)
	newrec = MutableSeq(str(record.seq))
	if regions:
		for a in regions:
			start = int(a[0]) - 1
			end = int(a[1])
			lenMASK = end - start
			newrec[start:end] = (args.symbol)*lenMASK
	seqALN.append(SeqRecord(Seq(str(newrec)), record.id, description=''))

# Write masked alignment to file
msg('Writing masked alignment to {} ... '.format(args.out))
SeqIO.write(seqALN, args.out, 'fasta')

# Write recombinant regions to file
if args.regions:
	with open(args.regions, 'w') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t')
		for k,v in list(d.items()):
			for b in v:
				x = b[0]
				y = b[1]
				z = b[2]
				csvwriter.writerow([k,x,y,z])

# Draw SVG
svgsize = args.svgsize.split('x',1)
width = int(svgsize[0])
height = int(svgsize[1])
numseqs = len(leafLIST)
h = height/float(numseqs)
s = width/float(int(seqlen))
fsize = int(h*0.9)
font = 'font-size:{}px; font-family:Arial'.format(fsize)
colour = args.svgcolour
main_colour = 'black'
ancestral_colour = '#A9A9A9'
interval = 500000		# Tick interval in bp

dwg = svgwrite.Drawing(args.svg)

def rect(x,p,w,colour):		# Draw regions
	dwg.add(dwg.rect(insert=((x*s)+5, (p*h)+5), size=(w, h*0.8), fill=colour))
	if args.consensus:
		dwg.add(dwg.rect(insert=((x*s)+5, 5), size=(w, h*0.8), fill=main_colour))

def label(n,p):			# Sequence ID labels
	dwg.add(dwg.text(n, insert=((seqlen*s)+15, (((p+1)*h)-(0.2*h))+5), fill=main_colour, style=font))

def ticks(q):			# Draw tick marks
	count = interval
	while count < seqlen:
		dwg.add(dwg.line((count*s,((q)*h)), (count*s,((q)*h)-3), stroke=main_colour))
		dwg.add(dwg.text((count/1000000.0), insert=((count*s)-6, ((q)*h)+13), fill=main_colour, style=font))
		count = count + interval

def box(width,h):		# Box with 5px margin
	dwg.add(dwg.line((0,0), ((width+10),0), stroke=main_colour))
	dwg.add(dwg.line(((width+10),0), ((width+10),((p+1)*h)), stroke=main_colour))
	dwg.add(dwg.line(((width+10),((p+1)*h)), (0,((p+1)*h)), stroke=main_colour))
	dwg.add(dwg.line((0,((p+1)*h)), (0,0), stroke=main_colour))

# Draw recombinant regions and labels
if args.svg:
	msg('Drawing SVG to {} ...'.format(args.svg))
	if args.consensus:
		p = 1
		dwg.add(dwg.text('Consensus', insert=((seqlen*s)+15, (0.8*h)+5), fill=main_colour, style=font))
	else:
		p = 0
	for seq in leafLIST:
		v = d.get(seq, None)
		if v:
			for b in v:
				x = int(b[0])
				y = int(b[1])
				z = b[2]
				w = (y-x)*s
				if w < 1:
					w = 1
				if z == 'extant':
					rect(x,p,w,colour)		# regions
				elif z == 'ancestral':
					rect(x,p,w,ancestral_colour)
		label(seq,p)			# labels
		p = p+1
	# Draw bounding box
	box(width,h)
	# Add tick marks
	ticks(p+1)
	dwg.save()

# Exit
msg('Done.')
sys.exit(0)
