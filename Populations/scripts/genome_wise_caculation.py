#!/usr/bin/env python
import sys
import re
import math
import argparse
from Bio import SeqIO


###### Argparse
parser = argparse.ArgumentParser(description = "This script estimates genome-wise pi and fst by using structre.r output. It ouputs three files: per-sample pi table, pairwise fst table and matrix", epilog = 'Please, report bugs to dpicazo@ifam.uni-kiel.de')
parser.add_argument("input_file",metavar="input_file",type=str,help="structure.r output")
parser.add_argument("genome",metavar="genome",type=str,help="Reference genome used for SNP calling")
parser.add_argument("output",metavar="output",type=str,help="Output given name")
args=parser.parse_args()


###########################################################################################
##### Here, we read the genes analyzed
genes_length = {}
genes_estimates= {}
genome_length = 0

genes = SeqIO.parse(open(args.genome),'fasta')
for fasta in genes:
    [name, sequence] = [fasta.id, fasta.seq]
    gene_length = len(sequence)
    genes_length[name] = gene_length
    genome_length += int(gene_length)
    genes_estimates[name] = []

print(f"genome:{genome_length}\n")
#print genes_length
############################################################################################
##### Here, we read the info in the output from structure.R

#out = open(pi_out, 'w')

#header_complete = []
header = []
start = -1
between = {}
within = {}
within_position = {}
between_position = {}
samples = []
genes_order = []
genome_wise = []
#counter = 0#number of pi's in the file (within or inter)
with open(args.input_file, 'r') as pis:
    for line in pis:
        line = line.strip().split('\t')
        if line[1] == 'Pos':
            for i in range(0, len(line)):
                if line[i] != 'Contig' and line[i] != 'Pos' and line[i] != 'Ref' and not line[i].endswith(('.A','.T','.G','.C')):
                    start = i
                    break
            print(f"start is {start}")
            header = line[start:-1]
            #print(header)
            counter = 0
            for element in header:
                spl = element.split('__')
                if len(spl) == 2:
                    within[spl[1]] = None
                    within_position[spl[1]] = counter
                    samples.append(spl[1])
                else:
                    if not spl[1] in between:
                        between[spl[1]] = {}
                        between[spl[1]][spl[2]] = None
                        between_position[spl[1]] = {}
                        between_position[spl[1]][spl[2]] = counter
                    else:
                        between[spl[1]][spl[2]] = None
                        between_position[spl[1]][spl[2]] = counter
                counter +=1
        else:
            gene = line[0]
            lst = [float(i) for i in line[start:-1]]
            genome_wise.append(lst)


## Here, we sum up for each sample all positions and divide by the genome length
sumindex = [sum(elts) for elts in zip(*genome_wise)]

new_list = [x/genome_length for x in sumindex]


if len(new_list) >0:
    new_genome_wise = new_list
	

#print(new_genome_wise)
#for gene in new_genes_estimates:
values = new_genome_wise

for i in range(0, len(values)):
	for key in within_position:
		if within_position[key] == i:
			within[key] = new_genome_wise[i]
	for key1 in between_position:
		for key2 in between_position[key1]:
			if between_position[key1][key2] == i:
				between[key1][key2] = new_genome_wise[i]

old_within = []
young_within = []
intermediate_within = []

old_between = []
young_between = []


out2 = open('05.statistics/%s_pi.txt'%args.output, 'w')
out2.write('sample\tintra_pi\n')


out = open('05.statistics/%s_fst.txt'%args.output, 'w')
out.write('samples\tfst\n')

out3 = open('05.statistics/%s_fst.mat'%args.output, 'w')

fst_comp = {}

for a in range(0, len(samples)):
    print("%s,%s"%(samples[a], within[samples[a]]))
    fst_comp[samples[a]] = {}
    for b in range(a+1, len(samples)):
        comparisons = []
        total_fst = 0
        d= float(within[samples[a]])
        e= float(within[samples[b]])
        f= float(between[samples[a]][samples[b]])
        if f != 0:
            fst = 1-((d+e)/2/f)
        else:
            fst = 0
        if fst >1:
            print(fst)
        fst_comp[samples[a]][samples[b]] = fst
		#comparison_type = sample_sizes[samples[a]] + '_' + sample_sizes[samples[b]]
        samp = samples[a] + '_' + samples[b]
		#spots_comp = spots[samples[a]] + '_' + spots[samples[b]]
        out.write(samp+ '\t' + str(fst) + '\n')
        out2.write(samples[a] + '\t' + str(within[samples[a]]) + '\n')

for i in range(0, len(samples)):
	out3.write(',%s'%samples[i])
out3.write('\n')


for i in range(0, len(samples)):
	out3.write(samples[i])
	for x in range(0, len(samples)):
		if samples[i] == samples[x]:
			out3.write(',0')
		elif samples[i] in fst_comp  and samples[x] in fst_comp[samples[i]]:
			out3.write(',%s'%fst_comp[samples[i]][samples[x]])
		else:
			out3.write(',%s'%fst_comp[samples[x]][samples[i]])
	out3.write('\n')


out.close()
out2.close()
out3.close()
