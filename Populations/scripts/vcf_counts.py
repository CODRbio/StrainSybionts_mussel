#!/usr/bin/env python
import sys
import os
import re
import argparse
import pandas as pd
import sys

parser = argparse.ArgumentParser(description = "Merge snv data from individual vcf files into a csv table")
parser.add_argument("count", metavar="count", type = str, help="path of all.genotyped.count.tsv")
args = parser.parse_args()


count_file = args.count

counts={}
all_positions={}
if not os.path.exists("05.statistics/snv_count/"):
    os.makedirs("05.statistics/snv_count/")


def read_file(table, all_positions, counts):
    df = pd.read_table(table,sep="\t")
    for column in list(df)[4:]:
        sample = column.strip('.AD')
        counts[sample] = {}

    col = df.shape[1]
    row = df.shape[0]
    for row_c in range(0,row):
        ref = df.iloc[row_c,2].strip() 
        alt = df.iloc[row_c,3].split(',')
        chr_pos = "--".join([str(df.iloc[row_c,0]),str(df.iloc[row_c,1])])
        for idx,data in df.iloc[row_c,4:].items():
            sample = idx.strip('.AD')
            ref_c = data.split(',')[0]
            counts[sample]
            counts[sample][chr_pos] = {'A':0, 'C':0, 'G':0, 'T':0, '*':0}
            all_positions[chr_pos] = ref
            lenData =  len(data.split(',')) 
            counts[sample][chr_pos][ref[0]] += int(ref_c)
            for id in range(1,lenData):
                counts[sample][chr_pos][alt[id - 1]] += int(data.split(',')[id])

nt_order='ACGT'

read_file(count_file, all_positions, counts)

for sample in counts: 
    with open(f"05.statistics/snv_count/{sample}.snv.count.csv",'w') as out:
        out.write("CHROM\tPOS\tRef\tA\tC\tG\tT\n")
        for snp in all_positions:
            snp_info = snp.split('--')
            if snp in counts[sample]:
                out.write(snp_info[0] + "\t" + snp_info[1] +"\t" + all_positions[snp])
                for nt in nt_order:
                    out.write("\t"+str(counts[sample][snp][nt]))
                out.write('\n')
            else:
                out.write(snp_info[0] + "\t" + snp_info[1] +"\t" + all_positions[snp])
                for nt in nt_order:
                    if nt == all_positions[snp]:
                        out.write('\t1')
                    else:
                        out.write('\t0')
                out.write('\n')

