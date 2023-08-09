#!/nfs_genome/anaconda3/envs/rnaseq/bin/python
# -*- coding: utf-8 -*-
import getopt
import pandas as pd
import re
import os
import sys
import pysam
def usage():
    print("merge the RPKG restuls for further binning)\n"
          " -h help\n" 
          " -s <fasta> a fasta sequences to get the order of the dataset\n"
          " -f <folder> Folder containing _DEPTH.tsv files (optional)\n" 
          " -o <output> merged RPKG results, default mergedDEPTH.tsv\n"
          " -v the columns for pseudo_variations are included with this option, default True that two forms of tables will be given\n"
          "")

def get_rpkg(path):
    return [os.path.join(path, f) for f in os.listdir(path) if f.endswith("_DEPTH.tsv")]
output = "mergedDEPTH"
if len(sys.argv) == 1:
    usage()
    sys.exit()

from Bio import SeqIO

def get_sequence_lengths(fasta_file):
    sequence_lengths = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence_lengths[record.id] = len(record.seq)
    return sequence_lengths

opts,args=getopt.getopt(sys.argv[1:], 'hb:s:f:o:v', ['help'])
for opt, arg in opts:
    if opt in ("-h", "--help"):
        usage()
        sys.exit()
    elif opt in ("-o"):
        output = arg
    elif opt in ("-f"):
        path = arg
    elif opt in ("-s"):
        seq = arg
    elif opt in ("-v"):
        variation = True
    else:
        usage()
        sys.exit()
rpkg_list = get_rpkg(path)

rpkg_list.sort(key=lambda x: int(re.search(r'_(\d+)_DEPTH', x).group(1)))

length_contig = get_sequence_lengths(seq)
id_index =  sorted(length_contig.keys())

for rpkg in rpkg_list:
    base=os.path.basename(rpkg)
    if rpkg_list.index(rpkg) == 0:
        rpkgDF = pd.read_csv(rpkg, sep="\t", index_col=0, skiprows=1, header=None,names=[base])
    else: 
        rpkgDF = pd.merge(rpkgDF, pd.read_csv(rpkg, sep="\t", index_col=0, skiprows=1, header=None,names=[base]), how="outer", left_index=True, right_index=True)
rpkgDF.fillna(0,inplace=True)

rpkgDF = rpkgDF.reindex(id_index).fillna(0)
mean = rpkgDF.mean(axis=1)
rpkgDF.to_csv(f"{output}.tsv", sep="\t", index=True, header=True, index_label="contigName")
if(variation):
    for col in rpkgDF.columns:
        new_col_name = col + '-var'
        new_col_values = rpkgDF[col].apply(lambda x: x*2)
        rpkgDF.insert(rpkgDF.columns.get_loc(col) + 1, new_col_name, new_col_values)
    rpkgDF.insert(0, "contigLen", rpkgDF.index.map(length_contig))
    rpkgDF.insert(1, "totalAvgDepth", list(mean))
    print(rpkgDF)
    rpkgDF.to_csv(f"{output}-pseudoVariation.tsv", sep="\t", index=True, header=True, index_label="contigName")

