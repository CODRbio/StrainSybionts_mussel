#!/usr/bin/env python
import sys, getopt, os
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np

surfix=("fas","fa","faa","pep","fasta","aa")

def usage():
    print('''
    -h help
    -p folders containing aligned faa seqs in fasta format
    -c all cds seqs in fasta format
    -s SequenceID map table
    -t surfix of files containing Trimmed columns information from trimal, default columns
    -o output prefix
    ''')

def revAlign(cds_seq, faa_seqs): #align nucleotides based on peptides, new list with 3 polymer were constructed first
    mer3 = []
    cds_align = []
    pep_count = 0
    faa = faa_seqs.replace("-","")
    if len(faa)*3 > len(cds_seq):
        gap_add = "-" * (len(faa)*3 -len(cds_seq))
        cds_seq = cds_seq + gap_add

    for num in range(0,len(cds_seq),3):
        mer3.append(cds_seq[num:(num+3)])
    for num in range(len(faa_seqs)):
        if faa_seqs[num] == "-":
            cds_align.append("-"*3)
        else:
            cds_align.append(mer3[pep_count])
            pep_count += 1
    return(cds_align)

def trimAlign(cds_align,columns_retain):
    trimmed = np.array(cds_align)[columns_retain]
    trimmed = trimmed.tolist()
    flatten_list = [i for item in trimmed for i in item]
    seqs = "".join(flatten_list)
    return(seqs)

def get_pep(path): 
    name = []
    for files in os.listdir(path):
        if files.endswith(surfix):
            name.append(os.path.join(path,files))
    return name

def id_mapper(SeqMap):
    mapper = {}
    IN=open(SeqMap)
    for lines in IN.readlines():
        ID = lines.split()[0].replace(":","")
        SeqID = lines.split()[1]
        mapper[ID] = SeqID
    return(mapper)
trimSurfix = "columns"
opts,args=getopt.getopt(sys.argv[1:],'-h-c:-p:-s:-t:-o:',['help'])
for k, v in opts:
    if k in ("-h", "--help"):
        usage()
        sys.exit()
    elif k in ("-o"):
        output = v
    elif k in ("-c"):
        cds = v
    elif k in ("-p"):
        peptide = v
    elif k in ("-s"):
        SeqMap = v
    elif k in ("-t"):
        trimSurfix = v
    else:
        usage()
        sys.exit()

cds = SeqIO.index(cds,"fasta")
#surfix=["fas","fa","faa","pep","fasta","aa"]
surfix="aln.fa"
peps = get_pep(path=peptide)
try:
    IDmap = id_mapper(SeqMap=SeqMap)
except:
    print("No id translation was performed")
for aln_file in peps:
    print(aln_file)
    #AlnRetain = [ int(x)-1 for x in goodAlnCol ]
    align=AlignIO.read(aln_file, "fasta")
    basename=os.path.split(aln_file)[1]
    basename=os.path.splitext(basename)[0]
    colName = os.path.splitext(aln_file)[0]+"."+trimSurfix
    IN=open(colName)
    first = IN.readline()
    goodAlnCol = first.strip().split("\t")[1].split(",")
    goodAlnCol = [ int(x) for x in goodAlnCol]
    try:
        cdsname = output + basename + "_cds.fa"
        cdsadj = output + basename + "_adj_cds.fa"
    except:
        cdsname = basename + "_cds.fa"
        cdsadj = basename + "_adj_cds.fa"
    rec_list = []
    rec_adj = []
    for record in align:
        faa_seqs = record.seq
        #print(record.id)
        try:
            seq_id = IDmap[record.id]
        except:
            seq_id = record.id
        print(seq_id)
        cds_seqs = cds[seq_id].seq
        #print(cds_seqs)
        #print(faa_seqs)
        cds_aln = revAlign(cds_seq=str(cds_seqs), faa_seqs=str(faa_seqs))
        cds_trim = trimAlign(cds_align=cds_aln,columns_retain=goodAlnCol)
        # print(cds_aln)
        rec_new = SeqRecord(Seq(cds_trim),id=seq_id,description="aligned cds of " + seq_id)
        rec_tre = SeqRecord(Seq(cds_trim),id=record.id)
        rec_list.append(rec_new)
        rec_adj.append(rec_tre)
    SeqIO.write(rec_list, cdsname, "fasta")
    SeqIO.write(rec_adj, cdsadj, "fasta")
