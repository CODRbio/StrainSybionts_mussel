#!/nfs_genome/anaconda/envs/rnaseq/bin/python
import pandas as pd
import sys, getopt, os
from Bio import SeqIO
def usage():
    print('''
    use to statistic contig nums of each cluster
    -f feasta of the references
    -t Tau_star files, default ./Filtered_Tau_star.csv
    -o Prefix of the ouput file, default strain
    ''')

opts,args=getopt.getopt(sys.argv[1:],'-h-f:-t:-o:',['help'])
output = "strain"
Tau = "./Filtered_Tau_star.csv"
for k, v in opts:
    if k in ("-h", "--help"):
        usage()
        sys.exit()
    elif k in ("-f"):
        fasta = v
    elif k in ("-o"):
        output = v
    elif k in ("-t"):
        Tau = v
    else:
        usage()
        sys.exit()
if len(opts)>0:
    pass
else:
    usage()
    sys.exit()

df =  pd.read_csv(Tau,header=0,index_col=0)
acgt = "ACGT"
strainseqs = {}

population_idx = range(int((df.shape[1]-1)/4))
ortholog_lst = list(df.index)


def get_strain_seqs(ortholog_lst, seq_index):
    strainseqs = {}
    for strain_id in population_idx:
        strainseqs[strain_id]={}
        for seqID in ortholog_lst:
            strainseqs[strain_id][seqID] = str(seq_index[seqID].seq)
    return(strainseqs)


def update_strain_seqs(line,index,strainseqs):
    idx = [i for i in range(len(line[1:])) if line[i+1] == 1]
    seqStrain = list(map(lambda x:acgt[int(int(x)%4)],idx))
    position = line[0]
    gene_id = index
    for strain_id in population_idx:
        strainseqs[strain_id][gene_id] = strainseqs[strain_id][gene_id][:position] + seqStrain[strain_id] + strainseqs[strain_id][gene_id][(position+1):]
    return(strainseqs)

def write_strains(strainseqs, output):
    for strain in population_idx:
        seq = ""
        with open("{}.strain.{}.fa".format(output,strain),'w') as fh:
            for gene in sorted(strainseqs[strain]):
                seq = seq + strainseqs[strain][gene]
            fh.write(">strain_{strain}\n{seq}\n".format(strain=strain,seq=seq))

seq_index = SeqIO.index(fasta,'fasta')

strainseqs = get_strain_seqs(ortholog_lst, seq_index)
for index,line in df.iterrows():
    strainseqs = update_strain_seqs(line=line,index=index,strainseqs=strainseqs)
write_strains(strainseqs, output)

