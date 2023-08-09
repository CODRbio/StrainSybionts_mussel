#!/usr/bin/env python
import sys
import argparse
import os
import pandas as pd
from Bio import SeqIO

gencode = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

parser = argparse.ArgumentParser(description = "sample-specific pN/pS estimation from vcf files and genome fasta file") 
parser.add_argument("snp_table",metavar="snp_table",type=str,help="Path to the snp AD tables")
parser.add_argument("genome",metavar="genome",type=str,help="Fasta file of the genome analyzed")
parser.add_argument("sample_sel",metavar="sample_sel",type=str,help="Samples to be inclued in the analysis")

#parser.add_argument("freq",metavar="freq",type=str,help="frequency threshold for fixed snps")

args=parser.parse_args()

########################## THIS PIECE OF CODE IS MEANT TO PREPARE THE DATA STRUCTURE

# Here we create 'possible_mutations'dict, which is saving ALL possible syn and non-syn changes for each codon in the genecode
possible_mutations = {}

nts = 'AGCT'

for codon in gencode:
    codon_list = []
    possible_mutations[codon] = {}
    possible_mutations[codon]['syn']= 0
    possible_mutations[codon]['non_syn']= 0
    codon_list = [e[0] for e in codon]
    for i in range(0, len(codon_list)):
        for nt in nts:
            new_codon_list =[e[0] for e in codon_list]
            new_codon_list[i] = nt
            if new_codon_list != codon_list:
                new_codon = ''.join(new_codon_list)
                if gencode[codon] == gencode[new_codon]:
                    possible_mutations[codon]['syn'] +=1
                else:
                    possible_mutations[codon]['non_syn']+= 1

#This dictionary is meant to store the two states that we can find for each of the mutated positions
# 1st_key: 'codon_position*position in the codon';  2nd.key: state; value: codon
genome_positions = {}
#Here, we store the REFERENCE GENOME
genome = {}
#Here, we store the possible mutations for each position in our genome, for genes or genomewise
genome_wise_possible = {'syn':0, 'non_syn':0}
per_gene_possible = {}
bad_genes = []
bad_position = []
genes = SeqIO.parse(open(args.genome),'fasta')
for fasta in genes:
    [header, seq] = [fasta.id, fasta.seq]
    genome[header] = seq
    if len(seq) % 3:
        bad_genes.append(header)
        continue
    per_gene_possible[header] = {}
    per_gene_possible[header]['syn'] = 0
    per_gene_possible[header]['non_syn'] =0
    position = 1
    genome_positions[header] = {}
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if codon not in gencode:
            for cnt in range(1,4):
                bad_position.append(str(position) + '*' + str(cnt))
            continue
        genome_wise_possible['syn'] += possible_mutations[codon]['syn']
        genome_wise_possible['non_syn'] += possible_mutations[codon]['non_syn']
        per_gene_possible[header]['syn'] += possible_mutations[codon]['syn']
        per_gene_possible[header]['non_syn'] += possible_mutations[codon]['non_syn']
        for x in range(0, 3):
            genome_positions[header][str(position) + '*' + str(x+1)] = {}
            genome_positions[header][str(position) + '*' + str(x+1)]['state1'] = codon
        position +=1



########################################################### FROM NOW ON WE ARE READING OUR DATA



def code_position(POS):  #getting the codon annotation tags
    codon_position = int(int(POS-1)/3)+1
    if int(POS)%3 == 0:
        code_tag = str(codon_position) + "*3"
    else:
        code_tag = str(codon_position) + "*" + str(int(POS)%3)
    return(code_tag)


def pos2codon_ref(x):    #getting the codons in refs
    header = x["CHROM"]
    tag = x["code_tag"]
    if header in genome_positions:
        if tag in genome_positions[header]:
            return("".join(genome_positions[header][tag]['state1']))
        else:
            return(pd.NaT)
    else:
        return(pd.NaT)

def pos2codon_alt(x):    #getting the codons in alts
    alt_codon = []
    header = x["CHROM"]
    ref_codon = x["ref_codon"]
    ref_list = list(ref_codon)
    alt = x["ALT"].split(",")
    pos = int(x["POS"])%3 -1
    for item in alt:
        if item in "ATGC":
            ref_list[pos] = item
            alt_codon.append("".join(ref_list))
    return(",".join(alt_codon))

def nonsym_check(x):   #checking wheter synonymous mutations
    non_sym = []
    alt_codon = x["alt_codon"].split(",")
    ref_codon = x["ref_codon"]
    for item in alt_codon:
        if gencode[item] == gencode[ref_codon]:
            non_sym.append('0')
        else:
            non_sym.append('1')
    return(",".join(non_sym))


# Sum the count of REF and ALT positions, which will be added to a new column "count_sum"
def count_AD(x):       
    data = x[4:-1]
    total = {}
    for item in data:
        count = 0
        num = item.split(",")
        for alt_item in num:
            if len(total) > count:
                total[count] += int(alt_item)
            else:
                total[count] = int(alt_item)
            count+=1
    res=str(total[0])
    for cnt in range(1, len(total)):
        res=res+","+str(total[cnt])
    return(res)

# Recorder the alternate snp for positons with more than 1 variations
def adj_data(x):
    countData = list(map(int,x["count_sum"].split(",")))
    refO = x["REF"]
    altO = x["ALT"]
    countO = x["count_sum"]
    if len(countData) > 2:
        if float(countData[0]) / float(sum(countData)) < 0.01:
            max = 0 
            for i in range(0, len(countData)):
                if countData[i] > countData[max]:
                    max = i
            tmp = countData[0]
            countData[0] = countData[max]
            countData[max] = tmp
            countDataP = [str(i) for i in countData]
            countN = ",".join(countDataP)
            refO = x["REF"]
            altO = x["ALT"]
            altL = altO.split(",")
            refN = altL[max-1]
            altL[max - 1] = refO
            altN = ",".join(altL)
            return(countN, refN, altN)
        else:
            return(countO, refO, altO)
    else:
        return(countO, refO, altO)


# filter snps occuring only in  one sample
def discard_row(x):
    data = list(map(int,x.split(",")))
    total = sum(data)
    select = 0
    for item in data:
        if float(item)/float(total) >0.025:
            select+=1
    if select<2:
        return(True)
    else:
        return(False)

# filter polymorphism sites, if the snp frequency is very low, we will delete such snps from our analysis
def good_whether(x):
    data = list(map(int,x.split(",")))
    total = sum(data) 
    #print(x)
    #print(total)
    select = 0
    select_list = []
    for item in data:
        if float(item)/float(total) >0.01:
            select+=1
    if select>=2:
        for item in data[1:]:
            if item / (item + data[0]) >0.01:
                select_list.append("1")
            else:
                select_list.append("0")
    else:
        for item in data[1:]:
            select_list.append("0")
    return(",".join(select_list))
# find the genes with only fixed snps in each individual

def divergent_whether(row):
    polymer_count = list(map(lambda element:len([x for x in element.split(",") if int(x)>0]),row))
    fixed_check = list(map(lambda x: x<2,polymer_count))
    return(all(fixed_check))


def nonsym_count(x):
    non = 0
    sym = 0
    status = x['no_sym'].split(",")
    good = x['good'].split(",")
    for cnt in range(0, len(status)):
        if int(status[cnt]):
            if int(good[cnt]):
                non += 1
    return(non)

def sym_count(x):
    non = 0
    sym = 0
    status = x['no_sym'].split(",")
    good = x['good'].split(",")
    for cnt in range(0, len(status)):
        if int(status[cnt]) == 0:
            if int(good[cnt]):
                sym += 1
    return(sym)


samples = []


# This one is obvious...
df = pd.read_table(args.snp_table,sep="\t")
title_com = ['CHROM','POS','REF','ALT']
try:
    sampleDF = pd.read_csv(args.sample_sel,sep="\t",header=None)
    title = title_com + list(sampleDF.iloc[0])
    newData = df[title]
    prefix = os.path.basename(args.sample_sel)
except:
    newData = df
    prefix = 'all'
print(prefix)
print(newData)
#print(newData.iloc[0:2])
newData["code_tag"] = newData["POS"].apply(code_position)

newData['count_sum'] = newData.apply(count_AD,axis=1)
newData['NOexpression'] = newData.apply(lambda x: sum(map(int,x['count_sum'].split(","))) <1,axis=1)
newData.drop(index = newData.index[newData['NOexpression']],inplace=True)
newData["tmp_adj"] = newData.apply(adj_data,axis=1)
newData["count_sum"] = newData["tmp_adj"].apply(lambda x: x[0])

newData["ALT"] = newData["tmp_adj"].apply(lambda x: x[2])
newData["REF"] = newData["tmp_adj"].apply(lambda x: x[1])

newData["ref_codon"] = newData.apply(pos2codon_ref,axis=1)
#print(newData.iloc[0:2])
newData.drop(index = newData[newData["ref_codon"].isnull()].index, inplace=True)
df.drop(index = newData[newData["ref_codon"].isnull()].index, inplace=True)
newData["alt_codon"] = newData.apply(pos2codon_alt,axis=1)
#print(newData.iloc[0:2])

newData['no_sym'] = newData.apply(nonsym_check,axis=1)
#print(newData.iloc[:,4:-6])
newData["fixed"] = newData.iloc[:,4:-7].apply(divergent_whether,axis=1)
newData['good'] = newData['count_sum'].apply(good_whether)
newData.drop(index=newData[newData['count_sum'].apply(discard_row)].index,inplace=True)
newData['nonSym_count'] = newData.apply(nonsym_count,axis=1)
newData['Sym_count'] = newData.apply(sym_count,axis=1)


print("The reconstructed data matrix\n")
print(newData.iloc[0:10])

sym_genome = newData['Sym_count'].sum() 
nonsym_genome = newData['nonSym_count'].sum()

pn_ps = (float(nonsym_genome)/float(genome_wise_possible['non_syn'])) / (float(sym_genome)/float(genome_wise_possible['syn']))


with open(f"05.statistics/{prefix}_pn_ps_genome.tsv", 'w') as pn:
    pn.write("pn/ps at genonme level\n")
    pn.write(str(pn_ps))

# gene-wise pN/pS

genewise_df = newData[["CHROM","REF","Sym_count","nonSym_count"]]
Gene_stat = genewise_df.groupby("CHROM").sum()
Gene_stat["Sym_count"] = Gene_stat["Sym_count"].apply(lambda x: x+1 if x==0 else x)


#Gene_stat["dn_ds"] = (float(Gene_stat["nonSym_count"])/float(genome_wise_possible['non_syn'])) / (float(Gene_stat["Sym_count"])/float(genome_wise_possible['syn']))


print("Gene wise matrix\n")
print(Gene_stat.iloc[0:5])
output = open(f'05.statistics/{prefix}_pn_ps_genes.tsv', 'w')
output.write('gene\tsyn_found\tsyn_possible\tnon_syn_found\tnon_syn_possible\tpN/pS_found\tpN/pS_possible\tgene_length\n')

for index,row in Gene_stat.iterrows():
    gene = index
    syn_found = row["Sym_count"]
    non_syn_found = row["nonSym_count"]
    pn_ps = (float(Gene_stat.loc[gene,'nonSym_count'])/float(per_gene_possible[gene]['non_syn'])) / (float(Gene_stat.loc[gene,'Sym_count'])/float(per_gene_possible[gene]['syn']))
    pn_ps_possible = float(per_gene_possible[gene]['non_syn'])/float(per_gene_possible[gene]['syn'])
    output.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(gene, str(syn_found), str(per_gene_possible[gene]['syn']), str(non_syn_found), str(per_gene_possible[gene]['non_syn']),str(pn_ps),str(pn_ps_possible),str(len(genome[gene]))))

output.close()



divData = newData.loc[newData['fixed']]
d_sym_genome = divData['Sym_count'].sum()
d_nonsym_genome = divData['nonSym_count'].sum()
dn_ds = (float(d_nonsym_genome)/float(genome_wise_possible['non_syn'])) / (float(d_sym_genome)/float(genome_wise_possible['syn']))
with open(f"05.statistics/{prefix}_dn_ds_genome.tsv", 'w') as pn:
    pn.write("dn/ds at genonme level\n")
    pn.write(str(dn_ds))
d_gene_df = divData[["CHROM","REF","Sym_count","nonSym_count"]]
Gene_stat = d_gene_df.groupby("CHROM").sum()
Gene_stat["Sym_count"] = Gene_stat["Sym_count"].apply(lambda x: x+1 if x==0 else x)

output = open(f'05.statistics/{prefix}_dn_ds_genes.tsv', 'w')
output.write('gene\tsyn_found\tsyn_possible\tnon_syn_found\tnon_syn_possible\tpN/pS_found\tpN/pS_possible\tgene_length\n')
for index,row in Gene_stat.iterrows():
    gene = index
    syn_found = row["Sym_count"]
    non_syn_found = row["nonSym_count"]
    pn_ps = (float(Gene_stat.loc[gene,'nonSym_count'])/float(per_gene_possible[gene]['non_syn'])) / (float(Gene_stat.loc[gene,'Sym_count'])/float(per_gene_possible[gene]['syn']))
    pn_ps_possible = float(per_gene_possible[gene]['non_syn'])/float(per_gene_possible[gene]['syn'])
    output.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(gene, str(syn_found), str(per_gene_possible[gene]['syn']), str(non_syn_found), str(per_gene_possible[gene]['non_syn']),str(pn_ps),str(pn_ps_possible),str(len(genome[gene]))))

output.close()

