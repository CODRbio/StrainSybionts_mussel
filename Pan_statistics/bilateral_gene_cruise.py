#!/usr/bin/env python
import gffutils
import os
import pybedtools 
import pandas as pd
import re
import sys, getopt
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from scipy import stats
from scipy import mean

window = 3000
def usage():
    print('''
    use to compare the gene_density around core genes and strain_specific genes
    -O Orthogroup table
    -f merged fasta file of all the refs 
    -s strain-specific gene list
    -g merged gff file of refs
    -c ortholog catagories included (a list of Ortholog of interest, for example strain-specific orhologs, no title should be included in the file)
    -w window size, default 3000
    ''')

opts,args=getopt.getopt(sys.argv[1:],'-h-f:-c:-t:-O:-s:-w:-g:',['help'])
for k, v in opts:
    if k in ("-h", "--help"):
        usage()
        sys.exit()
    elif k in ("-O"):
        Orthogroup = v
    elif k in ("-g"):
        gff = v
    elif k in "-f":
        fasta =v
    elif k in ("-s"):
        specific = v
    elif k in ("-c"):
        catagory = v
    elif k in ("-w"):
        window = int(v)
    else:
        usage()
        sys.exit()

if len(opts)>0:
    pass
else:
    usage()
    sys.exit()

StrinSpecific = []
with open(specific,'r') as SPE:
    for item in SPE.readlines():
        StrinSpecific.append(item.strip())



db_file = "AllRef.db"

if not os.path.exists(db_file):
    db = gffutils.create_db(gff, db_file)
else:
    db = gffutils.FeatureDB(db_file)
df = pd.read_table(Orthogroup,index_col="Orthogroup", sep ="\t")
core = df[df.count(axis=1)/df.shape[1]>=0.9].index

OrthGeneDic = {}

def Ortho2Gene(x):
    pattern = r', *'
    OrthGeneDic[x.name] = []
    for cnt in range(0,len(x)):
        if not pd.isnull(x[cnt]):
            ele = re.split(pattern,x[cnt])
            OrthGeneDic[x.name].extend(ele)

from Bio import SeqIO

def get_sequence_lengths(fasta_file):
    lengths = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        lengths[record.id] = len(record.seq)
    return lengths

def gene_info(geneID, db=db):
    gene = db[geneID]
    info = str(gene).split("\t")[8]
    location = (gene.chrom, gene.start, gene.end)
    return(location)

def get_chromsome_len(gff):
    gene_list = {}
    with open(gff,"r") as GFF:
        for line in GFF.readlines():
            if line.startswith("#"):
                continue
            info = line.strip().split("\t")
            if info[2] == "gene":
                Gid = info[8] 
                matchObj = re.match(r'ID=(\S*?)_gene', Gid) #ID=LOKAKALI_00001_gene
                TarGID = matchObj.group(1)
                gene_list[TarGID] = 1
    return(gene_list)

df.apply(Ortho2Gene, axis=1)

def gene_in_right(query, ref,windows=window):
    q = pybedtools.BedTool(query)
    r = pybedtools.BedTool(ref)
    q.window(r, l=0,r=windows).moveto(f"{query}_right.bed")
    
def gene_in_left(query,ref,windows=window):
    q = pybedtools.BedTool(query)
    r = pybedtools.BedTool(ref)
    q.window(r, l=windows,r=0).moveto(f"{query}_left.bed")


def FilterDB(lists):
    print(lists)
    GENEdb = []
    with open(lists, 'r') as IN:
        for item in IN.readlines():
            item = item.strip()
            GENEdb.extend(OrthGeneDic[item])
    return(GENEdb)

print("Getting Ortho2Gene Dict\n")

df.apply(Ortho2Gene, axis=1)

print("Getting the length of Chromsome\n")

chrom_len = get_sequence_lengths(fasta)
GeneList = get_chromsome_len(gff=gff)


if "catagory" not in locals().keys():
    ortho = "All"
else:
    ortho = "Tranpose"
    print(catagory)
    GeneList = FilterDB(catagory)
   


def write_boundary(OrthoList, mark, OrthGeneDic=OrthGeneDic):
    UPP = open(f"{mark}_right.bed",'w')
    LOW = open(f"{mark}_left.bed",'w') 
    for Ortho in OrthoList: 
        geneList = OrthGeneDic[Ortho]
        for gene in geneList:
            try:
                gene_chrom, gene_start, gene_end = gene_info(gene)
                #print(f"{gene_chrom}\t{gene_start}\t{gene_end}\t{chrom_len[gene_chrom]}")
                if (int(gene_start) <= 1000) or (gene_end > (int(chrom_len[gene_chrom])-1000)):
                    continue
                else:
                    upper_end = int(gene_end) + 1
                    lower_end = int(gene_start) - 2
                    UPP.write(f"{gene_chrom}\t{upper_end}\t{upper_end}\t{gene}\n")
                    LOW.write(f"{gene_chrom}\t{lower_end}\t{lower_end}\t{gene}\n")
            except:
                continue
            else:
                pass
    UPP.close()
    LOW.close()

print("Generating the boundary files for genes of interest")


write_boundary(core, "core")
write_boundary(StrinSpecific, "Strain_specific")

print("Generating the adjacent gene lists within the window")
gene_in_left("core_left.bed", gff)
gene_in_left("Strain_specific_left.bed", gff)
gene_in_right("core_right.bed", gff)
gene_in_right("Strain_specific_right.bed", gff)

Density_r = {}
WinFinal_r = {}
def calculate_density_right(x, window=window, filterDB=GeneList):
    gene_name = x[3]
    if gene_name not in Density_r:
        Density_r[gene_name] = []
        WinFinal_r[gene_name] = []
    chrom = x[0]
    point = int(x[1])
    start = x[7]
    end = x[8]
    Gid = x[12] 
    matchObj = re.match( r'ID=(\S*?)_gene', Gid)
    TarGID = matchObj.group(1)
    if TarGID in filterDB:
        if (point + window -1) > int(chrom_len[chrom]):
            window = int(chrom_len[chrom]) + 1 -point 
        if (end - point) > window:
            end = point + window -1
        Density_r[gene_name].extend(range(start,end))
    WinFinal_r[gene_name].append(window)

Density_l = {}
WinFinal_l = {}
def calculate_density_left(x, window=window, filterDB=GeneList):
    gene_name = x[3]
    if gene_name not in Density_l:
        Density_l[gene_name] = []
        WinFinal_l[gene_name] = []
    chrom = x[0]
    point = int(x[1])
    start = x[7]
    end = x[8]
    Gid = x[12]
    matchObj = re.match( r'ID=(\S*?)_gene', Gid)
    TarGID = matchObj.group(1)
    if TarGID in filterDB:
        if (point - window) < 0:
            window = point
        if (point - start +1) > window:
            start = point - window + 1
        Density_l[gene_name].extend(range(start,end))
    WinFinal_l[gene_name].append(window)
     
df_coreL = pd.read_table("core_left.bed_left.bed",sep="\t", header=None)
df_coreL.drop(index = df_coreL[df_coreL[6]!="gene"].index,inplace=True)
df_coreR = pd.read_table("core_right.bed_right.bed",sep="\t", header=None)
df_coreR.drop(index = df_coreR[df_coreR[6]!="gene"].index,inplace=True)
df_speL = pd.read_table("Strain_specific_left.bed_left.bed",sep="\t", header=None)
df_speL.drop(index = df_speL[df_speL[6]!="gene"].index,inplace=True)
df_speR = pd.read_table("Strain_specific_right.bed_right.bed",sep="\t", header=None)
df_speR.drop(index = df_speR[df_speR[6]!="gene"].index,inplace=True)

print("Caculating the density\n")
df_coreL.apply(calculate_density_left, axis=1)
Core_d_l = Density_l
Core_w_l = WinFinal_l
Density_l = {}
WinFinal_l = {}
df_coreR.apply(calculate_density_right, axis=1)
Core_d_r = Density_r
Core_w_r = WinFinal_r
Density_r = {}
WinFinal_r = {}
df_speL.apply(calculate_density_left, axis=1)
Special_d_l = Density_l
Special_w_l = WinFinal_l
df_speR.apply(calculate_density_right, axis=1)
Special_d_r = Density_r
Special_w_r = WinFinal_r

dens_core_res = {}
COREwrite = open(f"core_{window}_{ortho}.tsv",'w')
for gene in (set(Core_d_l.keys()) | set(Core_d_r.keys())):
    if gene not in Core_d_l:
        Core_d_l[gene]=[]
    if gene not in Core_d_r:
        Core_d_r[gene]=[]
    if gene not in Core_w_r:
        Core_w_r[gene]=[0]
    if gene not in Core_w_l:
        Core_w_l[gene]=[0]
    length = len(set(Core_d_l[gene])) + len(set(Core_d_r[gene]))
    den = float(length) / float(min(Core_w_l[gene]) + min(Core_w_r[gene]) )
    COREwrite.write(f"{gene}\t{den}\n")
    dens_core_res[gene] = den
COREwrite.close()

dens_speicify_res={}
SPECIwrite =  open(f"Speicified_{window}_{ortho}.tsv",'w')
for gene in (set(Special_d_r.keys()) | set(Special_d_l.keys())):
    if gene not in Special_d_l:
        Special_d_l[gene]=[]
    if gene not in Special_d_r:
        Special_d_r[gene]=[]
    if gene not in Special_w_r:
        Special_w_r[gene]=[0]
    if gene not in Special_w_l:
        Special_w_l[gene]=[0]
    length = len(set(Special_d_l[gene])) + len(set(Special_d_r[gene]))
    den = float(length) / float(min(Special_w_l[gene]) + min(Special_w_r[gene]) )
    SPECIwrite.write(f"{gene}\t{den}\n")
    dens_speicify_res[gene] = den
SPECIwrite.close()
CORE = list(dens_core_res.values())
SPECI = list(dens_speicify_res.values())
data = [CORE, SPECI]
fig, ax = plt.subplots()
y_major_locator=MultipleLocator(0.1)
x_major_locator=MultipleLocator(1)

ax.yaxis.set_major_locator(y_major_locator)
ax.xaxis.set_major_locator(x_major_locator)

plt.ylim(0,1.1)

ax.violinplot(data,showmeans=True)

ax.set_title("gene density of different gene sets")
ax.set_xticklabels(["","CORE","SPECIFY"])



plt.savefig(f"density_{window}_{ortho}_violin.pdf")
plt.show()

print("executing statistic analysis")

if stats.levene(CORE,SPECI)[1] <= 0.05:
    equal_var = False
else:
    equal_var = True
with open(f"density_{window}_{ortho}_test.tsv",'w') as TEST:
    pvalue = stats.ttest_ind(CORE,SPECI, equal_var = equal_var)[1]
    if pvalue > 0.05:
        TEST.write(f"Mean value of CORE {mean(CORE)} is not signifcantly different from the SPECIFIED {mean(SPECI)}, pvalue is {pvalue}\n")
        print(f"Mean value of CORE {mean(CORE)} is not signifcantly different from the SPECIFIED {mean(SPECI)}, pvalue is {pvalue}\n")
    else:
        TEST.write(f"Mean value of CORE {mean(CORE)} is signifcantly different from the SPECIFIED {mean(SPECI)}, pvalue is {pvalue}\n")
        print(f"Mean value of CORE {mean(CORE)} is signifcantly different from the SPECIFIED {mean(SPECI)}, pvalue is {pvalue}\n")
