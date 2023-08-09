#!/nfs_genome/anaconda/envs/lathe/bin/python3
import sys, getopt, os


def usage():
    print('''
    use to distincuish the core, soft core and shell genes and compare the pNpS values between the gene categories
    -O Orthogroup gene count table
    -s species_id.txt
    -S strain specific list
    -p pN/pS results (optional)
    -d dN/dS result (optional)
    -o output prefix
    ''')

opts,args=getopt.getopt(sys.argv[1:],'-h-d:-p:-t:-O:-s:-S:-o:',['help'])
pNpS = False
dNdS = False
specified = False
for k, v in opts:
    if k in ("-h", "--help"):
        usage()
        sys.exit()
    elif k in ("-O"):
        Orthogroup = v
    elif k in ("-s"):
        speMap = v
    elif k in ("-S"):
        specified = v
    elif k in ("-d"):
        dNdS = v
    elif k in ("-p"):
        pNpS  = v
    elif k in ("-o"):
        prefix  = v
    else:
        usage()
        sys.exit()
if len(opts)>0:
    pass
else:
    usage()
    sys.exit()

import pandas as pd
import numpy as np
import scipy.stats as stats
from numpy import mean
import seaborn as sns
def countOrthoFre(m):
    total = len(m) - 1
    present = len([i for i in m[0:-1] if i>0])
    return(float(present)/float(total))

df = pd.read_table(Orthogroup,sep="\t")
df.set_index("Orthogroup",inplace=True)
df["Frequcy"] = df.apply(countOrthoFre, axis=1)

Core_hard = list(df[df["Frequcy"] >= 0.99].index)
Core_soft = list(df[(df["Frequcy"]>=0.95) & (df["Frequcy"] <0.99)].index)
Core_all = list(df[df["Frequcy"]>=0.90].index)
shell = df[(df["Frequcy"]>=0.15) & (df["Frequcy"] <0.90)].index
cloud = list(df[df["Frequcy"] < 0.15].index)
StrinSpecific = []
if(specified):
    with open(specified,'r') as SPE:
        for item in SPE.readlines():
            StrinSpecific.append(item.strip())
if(pNpS):
    pn = pd.read_table(pNpS,index_col="gene", sep ="\t")
    pn_coreH = list(set([i+".fa.aln" for i in Core_hard]).intersection(set(list(pn.index))))
    pn_coreA = list(set([i+".fa.aln" for i in Core_all]).intersection(set(list(pn.index))))
    pn_speicify = list(set([i+".fa.aln" for i in StrinSpecific]).intersection(set(list(pn.index))))
    pn_coreS = list(set([i+".fa.aln" for i in Core_soft]).intersection(set(list(pn.index))))
    pn_shell = list(set([i+".fa.aln" for i in shell]).intersection(set(list(pn.index))))
    pn_cloud = list(set([i+".fa.aln" for i in cloud]).intersection(set(list(pn.index))))
    pn.loc[pn_coreH,"pN/pS_found"].to_csv(f"{prefix}_pnps_coreH.csv")
    pn.loc[pn_coreS,"pN/pS_found"].to_csv(f"{prefix}_pnps_coreS.csv")
    pn.loc[pn_shell,"pN/pS_found"].to_csv(f"{prefix}_pnps_shell.csv")
    pn.loc[pn_cloud,"pN/pS_found"].to_csv(f"{prefix}_pnps_cloud.csv")
    pn.loc[pn_coreA,"pN/pS_found"].to_csv(f"{prefix}_pnps_coreA.csv")
    pn.loc[pn_speicify,"pN/pS_found"].to_csv(f"{prefix}_pnps_strainSPE.csv")
if(dNdS):
    dn =pd.read_table(dNdS,index_col="Gene", sep ="\t")
    dn_coreH = list(set([i+".fa.aln_adj_cds" for i in Core_hard]).intersection(set(list(dn.index))))
    dn_coreA = list(set([i+".fa.aln_adj_cds" for i in Core_all]).intersection(set(list(dn.index))))
    dn_speicify = list(set([i+".fa.aln_adj_cds" for i in StrinSpecific]).intersection(set(list(dn.index))))
    dn_coreS = list(set([i+".fa.aln_adj_cds" for i in Core_soft]).intersection(set(list(dn.index))))
    dn_shell = list(set([i+".fa.aln_adj_cds" for i in shell]).intersection(set(list(dn.index))))
    dn_cloud = list(set([i+".fa.aln_adj_cds" for i in cloud]).intersection(set(list(dn.index))))
    dn.loc[dn_coreH,"w"].to_csv(f"{prefix}_dnds_coreH.csv")
    dn.loc[dn_coreS,"w"].to_csv(f"{prefix}_dnds_coreS.csv")
    dn.loc[dn_shell,"w"].to_csv(f"{prefix}_dnds_shell.csv")
    dn.loc[dn_cloud,"w"].to_csv(f"{prefix}_dnds_cloud.csv")
    dn.loc[dn_coreA,"w"].to_csv(f"{prefix}_dnds_coreA.csv")
    dn.loc[dn_speicify,"w"].to_csv(f"{prefix}_dnds_strainSPE.csv")


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
if(pNpS):
    COREA = list(pn.loc[pn_coreA,"pN/pS_found"])
    STRAIN =list(pn.loc[pn_speicify,"pN/pS_found"])
    COREH = list(pn.loc[pn_coreH,"pN/pS_found"])
    CORES = list(pn.loc[pn_coreS,"pN/pS_found"])
    SHELL = list(pn.loc[pn_shell,"pN/pS_found"])
    CLOUD = list(pn.loc[pn_cloud,"pN/pS_found"])
    data = [COREH, CORES, SHELL, CLOUD, COREA, STRAIN]
    fig, ax = plt.subplots()
    y_major_locator=MultipleLocator(1)
    x_major_locator=MultipleLocator(1)
    ax.yaxis.set_major_locator(y_major_locator)
    ax.xaxis.set_major_locator(x_major_locator)
    plt.ylim(0,max(COREH+CORES+SHELL+CLOUD+COREA+STRAIN)+1)
    ax.violinplot(data,showmeans=True)
    ax.set_title("pn/ps of different gene sets")
    ax.set_xticklabels(["","COREH","CORES","SHELL","CLOUD","COREA","STRAIN"])
    plt.savefig(f"{prefix}_pn_ps_violin.pdf")
    plt.show()
    if stats.levene(COREA,STRAIN)[1] <= 0.05:
        equal_var = False
    else:
        equal_var = True
    with open(f"{prefix}_pnps_COREAcSTRAIN_test.tsv",'w') as TEST:
        pvalue = stats.ttest_ind(COREA,STRAIN, equal_var = equal_var)[1]
        if pvalue > 0.05:
            TEST.write(f"Mean value of COREA {mean(COREA)} is not signifcantly different from the SPECIFIED {mean(STRAIN)}, pvalue is {pvalue}\n")
            print(f"Mean value of COREA {mean(COREA)} is not signifcantly different from the SPECIFIED {mean(STRAIN)}, pvalue is {pvalue}\n")
        else:
            TEST.write(f"Mean value of COREA {mean(COREA)} is signifcantly different from the SPECIFIED {mean(STRAIN)}, pvalue is {pvalue}\n")
            print(f"Mean value of COREA {mean(COREA)} is signifcantly different from the SPECIFIED {mean(STRAIN)}, pvalue is {pvalue}\n")

if(dNdS):
    COREA = list(dn.loc[dn_coreA,"w"])
    STRAIN =list(dn.loc[dn_speicify,"w"])
    COREH = list(dn.loc[dn_coreH,"w"])
    CORES = list(dn.loc[dn_coreS,"w"])
    SHELL = list(dn.loc[dn_shell,"w"])
    CLOUD = list(dn.loc[dn_cloud,"w"])
    data = [COREH, CORES, SHELL, CLOUD, COREA, STRAIN]
    fig, ax = plt.subplots()
    y_major_locator=MultipleLocator(1)
    x_major_locator=MultipleLocator(1)
    ax.yaxis.set_major_locator(y_major_locator)
    ax.xaxis.set_major_locator(x_major_locator)
    plt.ylim(0,max(COREH+CORES+SHELL+CLOUD+COREA+STRAIN)+1)
    ax.violinplot(data,showmeans=True)
    ax.set_title("dn/ds of different gene sets")
    ax.set_xticklabels(["","COREH","CORES","SHELL","CLOUD","COREA","STRAIN"])
    plt.savefig(f"{prefix}_dn_ds_violin.pdf")
    plt.show()
    if stats.levene(COREA,STRAIN)[1] <= 0.05:
        equal_var = False
    else:
        equal_var = True
    with open(f"{prefix}_dnds_COREAcSTRAIN_test.tsv",'w') as TEST:
        pvalue = stats.ttest_ind(COREA,STRAIN, equal_var = equal_var)[1]
        if pvalue > 0.05:
            TEST.write(f"Mean value of COREA {mean(COREA)} is not signifcantly different from the SPECIFIED {mean(STRAIN)}, pvalue is {pvalue}\n")
            print(f"Mean value of COREA {mean(COREA)} is not signifcantly different from the SPECIFIED {mean(STRAIN)}, pvalue is {pvalue}\n")
        else:
            TEST.write(f"Mean value of COREA {mean(COREA)} is signifcantly different from the SPECIFIED {mean(STRAIN)}, pvalue is {pvalue}\n")
            print(f"Mean value of COREA {mean(COREA)} is signifcantly different from the SPECIFIED {mean(STRAIN)}, pvalue is {pvalue}\n")

with open(f"Core_hard.tsv", "w") as CH:
    CH.write("\n".join(Core_hard))

with open(f"Core_soft.tsv", "w") as CS: 
    CS.write("\n".join(Core_soft))

with open(f"shell.tsv", "w") as SH:
    SH.write("\n".join(shell))

with open(f"cloud.tsv", "w") as CL:
    CL.write("\n".join(cloud))
with open(f"Core_all.tsv", "w") as CA:
    CA.write("\n".join(Core_all))
if(specified):
    if stats.levene(COREA,STRAIN)[1] <= 0.05:
        equal_var = False
    else:
        equal_var = True
    with open(f"{prefix}_pnps_COREAcSTRAIN_test.tsv",'w') as TEST:
        pvalue = stats.ttest_ind(COREA,STRAIN, equal_var = equal_var)[1]
        if pvalue > 0.05:
            TEST.write(f"Mean value of COREA {mean(COREA)} is not signifcantly different from the SPECIFIED {mean(STRAIN)}, pvalue is {pvalue}\n")
            print(f"Mean value of COREA {mean(COREA)} is not signifcantly different from the SPECIFIED {mean(STRAIN)}, pvalue is {pvalue}\n")
        else:
            TEST.write(f"Mean value of COREA {mean(COREA)} is signifcantly different from the SPECIFIED {mean(STRAIN)}, pvalue is {pvalue}\n")
            print(f"Mean value of COREA {mean(COREA)} is signifcantly different from the SPECIFIED {mean(STRAIN)}, pvalue is {pvalue}\n")

