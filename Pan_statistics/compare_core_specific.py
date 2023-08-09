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


def usage():
    print('''
    use to compare the selection pressures in core genes and strain_specific genes
    -p pn_ps result (optional) 
    -d filtered dn_ds result files (optional)
    -s strain-specific ortholog list
    -c core ortholog list
    -o output prefix
    ''')
pNpS = False
dNdS = False
opts,args=getopt.getopt(sys.argv[1:],'-h-p:-d:-c:-o:-s:',['help'])
for k, v in opts:
    if k in ("-h", "--help"):
        usage()
        sys.exit()
    elif k in ("-o"):
        outputP = v
    elif k in ("-p"):
        pNpS = v
    elif k in ("-d"):
        dNdS = v
    elif k in ("-s"):
        specific = v
    elif k in ("-c"):
        core = v
    else:
        usage()
        sys.exit()

if len(opts)>0:
    pass
else:
    usage()
    sys.exit()
if(pNpS):
    df = pd.read_csv(pNpS,sep="\t",header=0)
elif(dNdS):
    df = pd.read_csv(dNdS,sep="\t",header=0)
specificList = pd.read_csv(specific,header=None).iloc[:,0].tolist()
coreList = pd.read_csv(core,header=None).iloc[:,0].tolist()
coreDF = df.loc[df.iloc[:,0].isin(coreList)]
specificDF = df.loc[df.iloc[:,0].isin(specificList)]
print(coreDF)
if(pNpS):
    CORE = list(coreDF["pN/pS_found"])
    SPECI = list(specificDF["pN/pS_found"])
elif(dNdS):
    CORE = list(coreDF["w"])
    SPECI = list(specificDF["w"])
pd.DataFrame(CORE, columns=['Core']).to_csv(f"{outputP}_core.tsv",index=False)
pd.DataFrame(SPECI, columns=['SPECI']).to_csv(f"{outputP}_strainSpecific.tsv",index=False)
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

plt.savefig(f"{outputP}_violin.pdf")
plt.show()


if stats.levene(CORE,SPECI)[1] <= 0.05:
    equal_var = False
else:
    equal_var = True
with open(f"{outputP}_test.tsv",'w') as TEST:
    pvalue = stats.ttest_ind(CORE,SPECI, equal_var = equal_var)[1]
    if pvalue > 0.05:
        TEST.write(f"Mean value of CORE {mean(CORE)} is not signifcantly different from the SPECIFIED {mean(SPECI)}, pvalue is {pvalue}\n")
        print(f"Mean value of CORE {mean(CORE)} is not signifcantly different from the SPECIFIED {mean(SPECI)}, pvalue is {pvalue}\n")
    else:
        TEST.write(f"Mean value of CORE {mean(CORE)} is signifcantly different from the SPECIFIED {mean(SPECI)}, pvalue is {pvalue}\n")
        print(f"Mean value of CORE {mean(CORE)} is signifcantly different from the SPECIFIED {mean(SPECI)}, pvalue is {pvalue}\n") 
