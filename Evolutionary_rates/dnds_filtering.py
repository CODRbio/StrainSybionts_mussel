#!/usr/bin/env python3
import sys, getopt, os
import pandas as pd



def usage():
    print('''
    -h help
    -m M0 results of all the species as criteria for filtering
    -t target dnds results to filter
    -o output prefix
    ''')

output = "filtered_RES.tsv"
opts,args=getopt.getopt(sys.argv[1:],'-h-t:-m:-o:',['help'])
for k, v in opts:
    if k in ("-h", "--help"):
        usage()
        sys.exit()
    elif k in ("-o"):
        output = v
    elif k in ("-t"):
        target = v
    elif k in ("-m"):
        m0 = v
    else:
        usage()
        sys.exit()

path = os.getcwd()

m0DF = pd.read_table(m0,sep="\t",header=0,index_col=0)
tarDF = pd.read_table(target,sep="\t",header=0)
m0DF["filtered"] = (m0DF["dN"]*m0DF["N"] >0.5) & (m0DF["dS"]*m0DF["S"] >0.5) & (m0DF["dS"]<2)
print(m0DF)
tarDF["filtered"] = tarDF.apply(lambda x:x["w"]<4 or ((x["w"] >=4 and x["w"]<10) and (m0DF.loc[x["Gene"],"filtered"])) or (x["w"]>=10 and x["dN"]*x["N"]>0.5 and x["dS"]*x["S"]>0.5),1)
res=tarDF[tarDF["filtered"]]
print(res)
res.to_csv(output,sep='\t', index=None)
