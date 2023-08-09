#!/usr/bin/env python
import sys, getopt
import pandas as pd
from collections import Counter

def usage():
    print('''
    use to summarize the distribution pattern of the COG categories 
    -c blast2cog annotation file
    -s selected orthologs in a file, none means all the seqs were used 
    -o output prefix
    ''')
selected = False


opts,args=getopt.getopt(sys.argv[1:],'-h-c:-o:-s:',['help'])
for k, v in opts:
    if k in ("-h", "--help"):
        usage()
        sys.exit()
    elif k in ("-o"):
        outputP = v
    elif k in ("-c"):
        cog = v
    elif k in ("-s"):
        selected = v
    else:
        usage()
        sys.exit()

if len(opts)>0:
    pass
else:
    usage()
    sys.exit()

df = pd.read_csv(cog, sep='\t',header=0)
if selected:
    selList = pd.read_csv(selected,header=None).iloc[:,0].tolist()
    df_sel = df.loc[df.iloc[:,0].isin(selList)]
else:
    df_sel = df

# Break down categories into individual letters and count each category
category_counts = Counter(''.join(df_sel['COG_catagory']))

# Convert the Counter object to a DataFrame
count_df = pd.DataFrame.from_dict(category_counts, orient='index', columns=['Count'])

print(count_df)
count_df.to_csv(f"{outputP}_COGstat.tsv",sep="\t")
