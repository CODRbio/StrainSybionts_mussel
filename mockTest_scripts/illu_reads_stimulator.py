#!/nfs_genome/anaconda/envs/rnaseq/bin/python
#work in sberry environment
import pandas as pd
import sys, getopt, os
import pathos


def usage():
    print('''
    use to statistic contig nums of each cluster
    -t tables providing coverage information of each strain
    -d directory containing the ref genomes
    -o output folds
    -n total depth to generate, default 2000
    -g average genome size, default 4.5m
    ''')

opts,args=getopt.getopt(sys.argv[1:],'-h-d:-t:-o:-g:-n:',['help'])
output = "./"
genome = 4.5
depth = 2000
for k, v in opts:
    if k in ("-h", "--help"):
        usage()
        sys.exit()
    elif k in ("-t"):
        design = v
    elif k in ("-o"):
        output = v
    elif k in ("-d"):
        directory = v
    elif k in ("-n"):
        depth = int(v)
    elif k in ("-g"):
        genome = float(v)
    else:
        usage()
        sys.exit()
if len(opts)>0:
    pass
else:
    usage()
    sys.exit()

df =  pd.read_csv(design,header=0).fillna(0)
df = (df*depth*genome*1000000/100/300).round().astype(int)
genomeName = list(df.columns.values)

def execCMD(cmd):
    try:
        print(f"Now running cmd {cmd}") 
        os.system(cmd)
    except:
        print(f"Fail: {cmd}")

def command_gen(row):
    command_list = []
    #row=unlist(row)
    sample=row.name
    for cnt in range(len(row)):
        if row[cnt] > 0:
            command = f"wgsim -e 0.005 -R 0.05 -d 320 -1 150 -2 150 -N {row[cnt]} {directory}/{genomeName[cnt]} {output}/St{sample}_{cnt}.fq1 {output}/St{sample}_{cnt}.fq2\n gzip {output}/St{sample}_{cnt}.fq1&\n gzip {output}/St{sample}_{cnt}.fq2 &\n wait"
            command_list.append(command)
        else:
            print("Error with the input")
    p = pathos.multiprocessing.ProcessPool(nodes=6)
    p.map(execCMD,command_list)
    os.system(f"cat {output}/St{sample}_*.fq1.gz >{output}/mergedIND_{sample}.fq1.gz")
    os.system(f"cat {output}/St{sample}_*.fq2.gz >{output}/mergedIND_{sample}.fq2.gz")
    os.system(f"rm -f {output}/St{sample}_*.fq*.gz")

for idx,row in df.iterrows():
    command_gen(row)
    




