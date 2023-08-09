#!/nfs_genome/anaconda/envs/rnaseq/bin/python
import os
import sys, getopt
import fnmatch
import pandas as pd
import pathos
import re

def usage():
    print('''
    use to assemble the selected read bins into MAG
    -o  output path, default "./"
    -c  paths of the checkm tables
    -p  root paths of the folders containing samples
    -C  completeness threshold 95
    -D  comtamination threshold 1000
    -H  Heterogeneity threhold, default 80
    ''')

pattern="_checkm.tsv"

def find_files(directory, pattern):
    for root, dirs, files in os.walk(directory):
        for name in files:
            if name.endswith(pattern):
                yield os.path.join(root, name)

def execCMD(cmd):
    try:
        print(f"Now running cmd {cmd}")
        os.system(cmd)
    except:
        print(f"Fail: {cmd}")
completeness = 95
output = "./"
comtamination = 1000
hetero = 80

opts,args=getopt.getopt(sys.argv[1:],'-h-c:p:C:H:D:o:',['help'])
for k, v in opts:
    if k in ("-h", "--help"):
        usage()
        sys.exit()
    elif k in ("-c"):
        checkm = v
    elif k in ("-p"):
        rPath = v
    elif k in ("-C"):
        completeness = float(v)
    elif k in ("-o"):
        output =v
    elif k in ("-H"):
        hetero = float(v)
    elif k in ("-D"):
        comtamination = float(v)
    else:
        usage()
        sys.exit()

if len(opts)>0:
    pass
else: 
    usage()
    sys.exit()
commandLists = []
contigCheck= []
checkmList = find_files(checkm,pattern)
checkmList = [os.path.basename(x) for x in checkmList]
for checkm in checkmList:
    match = re.match(r"(.+)_(c\d+)_checkm\.tsv", checkm)
    if match:
        sample = match.group(1)  # 提取样品部分
        cluster = match.group(2)  # 提取聚类策略部分
        try:
            os.mkdir(f"{sample}_{cluster}")
        except:
            pass
        df = pd.read_csv(checkm,sep=r'\s+',header=None,skiprows=[0,1,2], skipinitialspace=True)
        
        selected_names = df.loc[(df.iloc[:,-1] > hetero) & (df.iloc[:, -2] > comtamination) & (df.iloc[:, -3] > completeness), 0]
        awk_command = "/^S/{print \">\"$2;print $3}"
        for bins in selected_names:
            command = f"/nfs_genome/anaconda/envs/pb_tools/bin/flye --min-overlap 8000 -g 5m -t 80 -i 2 --pacbio-hifi {rPath}/{sample}/{cluster}/{bins}.fasta -o {output}/{sample}_{cluster}_{bins}"
            #command = f"/nfs_genome/anaconda/envs/pb_tools/bin/hifiasm_meta -t 90 -o {output}/{sample}_{cluster}_{bins}_asm {rPath}/{sample}/{cluster}/{bins}.fasta; awk '{awk_command}' {sample}_{cluster}_{bins}_asm.p_ctg.gfa > {output}/{sample}_{cluster}_{bins}_asm.p_ctg.fa"
            commandLists.append(command)
            contigCheck.append(f"checkm taxonomy_wf family Enterobacteriaceae -x fasta -t 80 {output}/{sample}_{cluster}_{bins} {output}/C_{sample}_{cluster}_{bins}_checkm >{output}/C_{sample}_{cluster}_{bins}_checkm.tsv")
    else:
        print("Filename does not match the expected pattern.")

p = pathos.multiprocessing.ProcessPool(nodes=8)
p.map(execCMD,commandLists)
p.close()
p.join()
#execCMD(f"checkm taxonomy_wf family Enterobacteriaceae -x p_ctg.fa -t 90 {output} {output}/checkm >{output}/checkm.tsv")
p = pathos.multiprocessing.ProcessPool(nodes=10)
p.map(execCMD,contigCheck)
p.close()
p.join()

