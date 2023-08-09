#!/nfs_genome/anaconda/envs/rnaseq/bin/python
#work in sberry environment
import pandas as pd
import sys, getopt, os
import pathos


def usage():
    print('''
    use to stimulate hifi pacbio reads from the selected genomes
    -t table providing relative coverage (depth) information of each strain
    -d directory containing the ref genomes
    -o output folds
    -n totol depth coverage to generate (1,000 default)
    ''')

opts,args=getopt.getopt(sys.argv[1:],'-h-d:-t:-n:-o:',['help'])
output = "./"

num = 1000
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
        num = int(v)
    else:
        usage()
        sys.exit()
if len(opts)>0:
    pass
else:
    usage()
    sys.exit()

df =  pd.read_csv(design,header=0).fillna(0)
df = (df / 100 * num).round().astype(int)
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
            prefix = genomeName[cnt][0:3]
            command = f"pbsim --prefix {output}/I{sample}_{prefix} --id-prefix {prefix} --strategy wgs --genome {directory}/{genomeName[cnt]} --depth {row[cnt]} --length-min 500 --length-mean 14500 --pass-num 9 --errhmm /nfs_genome/BioInformatics/pbsim3-master/data/ERRHMM-SEQUEL.model --method errhmm\n"
            command = command + f"samtools view -b -@ 80 -o {output}/I{sample}_{prefix}.bam {output}/I{sample}_{prefix}_0001.sam\n"
            command = command + f"ccs --min-passes 7 {output}/I{sample}_{prefix}.bam {output}/I{sample}_{prefix}.fastq.gz"
            command_list.append(command)
        else:
            print("Error with the input, or the zero means the strain do not occur")
    p = pathos.multiprocessing.ProcessPool(nodes=10)
    p.map(execCMD,command_list)
    os.system(f"cat {output}/I{sample}*.fastq.gz >{output}/merged_{sample}.fastq.gz") 
    os.system(f"rm -f {output}/I{sample}_*")

for idx,row in df.iterrows():
    command_gen(row)
    




