#!/nfs_genome/anaconda/envs/rnaseq/bin/python
import os
import sys, getopt
import fnmatch
import pandas as pd
import pathos

def usage():
    print('''
    use to carry out binning reads and statistic assessment of the binning results
    -p folders containing the mock pacbio reads
    -i folders containing the mock illumina reads
    -o output prefix
    -P threads in Prafly,default = 2
    -m how to deal with multiple hitting, 0 zero, 1 all, 2 thresholding, 3 proportion
    -t thread used in parrelization, default = 1
    -T temporary folders, default /BioData/tmp
    ''')

from Bio import SeqIO
def delete_empty_files(directory):
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        if os.path.isfile(file_path) and os.path.getsize(file_path) == 0:
            os.remove(file_path)
            print(f'{file_path} removed.')


def split_fasta_by_binning(input_fasta, binning_csv, column_number, output_dir, min_size):
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read binning information
    binning = pd.read_csv(binning_csv, sep=',', header=None)
    bin_dict = binning.set_index(0)[column_number - 1].to_dict()

    # Initialize or clear output files and bin sizes
    output_files = {bin_name: open(os.path.join(output_dir, f"bin.{bin_name}.fasta"), 'w') for bin_name in bin_dict.values()}
    output_files['unbinned'] = open(os.path.join(output_dir, "unbinned.fasta"), 'w')
    bin_sizes = {bin_name: 0 for bin_name in bin_dict.values()}
    bin_sizes['unbinned'] = 0

    # Parse fasta file
    with open(input_fasta, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            bin_name = bin_dict.get(record.id, 'unbinned')
            bin_sizes[bin_name] += len(record)
            if bin_sizes[bin_name] >= min_size:
                SeqIO.write([record], output_files[bin_name], 'fasta')

    # Close output files
    for f in output_files.values():
        f.close()
multiple = 1
threads = 1
Temporary = "/BioData/tmp"
parallization = 2
opts,args=getopt.getopt(sys.argv[1:],'-h-p:P:m:-i:-o:-t:-T:',['help'])
for k, v in opts:
    if k in ("-h", "--help"):
        usage()
        sys.exit()
    elif k in ("-o"):
        output = v
    elif k in ("-p"):
        pacbio = v
    elif k in ("-m"):
        multiple = int(v)
    elif k in ("-P"):
        parallization = int(v)
    elif k in ("-i"):
        illumina = v
    elif k in ("-s"):
        specific = v
    elif k in ("-t"):
        threads = int(v)
    elif k in ("-T"):
        Temporary = v
    else:
        usage()
        sys.exit()

if len(opts)>0:
    pass
else:
    usage()
    sys.exit()

def list_gz_files(directory):
    gz_files = [f for f in os.listdir(directory) if fnmatch.fnmatch(f, '*.gz')]
    return gz_files

def list_fna_files(directory):
    fa_files = [f for f in os.listdir(directory) if fnmatch.fnmatch(f, '*.fna')]
    return fa_files


illureads = list_gz_files(illumina)
illureads = [os.path.join(illumina, f) for f in illureads]
illureads = " ".join(illureads)
pacnames = list_fna_files(pacbio)

def execCMD(cmd):
    try:
        print(f"Now running cmd {cmd}")
        os.system(cmd)
    except:
        print(f"Fail: {cmd}")

def is_file_nonempty(filepath):
    # 检查文件是否存在
    if not os.path.isfile(filepath):
        return False
    # 检查文件是否非空
    if os.stat(filepath).st_size > 0:
        return True
    return False

def align_illu(pacname):
    try:
        os.mkdir(pacname)
    except:
        print("The directory exists!")
    command = f"metawrap binning -T {Temporary} -P {parallization} --mapOnly --strain -t 180 -M {multiple} -o {pacname} -a {os.path.join(pacbio,pacname)} {illureads}"
    if not(is_file_nonempty(f"{pacname}/work_files/mergedDEPTH.tsv")):
        return command

command_align_illu = list(map(align_illu,pacnames))
command_align_illu = [x for x in command_align_illu if x is not None]
p = pathos.multiprocessing.ProcessPool(nodes=threads)
p.map(execCMD,command_align_illu)
p.close()
p.join()

cluster = [6,8,10,50,100,400]

commandList = []
for sample in pacnames:
    for cvalue in cluster:
        try:
            os.mkdir(f"{sample}/c{cvalue}")
        except:
            print("The directory exists!")
        if os.path.isfile(f"{sample}/c{cvalue}/clustering_gt1000.csv"):
            continue
        command = f"/nfs_genome/anaconda/envs/concoct-env/bin/concoct -c {cvalue} -k 4 --threads 60 -b {sample}/c{cvalue} -l 1000 -i 2000 --composition_file {sample}/work_files/assembly.fa --coverage_file {sample}/work_files/mergedDEPTH.tsv"
        commandList.append(command)
p = pathos.multiprocessing.ProcessPool(nodes=threads*3)
p.map(execCMD,commandList)
p.close()
p.join()

pattern="clustering_gt1000.csv"

def find_files(directory, pattern):
    for root, dirs, files in os.walk(directory):
        for name in files:
            if name.endswith(pattern):
                yield os.path.join(root, name)

import glob
current_directory = os.getcwd()
checkmCom = []
for sample in pacnames:
    concoctRes = list(find_files(f"./{sample}",pattern))

    for item in concoctRes:
        Dir = os.path.dirname(item)
        clu = Dir.split(os.sep)[-1]
        if int(clu.replace("c", "")) not in cluster:
            continue
        sample = Dir.split(os.sep)[-2]
        command = f"clustering_stat_conococt.py -i {item} -m 1000 -o {sample}_{clu}_clustering.tsv"
        execCMD(command)
        file_list = glob.glob(f"{sample}/{clu}/*.fasta")
        if file_list:
            print("bins already parsed")
        else:
            split_fasta_by_binning(f"{current_directory}/{Dir}/../work_files/assembly.fa", f"{current_directory}/{item}",2,f"{sample}/{clu}/",10000000)
            delete_empty_files(f"{sample}/{clu}/")
        if os.path.exists(f"{sample}_{clu}_checkm.tsv"):
            pass
        else:
            command = f"checkm taxonomy_wf family Enterobacteriaceae -x fasta -t 80 {sample}/{clu} {sample}_{clu}_checkm >{sample}_{clu}_checkm.tsv"
            checkmCom.append(command)
p = pathos.multiprocessing.ProcessPool(nodes=12)
p.map(execCMD, checkmCom)
p.close()
p.join()

