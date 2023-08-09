#!/nfs_genome/anaconda/envs/rnaseq/bin/python
import os
import sys, getopt
import fnmatch
import pandas as pd
import pathos

def usage():
    print('''
    use to caculate statistic assessment of the reads binning results, including the checm of each binned reads and the proportion of strains they are derived from
    -d position of the mergedDEPTH.csv
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


opts,args=getopt.getopt(sys.argv[1:],'-h-d:',['help'])
for k, v in opts:
    if k in ("-h", "--help"):
        usage()
        sys.exit()
    elif k in ("-o"):
        output = v
    elif k in ("-d"):
        depth = v
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

def execCMD(cmd):
    try:
        print(f"Now running cmd {cmd}")
        os.system(cmd)
    except:
        print(f"Fail: {cmd}")


Dir = os.path.dirname(depth)
sample = Dir.split(os.sep)[-2]


cluster = [10,50,100,400,800]
commandList = []

for cvalue in cluster:
    try:
        os.mkdir(f"{sample}/c{cvalue}")
    except:
        print("The directory exists!")
    if os.path.isfile(f"{sample}/c{cvalue}/clustering_gt1000.csv"):
        continue
    command = f"/nfs_genome/anaconda/envs/concoct-env/bin/concoct -c {cvalue} -k 4 --threads 60 -b {sample}/c{cvalue} -l 1000 -i 2000 --composition_file {sample}/work_files/assembly.fa --coverage_file {sample}/work_files/mergedDEPTH.tsv"
    commandList.append(command)
p = pathos.multiprocessing.ProcessPool(nodes=10)
p.map(execCMD,commandList)
p.close()
p.join()

pattern="clustering_gt1000.csv"

def find_files(directory, pattern):
    for root, dirs, files in os.walk(directory):
        for name in files:
            if name.endswith(pattern):
                yield os.path.join(root, name)

current_directory = os.getcwd()
concoctRes = list(find_files(f"./{sample}",pattern))
checkmCom =[]
import glob


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
print(checkmCom)

p = pathos.multiprocessing.ProcessPool(nodes=12)
p.map(execCMD, checkmCom)
p.close()
p.join()

