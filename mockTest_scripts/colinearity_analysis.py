#!/nfs_genome/anaconda/envs/rnaseq/bin/python
import os
import pathos
import pandas as pd
import sys
import getopt

current_dir = os.getcwd()

def nucmer_analysis(row):
    sample = row['Query'].split('_genomic', 1)[0]
    sample_name = os.path.basename(sample)
    ref_name = os.path.basename(row['Reference'])
    command = f"nucmer -c 15000 -p {sample_name}_{ref_name} {row['Query']} {row['Reference']}\nshow-tiling -l 1000 {sample_name}_{ref_name}.delta  >{sample_name}_{ref_name}_tiling.res\nmummerplot {sample_name}_{ref_name}.delta -f --fat --color -p {sample_name}_{ref_name} --postscript"
    return command

def main(argv):
    input_file = ''
    try:
        opts, args = getopt.getopt(argv,"hi:s:",["ifile="])
    except getopt.GetoptError:
        print('filter.py -i <inputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('filter.py -i <inputfile>\n -s <similarity>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            input_file = arg
        elif opt in ("-s"):
            similarity = float(arg)
    return input_file, similarity

file_name,similarity = main(sys.argv[1:])

df = pd.read_csv(file_name, sep="\t", header=None, names=['Query', 'Reference', 'ANI', 'Mapped', 'TotalMarkers'])

df['Query_is_ref'] = df['Query'].str.contains('genomic')
df['Ref_is_emulated'] = df['Reference'].str.contains('merge')

df_filtered = df[(df['Query_is_ref'] & df['Ref_is_emulated'] & (df['ANI'] > similarity))]

commandList = df_filtered.apply(nucmer_analysis, axis=1)

def execCMD(cmd):
    try:
        print(f"Now running cmd {cmd}")
        os.system(cmd)
    except:
        print(f"Fail: {cmd}")

p = pathos.multiprocessing.ProcessPool(nodes=12)
p.map(execCMD, commandList)
p.close()
p.join()

