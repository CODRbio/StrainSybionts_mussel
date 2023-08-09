#!/nfs_genome/anaconda3/envs/rnaseq/bin/python
# -*- coding: utf-8 -*-
import getopt
import os
import re
import sys
from Bio import SeqIO
from pathos.multiprocessing import ProcessingPool as Pool
from pathos.pools import ProcessPool
import pandas as pd
import pysam
import itertools
import concurrent.futures
import time


def usage():
    print("calculate the reads counts for contig (RPK for single sam file) or MAG (RPKG for bams in folders)\n"
          " -h help\n" 
          " -b <folder> Folder containing binned fasta files (optional)\n" 
          " -c <coverage> reads length coveragem, default 0.7\n"
          " -f <folder> Folder containing bam files (optional)\n" 
          " -S indicate output only the count table for binning\n"
          " -i <identity> Identity threshold for mapping reads, default 0.98\n"
          " -s indicate the format of the file is sam (optional, -f the input sam file name should be given as arg)\n"
          " -m <protocol> Protocol to deal with multiple mapping reads, 0 for discarding, 1 for count all,2 for count in threshold and 3 for count in proportion\n" 
          " -o <output> RPKG results of the draft genomes, should be ends with _RPKG\n"
          " -M <mum> maximum reads in the seperated sam files, default 150000\n"
          " -r <readsLength for RPKG calculation>, 0 for all the reads sequenced, 1 for the reads that mapped to references which is suggested for high contaminated samples such as gut microbiome and symbionts\n"
          "")

def get_bam(path):
    return [os.path.join(path, f) for f in os.listdir(path) if f.endswith(".bam")]

def split_sam(filename):
    contig_len = {}
    with open(filename, 'r') as f:
        header = []
        file_count = 0
        read_count = 0
        last_read = ""
        outfile = None
        dir_name = os.path.dirname(filename)
        if dir_name == "":
            dir_name = "./"
        base_name = os.path.splitext(os.path.basename(filename))[0]
        for line in f:
            # keep header lines separate
            if line.startswith('@'):
                header.append(line)
                if line.startswith('@SQ'):
                    parts = line.split("\t")
                    contig_name = parts[1][3:]
                    contig_length = int(parts[2][3:])
                    contig_len[contig_name] = contig_length
            else:
                fields = line.strip().split()
                read = fields[0]
                # If we have seen max_reads or this line is a new read, and we have an open file
                if outfile and ((read != last_read and read_count >= max_reads) or read_count >= 10 * max_reads):
                    outfile.close()
                    outfile = None
                # If we do not have an open file, open next file
                if not outfile:
                    file_count += 1
                    outfile = open(f'{dir_name}/{base_name}_part_{file_count:04d}.sam', 'w')
                    #print(f'{dir_name}/{base_name}_part_{file_count:04d}.sam')
                    # write header to each new file
                    for header_line in header:
                        outfile.write(header_line)
                    read_count = 0
                outfile.write(line)
                if read != last_read:
                    read_count += 1
                last_read = read
        # close last file
        if outfile:
            outfile.close()
    return contig_len

def get_bins(folder, suffix):
    bin_dict = {}
    length_dict = {}
    length_contig = {}
    for file in os.listdir(folder):
        if file.endswith(suffix):
            basename = os.path.splitext(file)[0]
            length_dict[basename] = 0
            for seq_record in SeqIO.parse(os.path.join(folder, file), "fasta"):
                id = seq_record.id.split("|")[-1]
                length_contig[id] += len(seq_record)
                bin_dict[id] = basename
                length_dict[basename] += len(seq_record)
    return [bin_dict, length_dict, length_contig]

def get_score(bam):  # function to get the score of each read
    basename=os.path.splitext(os.path.basename(bam))[0]
    bamfile = pysam.AlignmentFile(bam, operation, threads=10, check_sq=False)
    score = {}
    best = {}
    query_reads = {}
    mapped_reads = {}
    for read in bamfile:
        readc = 1 if read.is_read1 else (2 if read.is_read2 else 0)
        if read.is_mapped:
            identity = 1 - read.get_tag("NM")/read.query_alignment_length
            #print(f"{read.qname}\t{read.query_alignment_length}\t{read.infer_read_length()}")
            coverage = int(read.query_alignment_length) / int(read.infer_read_length())
            #print(f"{read.query_name}\t{coverage}")
            if identity >= threshold and coverage >= cov_threshold:
                mark = "__".join([read.qname, read.reference_name, str(read.pos), str(readc)]) 
                score[mark] = read.get_tag("AS")
                if not read.is_secondary:
                    best["__".join([read.qname, str(readc)])] = mark
                    mapped_reads["__".join([read.query_name,str(readc)])] = read.query_alignment_length
        query_reads["__".join([read.query_name,str(readc)])] = read.infer_read_length()
    bamfile.close()
    mapped_reads = sum(mapped_reads.values())
    readsLength = sum([x for x in query_reads.values() if x is not None])
    print(mapped_reads)
    print(readsLength)
    return([basename, score, best, mapped_reads, readsLength])

def read_len_mean(bam):
    bamfile = pysam.AlignmentFile(bam, operation, threads=10, check_sq=False)
    rL = []
    for read in bamfile:
        length = read.infer_read_length()
        rL.append(length)
    rL = [x for x in rL if x is not None] 
    return sum(rL) / len(rL)


def update_parameter(score_args): 
    print("Updating parameters\n")
    score, best = score_args[1], score_args[2]
    highest_score = {}
    multiple_forRead = {}
    equally_best = {}
    for queryID in best:
        highest_score.update({queryID:score.get(best[queryID],300)})
    for alignID in score:
        readname, readc = alignID.split("__")[0], alignID.split("__")[3]
        mark = "__".join([readname,readc])
        if score[alignID] >= highest_score.get(mark, 300):
            equally_best[alignID] = True
            multiple_forRead[mark] = multiple_forRead.get(mark, 0) + 1
    return [multiple_forRead, equally_best]

def update_bam(bam, equally_best):
    basename=os.path.splitext(os.path.basename(bam))[0]
    path = os.path.dirname(os.path.abspath(bam))
    if path == "":
        path = "./"
    with pysam.AlignmentFile(bam, operation, threads=10, check_sq=False) as infile, pysam.AlignmentFile(f"{path}/filtered_{basename}.bam", "wb", template=infile) as outfile:
        for read in infile:
            if read.is_mapped:
                readc = 1 if read.is_read1 else (2 if read.is_read2 else 0)
                mark = "__".join([read.qname, read.reference_name, str(read.pos), str(readc)])
                if equally_best.get(mark,False):
                    outfile.write(read)

def contig_count(ParameterUpdate, multiple,quantile):
    Reads_in_NumContigs, equally_best = ParameterUpdate[0], ParameterUpdate[1]
    contig_counts = {}
    print("generating adjusted counts for each contig")
    for alignID in equally_best:
        readname = alignID.split("__")[0]
        refID = alignID.split("__")[1]
        readc = alignID.split("__")[3]
        mark = "__".join([readname, readc])
        count = 0
        if multiple == 0:
            if Reads_in_NumContigs[mark] > 1:
                count = 0
            else:
                count = 1
        elif multiple == 1:
            count = 1
        elif multiple == 2:
            Founds = Reads_in_NumContigs.get(mark,1)
            if Founds >= quantile[4]:
                #count = 1 / Reads_in_NumContigs.get(mark,1)
                count = 1/Founds
            elif Founds >= quantile[3]:
                #count = 2 / Reads_in_NumContigs.get(mark,1)
                count = 2/Founds
            elif Founds >= quantile[2]:
                #count = 4 / Reads_in_NumContigs.get(mark,1)
                count = 4/Founds
            elif Founds >= quantile[1]:
                #count = 8 / Reads_in_NumContigs.get(mark,1)
                count = 8/Founds
            else:
                count = 10/Founds
                #count = 10 / Reads_in_NumContigs.get(mark,1)
        elif multiple == 3:
            Founds = Reads_in_NumContigs.get(mark,1)
            count = 1 / Founds
        contig_counts[refID] = contig_counts.get(refID, 0) + count
    return contig_counts

import numpy as np

def calculate_quantiles(data):
    sorted_data = sorted(data) # 从小到大排序
    quantiles = [0, 0.2, 0.4, 0.6, 0.8, 1] # 分位数
    return [np.quantile(sorted_data, q) for q in quantiles]


def wrap_pocessing(bam):
    print("parsing scores")
    scoreArgs = get_score(bam)
    basename = scoreArgs[0]
    prefix = os.path.splitext(bam)[0]
    print(f"processing {basename}")
    updated_score = update_parameter(scoreArgs)
    reads_IN = updated_score[0]
    contig_countRes = []
    quantile = calculate_quantiles(reads_IN.values()) 
    for id in range(4):
        contig_countRes.append(contig_count(updated_score, id,quantile))
    print("updating sam to give bam")
    if formatIn == "sam":
        print(f"generating Coverages for {basename}")
        for id in range(4):
            df=pd.DataFrame(list(contig_countRes[id].items()), columns=['Contigs', 'Counts_adj'])
            df.to_csv(f"{prefix}_{id}.tsv",index=False,sep="\t")
        equally_best = updated_score[1]
        print("filtering bams")
        mapped_reads = scoreArgs[3]
        readsLength = scoreArgs[4]
        with open(f"{prefix}.reads", 'w') as outfile:
            outfile.write(str(mapped_reads))
        with open(f"{prefix}.rLen", 'w') as outfile:
            outfile.write(str(readsLength))
        if not simple:
            update_bam(bam, equally_best)
            os.remove(bam)
    else:
        print("pocessing bam")
        return [{basename:contig_countRes},{basename:scoreArgs[3]}, {basename:scoreArgs[4]}]

def RPKG_MAGs(contigCountAll):
    contig_countDict = {}
    fqmapped_dict = {}
    for ind in contigCountAll:
        contig_countDict.update(ind[0])
        if readsUse == 1:
            fqmapped_dict.update(ind[1])
        else:
            fqmapped_dict.update(ind[2])
    mergedDF = pd.DataFrame.from_dict(contig_countDict[list(contig_countDict.keys())[0]], orient='index',
                                          columns=[list(contig_countDict.keys())[0]])
    for sample in list(contig_countDict.keys())[1:]:
        DF = pd.DataFrame.from_dict(contig_countDict[sample], orient='index', columns=[sample])
        mergedDF = pd.merge(mergedDF, DF, left_index=True, right_index=True, how='outer')
    mergedDF = mergedDF.fillna(0)
    mergedDF = mergedDF.astype(int)
    suffix = ("fa", "fna", "fasta")
    [BinDict, lengthDict] = get_bins(bin_folder, suffix)
    mergedDF["MAG"] = mergedDF.index.map(lambda x: BinDict[x])
    RPK = mergedDF.groupby("MAG").sum()
    RPK = RPK.apply(lambda x: x / lengthDict[x.name] * 1e3, 1)
    RPKG = RPK.apply(lambda x: x / fqmapped_dict[x.name] * 1e9, axis=0)
    RPKG = RPKG.fillna(0)
    return RPKG

def RPKG_contigs(contigCountAll,output):
    fqmapped_dict = {}
    if readsUse == 1:
        fqmapped_dict = contigCountAll[1]
    else:
        fqmapped_dict = contigCountAll[2]
    contig_countDict = contigCountAll[0]
    basename = list(contig_countDict.keys())[0]
    mergedDF = pd.DataFrame.from_dict(contig_countDict[basename], orient='index',
                                      columns=[basename]) 
    RPK = mergedDF.apply(lambda x: x / contig_len[x.name] * 1e3, 1)
    RPKG = RPK.apply(lambda x: x / fqmapped_dict[x.name] * 1e9, axis=0)
    RPKG.to_csv(output, sep="\t")


def sort_bam(file):
    sorted_file = os.path.splitext(file)[0] + ".sorted.bam"
    sorted_file = re.sub(r'filtered_',"",sorted_file)
    pysam.sort("-o", sorted_file, file)
    return sorted_file

def sum_lengths(file_list):
    total_length = 0
    for file in file_list:
        with open(file, 'r') as f:
            for line in f:
                length = int(line.strip())  # assuming each line is a single integer
                total_length += length
    return total_length

def calculate_depth(row):
    depths = {}
    name = row["Contigs"]
    contig_length = contig_len[name]
    total_bp = row["Counts_adj"] * readL 
    depth = total_bp / contig_length
    return depth

def process_id(id_tsv_list):
    id, tsv_list = id_tsv_list  # 解包元组
    dfs = [pd.read_csv(table, header=0, sep="\t") for table in tsv_list]
    combined_df = pd.concat(dfs)
    summarized_df = combined_df.groupby('Contigs')['Counts_adj'].sum().reset_index()
    depths = summarized_df.apply(calculate_depth, axis=1)
    depth_df = pd.DataFrame({'Contigs': summarized_df['Contigs'], base_name: depths})
    summarized_df.to_csv(f"{output_dir}/{count_base}_{id}.tsv",index=False,header=True,sep="\t")
    depth_df.to_csv(f"{output_dir}/{depth_base}_{id}.tsv",index=False,header=True,sep="\t")
    return summarized_df

# Parsing command line arguments
if len(sys.argv) == 1:
    usage()
    sys.exit()
threshold = 0.98
cov_threshold = 0.7
max_reads = 150000
simple = False

opts,args=getopt.getopt(sys.argv[1:], 'hsSr:o:b:f:m:M:i:c:', ['help'])
for opt, arg in opts:
    if opt in ("-h", "--help"):
        usage()
        sys.exit()
    elif opt in ("-o"):
        output = arg
    elif opt in ("-s"):
        formatIn = "sam"
    elif opt in ("-b"):
        bin_folder = arg
    elif opt in ("-f"):
        bamIN = arg
    elif opt in ("-m"):
        multiple = int(arg)
    elif opt in ("-i"):
        threshold = float(arg)
    elif opt in ("-r"):
        readsUse = int(arg)
    elif opt in ("-c"):
        cov_threshold = float(arg)
    elif opt in ("-M"):
        max_reads = int(arg)
    elif opt in ("-S"):
        simple = True
    else:
        usage()
        sys.exit()

if formatIn == "sam":
    operation = "r" 
    dir_name = os.path.dirname(bamIN)
    if dir_name == "":
        dir_name = "./"
    output_dir = os.path.dirname(output)
    if output_dir == "":
        output_dir = "./"
    output_base = os.path.basename(output)
    depth_base = re.sub(r'_RPKG$',"_DEPTH",output_base)
    count_base = re.sub(r'_RPKG$',"_COUNT",output_base)
    contig_len = split_sam(bamIN) 
    time.sleep(1)
    if simple:
        try:
            os.remove(bamIN)
        except:
            pass
    base_name = os.path.splitext(os.path.basename(bamIN))[0]
    splitedSams = [os.path.join(dir_name, f) for f in os.listdir(dir_name) if f.startswith(f"{base_name}_part_")] 
    readL = read_len_mean(splitedSams[0])
    with Pool(ncpus=50) as pool:
        pool.map(wrap_pocessing, splitedSams)
    tsv_list = []
    for id in range(4):
        tsv_list.append([re.sub(r'\.sam$', f"_{id}.tsv", f) for f in splitedSams])
    bam_list = [re.sub(r'\.sam$', '.bam', f) for f in splitedSams]
   
    with Pool(ncpus=5) as pool:
        summarized_df = pool.map(process_id, enumerate(tsv_list))
    for tsv in tsv_list:
        for file in tsv:
            try:
                os.remove(file)
            except:
                print("No files\n")
    rLen_list = [re.sub(r'\.sam$', '.rLen', f) for f in splitedSams]
    readsList = [re.sub(r'\.sam$', '.reads', f) for f in splitedSams]
    readsT = sum_lengths(readsList)
    rLenT = sum_lengths(rLen_list)
    cwd =os.getcwd()
    for id in range(4):
        print(summarized_df[id])
        contig_countRes = summarized_df[id].set_index('Contigs')['Counts_adj'].to_dict()
        contigCountAll = [{base_name:contig_countRes},{base_name:readsT}, {base_name:rLenT}]
        RPKG_contigs(contigCountAll,f"{output_dir}/{output_base}_{id}.tsv")
    try:
        if os.path.isabs(output_dir):
            os.symlink(f"{output_dir}/{output_base}_{multiple}.tsv", f"{output_dir}/{output_base}.tsv")
            os.symlink(f"{output_dir}/{depth_base}_{multiple}.tsv",f"{output_dir}/{depth_base}.tsv")
        else:
            os.symlink(f"{cwd}/{output_dir}/{output_base}_{multiple}.tsv", f"{output_dir}/{output_base}.tsv") 
            os.symlink(f"{cwd}/{output_dir}/{depth_base}_{multiple}.tsv",f"{output_dir}/{depth_base}.tsv")
    except Exception as e:
        print(f"Failed to create symbolic link. Error: {str(e)}")

    for file in rLen_list:
        try:
            os.remove(file)
        except:
            pass
    for file in readsList:
        try:
            os.remove(file) 
        except:
            pass
    for file in splitedSams:
        try:
            os.remove(file)
        except:
            pass
    if not simple:
        filtered_list = []
        for file in bam_list:
            basename = os.path.splitext(os.path.basename(file))[0]
            basename =f"filtered_{basename}"
            dir=os.path.dirname(bamIN)
            if dir == "":
                dir = "./"
            path=f"{dir}/{basename}.bam"
            filtered_list.append(path)
        bam_list = [f"filtered_{f}" for f in bam_list] 
        with concurrent.futures.ThreadPoolExecutor(max_workers=50) as executor: 
            sorted_bam_files = list(executor.map(sort_bam, filtered_list))
        pysam.merge(f"{output_dir}/{base_name}.sorted.bam", *sorted_bam_files)
        for file in filtered_list:
            try:
                os.remove(file)
            except:
                print("No files\n")
        for file in sorted_bam_files:
            try:
                os.remove(file)
            except:
                print("No files\n")
        try:
            os.remove(bamIN)
        except:
            print(f"The {bamIN} cannot be deleted\n")
else:
    operation = "rb"
    bam_list = get_bam(bamIN)
    with Pool(ncpus=10) as pool:
        contigCountAll = pool.map(wrap_pocessing, bam_list)
    RPKG_MAGs = RPKG_MAGs(contigCountAll)
    RPKG_MAGs.to_csv(output, sep="\t")
