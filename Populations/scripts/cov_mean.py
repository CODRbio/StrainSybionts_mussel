#!/usr/bin env python3
import os
import statistics

def do_avg_calculation(input,output):
    lists=[]
    for files in input:
        with open(files, 'r') as f:
            line = f.readline()
            item = float(line.strip().split()[-1])
            lists.append(item)
    with open(output, 'w') as f:
        f.write(str(statistics.mean(lists)))
    # python code
    print(statistics.mean(lists))
do_avg_calculation(snakemake.input, snakemake.output[0])
