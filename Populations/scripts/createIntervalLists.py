#!/usr/bin/env python3

# Usage:
# ./createIntervalLists.py seq.fa.gz out.dir min.length interval.number

import re
import gzip
import sys
import os

# The initial order of the input sequences has to be preserved. Python starting
# with CPython 3.6 guarantees that dictionary order is identical to insertion
# order. This feature is used here. So the Python version has to be checked.
# Not preserving the order causes trouble with gatk GatherVcfs.
PythonMajor = sys.version_info.major
PythonMinor = sys.version_info.minor

if PythonMajor < 3 and PythonMinor < 6:
    print("Python version has to be at least 3.6. Found the following version:",
          str(PythonMajor) + "." + str(PythonMinor))
    print("Program EXIT.")
    sys.exit()

# check that the number of supplied arguments is sufficient
if len(sys.argv) != 5:
    print("\nThe number of supplied aguments is not sufficient.",
          "In additon to the program name, four arguments must begiven:",
          "  1. fasta sequence name",
          "  2. name of output directory",
          "  3. minimum length of contigs to be considered",
          "  4. number of intervals to create",
          "Supplied number of arguments:" + str( len(sys.argv) ) + "\n",
          sep="\n", file=sys.stderr)
    print("The following arguments were supplied:")
    for n,m in enumerate(sys.argv):
        if n == 0:
            print("  Script name:", m, file=sys.stderr)
        else:
            print("  ", n, ": ", m, sep='', file=sys.stderr)
    sys.exit()

DEBUGGING = False #True

# fasta file to be analyzed
FASTA_SEQ = sys.argv[1]
# output directory
OUT_DIR = sys.argv[2]
# min length for sequence to be included
TO_SHORT_SEQ = int(sys.argv[3])
# number of intervals to create
INTERVAL_NUM = int(sys.argv[4])

# does the output directory exist?
if not os.path.isdir(OUT_DIR):
    print("Output directory is not existent.\n"
          "Program EXIT.")
    sys.exit()

#************************************************
# store all sequences and their length in a dict
#************************************************
seqDict = {}
suffix =  os.path.splitext(FASTA_SEQ)[-1][1:]
if (suffix in ["gz","gzip","GZ","GZIP"]):
    with gzip.open(FASTA_SEQ, 'r') as fasta:
        seqLen = 0
        regex = re.compile('^>')
        for line in fasta:
            line = line.decode('utf-8')
            line = line.rstrip()
            if regex.match(line):
                if seqLen:
                    seqDict[key] = seqLen
                    seqLen = 0
                key = line.split(' ')[0].lstrip('>')
            else:
                seqLen += len(line)
        else:
            # print the length of the last sequence
            # after no line is left
            seqDict[key] = seqLen
else:
    with open(FASTA_SEQ, 'r') as fasta:
        seqLen = 0
        regex = re.compile('^>')
        for line in fasta:
            #line = line.decode('utf-8')
            line = line.rstrip()
            if regex.match(line):
                if seqLen:
                    seqDict[key] = seqLen
                    seqLen = 0
                key = line.split(' ')[0].lstrip('>')
            else:
                seqLen += len(line)
        else:
                # print the length of the last sequence
                # after no line is left
            seqDict[key] = seqLen

if DEBUGGING:
    for key in seqDict:
        print(key, seqDict[key], sep='\t')

totalGenomeLength = sum(seqDict.values())

DEBUGGING and print("total Length:", totalGenomeLength)

#***********************
# remove short sequences
#***********************
longSeqDict = {}
for key in seqDict:
    if seqDict[key] > TO_SHORT_SEQ:
        longSeqDict[key] = seqDict[key]

if DEBUGGING:
    print("number of remaining sequences:", len(longSeqDict))

finalGenomeLength = sum(longSeqDict.values())

#*********************************
# create a dict for the intervals
#*********************************
intervalLength = finalGenomeLength / INTERVAL_NUM

intervalDict = {}
intervalCounter = 1
intervalDict = {}
#intervalDict["interval_"+str(intervalCounter)] = []
currentInterval = "interval_" + str(intervalCounter)
# add the first interval to the intervalDict
intervalDict[ currentInterval ] = []
cumulativeLen = 0
# Python >= 3.6 uses ordered dicts; initial order is preserved
for key in longSeqDict.keys():
    DEBUGGING and print("\t",cumulativeLen)
    # the first contig has to be added to the interval in any case
    # the length of the interval doesn't matter
    if len( intervalDict[ currentInterval ] ) == 0:
        intervalDict[ currentInterval ].append(key)
        cumulativeLen += longSeqDict[key]
        continue
    # add the current sequence if space is left and
    # make sure that the added contig won't exceed the max interval length by more than 5 %
    if cumulativeLen < intervalLength and \
       (cumulativeLen + longSeqDict[key]) < ( intervalLength * 1.05 ):
        intervalDict[ currentInterval ].append(key)
        cumulativeLen += longSeqDict[key]
    # if no space left: create new entry in intervalDict
    else:
        intervalCounter += 1
        currentInterval = "interval_" + str(intervalCounter)
        cumulativeLen = longSeqDict[key]
        # add the new interval to the intervalDict
        intervalDict[ currentInterval ] = []
        # add the current contig to the interval
        intervalDict[ currentInterval ].append(key)

if DEBUGGING:
    for key in intervalDict:
        print(key, ":", intervalDict[key])
    print( intervalDict.keys() )

#**************************************
# control: all sequences still present?
#**************************************

# count entries in intervalDict
intervals = 0
for key in intervalDict:
    intervals += len(intervalDict[key])

if len(longSeqDict) != intervals:
    print("Potential data corruption.\n"
          "Number of intervals has changed.\n"
          "Program exit!")
    sys.exit()

#*****************************
# write the intervals to files
#*****************************
for key in intervalDict:
    with open(OUT_DIR +'/'+ key+'.list', 'w') as interval:
        for entry in intervalDict[key]:
            interval.writelines(entry + '\n')

if DEBUGGING:
    for key in intervalDict:
        print(key)
        for entry in intervalDict[key]:
            print("\t", entry)

