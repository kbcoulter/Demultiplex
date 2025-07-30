#!/usr/bin/env python

##### IMPORTS #####
import argparse
import numpy as np
import matplotlib.pyplot as plt
import gzip
import bioinfo # IMPORTS convert_phred (Set to 33 encoding)

##### COLLECT ARGS #####
def get_args():
     parser = argparse.ArgumentParser(description="A script to average the quality scores at each position for all reads and generate a per nucleotide mean distribution of quality scores for read1, read2, index1, and index2.")
     parser.add_argument("-f", "--filename", help="An input FASTQ file", required=True, type=str)
     parser.add_argument("-k", "--readsize", help="Size (length) of the reads", required=True, type=int)
     parser.add_argument("-o", "--outname", help="A name for the output histogram: <outname>.png", required=True, type=str)
     return parser.parse_args()

args = get_args()

##### SCRIPT BODY #####
qscores = np.zeros(args.readsize, dtype = float) # Will hold summed qscores then the mean after div.
num_records = 0

with open(args.filename, "rt") as fastq: # FOR TESTING FILES
#with gzip.open(args.filename, "rt") as fastq: # FOR REAL FILES
     for index_line, line_data in enumerate(fastq):
            if index_line % 4 == 3:
                num_records += 1
                for ascii_ind, ascii_char in enumerate(line_data.strip()):
                    qscores[ascii_ind] += bioinfo.convert_phred(ascii_char) # Sum vals within the qscores array

qscores = (qscores / num_records) # Div sum by num records and update the qscores arr

##### PLOTS #####
plt.bar(range(args.readsize), qscores, color = "green", edgecolor = "green")
plt.xlim(-1,args.readsize)
#plt.ylim(0, np.max(qscores))
plt.xlabel("Base Position in Read (0 Ind.)", fontsize = 14)
plt.ylabel("Average Quality Score", fontsize = 14)
plt.title(f"Mean Quality Score by {args.filename} Base Position", fontsize = 16, fontweight = "bold")
plt.grid(axis='y', linestyle='--', color = "grey", alpha=0.5)
plt.savefig(f"{args.outname}.png")