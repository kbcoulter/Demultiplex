#!/usr/bin/env python

##### IMPORTS #####
import numpy as np 
import argparse

def reverse_comp(sequence: str) -> str:
    '''Takes a sequence (string) with or without newline characters and returns the reverse compliment of the string without newline characters in uppercase'''
    sequence = sequence.upper().strip()
    map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',}
    rev_comp = ("".join([map.get(entry, entry)for entry in sequence]))[::-1] # BRACKET REVERSES, OTHER COMPS
    return rev_comp

# Input: CATG
# Expected Output: GTAC -> Add Tests, or is this okay? 

def convert_phred(letter: str) -> int:
    '''Takes a single ASCII character (string) encoded in Phred+33 and returns the quality score value as an integer.'''
    return (ord(letter) - 33)

#Input: I
#Expected output: 40 -> Add Tests, or is this okay? 

def mean_qual_score(sequence: str) -> float:
    '''Takes a quality score sequence (string), uses convert_phred function to convert each phred+33 encoding into a integer quality score, then averages the quality score and returns the mean values. '''
    holder = []
    for char in sequence:
        holder.append(convert_phred(char))
    return (np.mean(holder))

def collect_record(file):
    record = []
    for line_counter in range(4):
        line = file.readline()
        if line:
            record.append(line.strip())

    return record


#Input: IIIIIIIIII
#Expected output: 40

# !!! ADD ARGPARSE AFTER TESTING !!!
"""
def get_args():
     parser = argparse.ArgumentParser(description="A script to average the quality scores at each position for all reads and generate a per nucleotide mean distribution of quality scores for read1, read2, index1, and index2.")
     parser.add_argument("-r1", "--Read1", help="An input FASTQ Read 1 file", required=True, type=str)
     parser.add_argument("-k", "--readsize", help="Size (length) of the reads", required=True, type=int)
     parser.add_argument("-o", "--outname", help="A name for the output histogram: <outname>.png", required=True, type=str)
     return parser.parse_args()

args = get_args()
"""

##### PSEUDOCODE OUTLINE FOR THE SCRIPT #####

# TO BE SET BY ARGPARSE LATER:
quality_score_cutoff = 5 #Create a quality_score_cutoff integer where every index/barcode read with quality score below goes into unknown. 
index_file = "/projects/bgmp/shared/2017_sequencing/indexes.txt"

### REAL DATA ###
#r1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
#r2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
#r3 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
#r4 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"

### TESTING DATA ### 
r1 = "/TEST-input_FASTQ/..."



##### SET GLOBAL VARS #####
match_dict = {} #Create a match_dict dictionary with the barcode or indices as keys, and match numbers initilizaed to 1 upon thier addition.
hopping_counter = 0 #Create a hopping_counter initialized to 0.
hopping_dict = {} #Create a hopping_dict dicitonary with nothing in it. The keys will be our Header_additions and the values will be counters initilaized to 1 upon their addition to the list. 
unknown_counter = 0 #Create a unkown_counter initialized to 0. 

# Lets create a file bank with barcodes and the filenames assosciated with those barcodes

##### PARSE FILES #####

##### CREATE A BARCODE DICT, WITH BARCODE FILES AS TUPLED VALUES, OPEN ALL FILES #####
with open(index_file, 'r') as file:
    header = file.readline().rstrip('\n').split('\t')
    seq_col = header.index('index sequence') # SEQUENCE COL, ONLY COL OF INTEREST

    barcode_bank  = {}
    for line in file:
        parts = line.rstrip('\n').split('\t')
        R1filename = f"{parts[seq_col]}_R1_match.fastq.txt"
        R2filename = f"{parts[seq_col]}_R2_match.fastq.txt"

        r1file = open(R1filename, 'x') # THESE ARE BEING OVERWRITTEN
        r2file = open(R2filename, 'x') # THESE ARE BEING OVERWRITTEN
        barcode_bank[parts[seq_col]] = (r1file, r2file)

### CREATE NON-BARCODE FILES ###
unknown_R1_fastq = open("unknown_R1_fastq.txt", "x")# OPEN OUR UNKNOWN FILES
unknown_R2_fastq = open("unknown_R2_fastq.txt", "x")
hopping_R1_fastq = open("hopping_R1_fastq.txt", "x") # OPEN OUR INDEX HOPPING FILES
hopping_R2_fastq = open("hopping_R2_fastq.txt", "x")

with (
    open(r1, "rt") as file1,
    open(r2, "rt") as file2,
    open(r3, "rt") as file3,
    open(r4, "rt") as file4
):
    file1.readlines()
    file2.readlines()
    file3.readlines()
    file4.readlines()

    n = 0
    bioread_lines = list(range(n, n+4))

    Read1_Bio = [line.strip() for idx, line in enumerate(file1) if idx in bioread_lines]
    Read2_Bio = [line.strip() for idx, line in enumerate(file4) if idx in bioread_lines]
    Read1_Barcode = [line.strip() for idx, line in enumerate(file2) if idx == (n+1)][0] # Second Line
    Read2_Barcode = reverse_comp([line.strip() for idx, line in enumerate(file3) if idx == (n+1)][0]) # Second Line, NOT REVERSE COMP YET
    quality_score = max(mean_qual_score(Read1_Barcode), mean_qual_score(Read2_Barcode))
    # ???? barcode_key = barcode_bank value at key Read1_Barcode ???? 


    Read1_Bio[0] = Read1_Bio[0] + (f" {Read1_Barcode}-{Read2_Barcode}") # ADD Barcodes to the end of the headers
    Read2_Bio[0] = Read2_Bio[0] + (f" {Read1_Barcode}-{Read2_Barcode}")
    # CREATE unknown_R2_fastq.txt
    # CREATE unknown_R1_fastq.txt

    if Read1_Barcode or Read2_Barcode not in barcode_bank or quality_score < quality_score_cutoff:
        unknown_counter += 1
        for line in Read1_Bio:
            unknown_R1_fastq.write(line)
        for line in Read2_Bio:
            unknown_R2_fastq.write(line)
    
    elif Read1_Barcode == Read2_Barcode:
        barcode_bank[Read1_Barcode][0].write(line for lines in (Read1_Bio))
        barcode_bank[Read1_Barcode][1]

 """


        Else, If the R2 barcode matches the R3 barcodes: we know that the reads match valid barcodes:

            For this barcode key, we will add one to the match_dict dictionary score. 

            If <barcode_key>_R1_fastq.txt and <barcode_key>_R2_fastq.txt exist: (each barcodes R1 and R2 fastq files)
                Write Read1_Bio record into the fastq file <barcode_key>_R1_fastq.txt
                Write Read2_Bio record into the fastq file <barcode_key>_R2_fastq.txt
            else:
                Write Read1_Bio record into a new fastq file titled <barcode_key>_R1_fastq.txt
                Write Read2_Bio record into a new fastq file titled <barcode_key>_R2_fastq.txt


        Else: the R2 and R3 Barcodes must not match, both barcodes are in barcode_bank, and quality_score >= quality_score_cutoff, we know that index hopping occured:

            We will add one to the hopping_counter.
            For our Header_Addition, we will add one to the hopping_dict dictionary score or initialize it to 1 if it does not currently exist.

            If hopping_R1_fastq.txt and hopping_R2_fastq.txt exist:
                Write Read1_Bio record into the fastq file hopping_R1_fastq.txt
                Write Read2_Bio record into the fastq file hopping_R2_fastq.txt
            else:
                Write Read1_Bio record into a new fastq file titled hopping_R1_fastq.txt
                Write Read2_Bio record into a new fastq file titled hopping_R2_fastq.txt

##### Output #####

output (match_dict_counter, unknown_counter, hopping_counter, hopping_dict) to the user 
"""