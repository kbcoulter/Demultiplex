#!/usr/bin/env python

##### IMPORTS #####
import numpy as np 
import argparse
import gzip
import os

def reverse_comp(sequence: str) -> str:
    '''Takes a sequence (string) with or without newline characters and returns the reverse compliment of the string without newline characters in uppercase'''
    sequence = sequence.upper().strip()
    map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',}
    rev_comp = ("".join([map.get(entry, entry)for entry in sequence]))[::-1] # BRACKET REVERSES, OTHER COMPS
    return rev_comp
# Input: CATG
# Expected Output: CATG

def convert_phred(letter: str) -> int:
    '''Takes a single ASCII character (string) encoded in Phred+33 and returns the quality score value as an integer.'''
    return (ord(letter) - 33)
#Input: I
#Expected output: 40

def mean_qual_score(sequence: str) -> float:
    '''Takes a quality score sequence (string), uses convert_phred function to convert each phred+33 encoding into a integer quality score, then averages the quality score and returns the mean values. '''
    holder = []
    for char in sequence:
        holder.append(convert_phred(char))
    return (np.mean(holder))

def collect_record(file, indexread) -> list:
    '''Takes a FASTQ file and returns the entire record if indexread is False, or the index line only if indexread is True. Strips newline chars'''
    if indexread:
        line_junk = file.readline()
        if not line_junk:
            return ([])
        line2 = file.readline()
        line_junk = file.readline()
        line_junk = file.readline() # Junk line 3 and 4
        return (line2.strip())

    else:
        record = []
        for line_counter in range(4):
            line = file.readline()
            if not line:
                return ([])
            record.append(line.strip())
        return record


# !!! ADD ARGPARSE AFTER TESTING !!!
"""
def get_args():
     parser = argparse.ArgumentParser(description="A script to average the quality scores at each position for all reads and generate a per nucleotide mean distribution of quality scores for read1, read2, index1, and index2.")
     parser.add_argument("-r1", "--read1", help="An input FASTQ Read 1 file", required=True, type=str)
     parser.add_argument("-r2", "--read2", help="An input FASTQ Read 2 file", required=True, type=str)
     parser.add_argument("-r3", "--read3", help="An input FASTQ Read 3 file", required=True, type=str)
     parser.add_argument("-r4", "--read4", help="An input FASTQ Read 4 file", required=True, type=str)
     parser.add_argument("-i", "--index", help="An index tsv containing all indexes(barcodes) under the column 'index sequence'", required=True, type=str)
     #parser.add_argument("-k", "--readsize", help="Size (length) of the reads", required=True, type=int)
     parser.add_argument("-o", "--outname", help="A name for the output dir, defaults to 'demultiplexed'", required=False, type=str)
     parser.add_argument("-q", "--qualcutoff", help="An optional quality score cutoff, defaults to 5", required=False, type= int)
     #parer.add_argument("-z", "--unzipped", help=" Are FASTQs unzipped (not gzipped)? Defaults to zipped FASTQ (False), otherwise, set to True.", required=False, type= bool, default = True)
     return parser.parse_args()

args = get_args()
"""
### REAL DATA ###
#r1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
#r2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
#r3 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
#r4 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"

### TESTING DATA ### 
quality_score_cutoff = 5 #Create a quality_score_cutoff integer where every index/barcode read with quality score below goes into unknown. 
index_file = "/projects/bgmp/shared/2017_sequencing/indexes.txt"
r1 = "/TEST-input_FASTQ/S1_TEST_R1_001.fastq"
r2 = "/TEST-input_FASTQ/S1_TEST_R1_002.fastq"
r3 = "/TEST-input_FASTQ/S1_TEST_R1_003.fastq"
r4 = "/TEST-input_FASTQ/S1_TEST_R1_004.fastq"
outname = "testing_output"

##### SET GLOBAL VARS #####
match_dict = {} # Dictionary with the barcode or indices as keys, and match numbers initilizaed to 1 upon thier addition.
match_counter = 0 # Create a _counter initialized to 0.
hopping_counter = 0
hopping_dict = {}
unknown_counter = 0

os.makedirs(outname, exist_ok=True)

##### CREATE A BARCODE DICT, WITH BARCODE FILES AS TUPLED VALUES, OPEN ALL FILES #####
with open(index_file, 'r') as file:
    header = file.readline().rstrip('\n').split('\t')
    seq_col = header.index('index sequence') # SEQUENCE COL, ONLY COL OF INTEREST

    barcode_bank  = {}
    for line in file:
        parts = line.rstrip('\n').split('\t')
        R1filename = f"{outname}/{parts[seq_col]}_R1_match.fastq.txt"
        R2filename = f"{outname}/{parts[seq_col]}_R2_match.fastq.txt"

        r1file = open(R1filename, 'x') # THESE ARE BEING OVERWRITTEN
        r2file = open(R2filename, 'x') 
        barcode_bank[parts[seq_col]] = (r1file, r2file)

### CREATE NON-BARCODE FILES ###
unknown_R1_fastq = open(f"{outname}/unknown_R1_fastq.txt", "x")# OPEN OUR UNKNOWN FILES
unknown_R2_fastq = open(f"{outname}/unknown_R2_fastq.txt", "x")
hopping_R1_fastq = open(f"{outname}/hopping_R1_fastq.txt", "x") # OPEN OUR INDEX HOPPING FILES
hopping_R2_fastq = open(f"{outname}/hopping_R2_fastq.txt", "x")

#file1, = open(r1, "rt") #file2 = open(r2, "rt") #file3 = open(r3, "rt") #file4 = open(r4, "rt")

#with gzip.open(r1, 'rt') as file1, gzip.open(r2, 'rt') as file2, gzip.open(r3, 'rt') as file3, gzip.open(r4, 'rt') as file4: #Zipped
with open(r1, 'rt') as file1, open(r2, 'rt') as file2, open(r3, 'rt') as file3, open(r4, 'rt') as file4: #Unzipped
    while True:
        Read1_Bio = collect_record(file1, False)
        Read2_Bio = collect_record(file4, False)
        Read1_Barcode = (collect_record(file2, True)[0])
        Read2_Barcode = reverse_comp(collect_record(file3, True)[0]) # TAKE THE REVERSE COMP HERE
        indextoindexkey = (f" {Read1_Barcode}-{Read2_Barcode}")

        if not (Read1_Barcode or Read2_Barcode or Read1_Bio or Read2_Bio): # Break as soon as we get nothing returned -> Read all records
            break
        
        quality_score = min(mean_qual_score(Read1_Barcode), mean_qual_score(Read2_Barcode))

        Read1_Bio[0] = Read1_Bio[0] + (indextoindexkey) # ADD Barcodes to the end of the headers
        Read2_Bio[0] = Read2_Bio[0] + (indextoindexkey)

        if Read1_Barcode or Read2_Barcode not in barcode_bank or quality_score < quality_score_cutoff:
            unknown_counter += 1
            for line in Read1_Bio:
                unknown_R1_fastq.write(line)
            for line in Read2_Bio:
                unknown_R2_fastq.write(line)
        
        elif Read1_Barcode == Read2_Barcode:
            match_counter += 1
            match_dict[indextoindexkey] = match_dict.get(indextoindexkey, 0) + 1

            barcode_bank[Read1_Barcode][0].write(line for lines in (Read1_Bio)) # HMM
            barcode_bank[Read1_Barcode][1].write(line for lines in (Read2_Bio))

        else:
            hopping_counter += 1
            hopping_dict[indextoindexkey] = hopping_dict.get(indextoindexkey, 0) + 1 # incr hopping_dict counter

            for line in Read1_Bio:
                hopping_R1_fastq.write(line)
            for line in Read2_Bio:
                hopping_R2_fastq.write(line)

        
##### Output #####
total = match_counter + hopping_counter + unknown_counter
report = open(f"{outname}/Demultiplex_Report.txt", "w")
report.write("#######################################################################")
report.write(f"### DEMULTIPLEX REPORT SUMMARY (Output Contained in: {outname} ###")
report.write("#######################################################################\n")

report.write("########## GENERAL SUMMARY ##########")
report.write(f"Matching Indexes:\t{match_counter} Matching Records:\t{round((match_counter * 100 / total),3)}% of Total")
report.write(f"Hopping Indexes:\t{hopping_counter} Hopping Records:\t{round((hopping_counter * 100 / total),3)}% of Total")
report.write(f"Unknown Indexes:\t{unknown_counter} Unknown Records:\t{round((unknown_counter * 100 / total),3)}% of Total")
report.write(f"Total Reads:\t{total}")
report.write(f"Selected Quality Score:\t{quality_score_cutoff}")
report.write(f"Indexes Colleted From:\t{index_file}\n")

report.write("########## READ PAIR RESULTS ##########\n")
report.write("### MATCH RESULTS ###")
report.write("Index_Pair\tNum_Pairs\tPerc_of_Matched\tPerc_of_Total")
for entry in match_dict:
    report.write(f"{entry}:\t{match_dict[entry]}\t{round((match_dict[entry] * 100 / match_counter),3)}\t{round((match_dict[entry] * 100 / total),3)}")

report.write("### HOPPING RESULTS ###")
report.write("Index_Pair\tNum_Pairs\tPerc_of_hoppinged\tPerc_of_Total")
for entry in hopping_dict:
    report.write(f"{entry}:\t{hopping_dict[entry]}\t{round((hopping_dict[entry] * 100 / hopping_counter),3)}\t{round((hopping_dict[entry] * 100 / total),3)}")

##### CLOSE ALL FILES #####
# ITER ALL FILES AND CLOSE (FUNCTION UNNECESSARY)