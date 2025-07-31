#!/usr/bin/env python

##### IMPORTS #####
import numpy as np 
import argparse

def reverse_comp(sequence: str) -> str:
    '''Takes a sequence (string) with or without newline characters and returns the reverse compliment of the string without newline characters in uppercase'''
    sequence = sequence.,upper().strip()
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

#Input: IIIIIIIIII
#Expected output: 40

# ADD ARGPARSE AFTER TESTING
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
barcode_bank = "barcode bank file" #Create a barcode_bank with barcode records (known indices) into a dictionary for faster lookup, with key sequence (GTAGCGTA) and value code (B1). 
r1 = 
r2 = 
r3 = 
r4 = 

##### SET GLOBAL VARS #####
match_dict = {} #Create a match_dict dictionary with the barcode or indices as keys, and match numbers initilizaed to 1 upon thier addition.
hopping_counter = 0 #Create a hopping_counter initialized to 0.
hopping_dict = {} #Create a hopping_dict dicitonary with nothing in it. The keys will be our Header_additions and the values will be counters initilaized to 1 upon their addition to the list. 
unkown_counter = 0 #Create a unkown_counter initialized to 0. 

##### PARSE FILES #####
"""
Open all 4 files, R1, R2, R3, R4 at the same time. We will iteratively take 1 record (4 lines) from R1 and R4, while collecting the respective indices ONLY from R2 and R3 (line 2) 

    ### PER RECORD VALUES ###
    Read1_Bio = R1 -> Bio Read 1
    Read1_Barcode = R2 with stripped newline characters -> Index Read 1
    Read2_Bio = R4 -> Bio Read 2
    Read2_Barcode = reverse_comp R3 with stripped newline characters -> Index Read 2 (we must reverse compliment the R3 barcode... as this read pairs with the reverse compliment).
    Header_Addition = "<Read1_Barcode>-<Read2_Barcode>"
    quality_score = minimum of (mean_quality_score(Read1_Barcode), mean_quality_score(Read2_Barcode))
    barcode_key = barcode_bank value at key Read1_Barcode

    Append Header_Addition to the 1st line of Read1_Bio and Read2_Bio. 

        If 1 or both of the barcodes are not in barcode_bank or quality_score < quality_score_cutoff: something went wrong and the barcode/index is unknown, or quality score it too poor to be confident
    
            We will add one to the unknown_counter.

            If unknown_R1_fastq.txt anf unknown_R2_fastq.txt exist:
                Write Read1_Bio record into the fastq file unknown_R1_fastq.txt
                Write Read2_Bio record into the fastq file unknown_R2_fastq.txt
            else:
                Write our Read1_Bio record into a new fastq file titled unknown_R1_fastq.txt
                Write our Read2_Bio record into a new fastq file titled unknown_R2_fastq.txt


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