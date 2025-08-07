#!/usr/bin/env python

##### IMPORTS #####
import numpy as np 
import argparse
import gzip
import os
import sys
sys.path.append(os.path.abspath("..")) # ONLY TO GET BIOINFO AND MAINTAIN VERSION W/ BIOINFO USED IN ASSIGN. 1ST
import bioinfo

##### FUNCTIONS #####
def get_args():
     parser = argparse.ArgumentParser(description="A script to average the quality scores at each position for all reads and generate a per nucleotide mean distribution of quality scores for read1, read2, index1, and index2.")
     parser.add_argument("-r1", "--read1", help="An input FASTQ Read 1 file", required=True, type=str)
     parser.add_argument("-r2", "--read2", help="An input FASTQ Read 2 file", required=True, type=str)
     parser.add_argument("-r3", "--read3", help="An input FASTQ Read 3 file", required=True, type=str)
     parser.add_argument("-r4", "--read4", help="An input FASTQ Read 4 file", required=True, type=str)
     parser.add_argument("-i", "--index", help="An index tsv containing all indexes(barcodes) under the column 'index sequence'", required=True, type=str)
     parser.add_argument("-o", "--outname", help="A name for the output dir, defaults to 'Demultiplexed'", required=False, type=str, default = "Demultiplexed")
     parser.add_argument("-q", "--qualcutoff", help="An optional quality score cutoff, defaults to 5", required=False, type=int, default = 5)
     ### ADD LATER ### parer.add_argument("-z", "--unzipped", help=" Are FASTQs unzipped (not gzipped)? Defaults to zipped FASTQ (False), otherwise, set to True.", required=False, type= bool, default = True)
     return parser.parse_args()

def reverse_comp(sequence: str) -> str | None:
    '''Takes a sequence (string) with or without newline characters and returns the reverse compliment of the string without newline characters in uppercase'''
    sequence = sequence.upper().strip()
    map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'} # MAP Ns to maintain them to break.
    rev_comp = ("".join([map.get(entry, entry)for entry in sequence]))[::-1] # BRACKET REVERSES, OTHER COMPS
    return rev_comp

def mean_qual_score(sequence: str) -> float:
    '''Takes a quality score sequence (string), uses convert_phred function to convert each phred+33 encoding into a integer quality score, then averages the quality score and returns the mean values. '''
    holder = []
    for char in sequence:
        holder.append(bioinfo.convert_phred(char))
    return (np.mean(holder))

def collect_record(file, indexread, headeradd="") -> str | None | list:
    '''Takes a FASTQ file and returns the entire record if indexread is False, or the index line only if indexread is True. Strips newline chars'''
    if indexread:
        line_junk = file.readline()
        if not line_junk:
            return (None)
        line2 = file.readline()
        line_junk = file.readline()
        line_junk = file.readline() # Junk line 3 and 4
        return (line2.strip())

    else:
        record = []
        for line_counter in range(4):
            line = file.readline()
            if not line:
                return (None)
            record.append(line.strip())
        record[0] += (f" {headeradd}")
        return record
    
def main():
    args = get_args()
    outname = args.outname
    index_file = args.index
    quality_score_cutoff = args.qualcutoff
    r1 = args.read1
    r2 = args.read2
    r3 = args.read3
    r4 = args.read4
    
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
        seq_col = header.index('index sequence')

        barcode_bank  = {}
        for line in file:
            parts = line.rstrip('\n').split('\t')
            R1filename = f"{outname}/{parts[seq_col]}_R1.fastq.txt"
            R2filename = f"{outname}/{parts[seq_col]}_R2.fastq.txt"
            r1file = open(R1filename, 'x') # THESE ARE BEING OVERWRITTEN
            r2file = open(R2filename, 'x') 
            barcode_bank[parts[seq_col]] = (r1file, r2file)

    ### CREATE NON-BARCODE FILES ###
    unknown_R1_fastq = open(f"{outname}/unknown_R1_fastq.txt", "x")# OPEN OUR UNKNOWN FILES
    unknown_R2_fastq = open(f"{outname}/unknown_R2_fastq.txt", "x")
    hopping_R1_fastq = open(f"{outname}/hopping_R1_fastq.txt", "x") # OPEN OUR INDEX HOPPING FILES
    hopping_R2_fastq = open(f"{outname}/hopping_R2_fastq.txt", "x")

    with gzip.open(r1, 'rt') as file1, gzip.open(r2, 'rt') as file2, gzip.open(r3, 'rt') as file3, gzip.open(r4, 'rt') as file4: #Zipped
    ### ADD BACK WIT UNZIPPED SUPPORT ### with open(r1, 'rt') as file1, open(r2, 'rt') as file2, open(r3, 'rt') as file3, open(r4, 'rt') as file4: #Unzipped
        # print("\nDemultiplexing...\n") # Commented out temporarily to avoid writing over output
        while True:
            Read1_Barcode = (collect_record(file2, True))
            if not Read1_Barcode:
                break # We can break when we don't have barcode
            Read2_Barcode = reverse_comp(collect_record(file3, True)) # TAKE THE REVERSE COMP HERE
            indextoindexkey = (f"{Read1_Barcode}-{Read2_Barcode}")
            Read1_Bio = collect_record(file1, False, indextoindexkey)
            Read2_Bio = collect_record(file4, False, indextoindexkey)
            quality_score = min(mean_qual_score((Read1_Barcode)), mean_qual_score(Read2_Barcode)) #Ignore Warn.

            if (Read1_Barcode not in barcode_bank or Read2_Barcode not in barcode_bank or quality_score < quality_score_cutoff):
                unknown_counter += 1
                for line in Read1_Bio: #Ignore Warn.
                    unknown_R1_fastq.write(line + '\n')
                for line in Read2_Bio:
                    unknown_R2_fastq.write(line + '\n')

            elif Read1_Barcode == Read2_Barcode:
                match_counter += 1
                match_dict[indextoindexkey] = match_dict.get(indextoindexkey, 0) + 1
                for line in Read1_Bio: #Ignore Warn.
                    barcode_bank[Read1_Barcode][0].write(line + '\n')
                for line in Read2_Bio:
                    barcode_bank[Read1_Barcode][1].write(line + '\n')

            else:
                hopping_counter += 1
                hopping_dict[indextoindexkey] = hopping_dict.get(indextoindexkey, 0) + 1
                for line in Read1_Bio: #Ignore Warn.
                    hopping_R1_fastq.write(line + '\n')
                for line in Read2_Bio:
                    hopping_R2_fastq.write(line + '\n')

    ##### OUTPUT REPORT #####
    total = match_counter + hopping_counter + unknown_counter
    report = open(f"{outname}/_Demultiplex_Report.txt", "w")
    report.write((46 + len(outname)) * "#" + "\n" + f"### DEMULTIPLEX REPORT SUMMARY (FASTQ In {outname}) ###\n" + (46 + len(outname)) * "#" + "\n\n")
    report.write("########## GENERAL SUMMARY ##########\n")
    report.write(f"Matching Indexes:\t{match_counter} Matching Records:\t{round((match_counter * 100 / total),3)}% of Total\n")
    report.write(f"Hopping Indexes:\t{hopping_counter} Hopping Records:\t{round((hopping_counter * 100 / total),3)}% of Total\n")
    report.write(f"Unknown Indexes:\t{unknown_counter} Unknown Records:\t{round((unknown_counter * 100 / total),3)}% of Total\n")

    report.write(f"Total Reads:\t{total}\nSelected Quality Score:\t{quality_score_cutoff}\nIndexes Colleted From {index_file}:\n")
    for entry in barcode_bank:
        report.write(entry + "\n")

    report.write("\n########## READ PAIR RESULTS ##########\n##### MATCH RESULTS #####\nIndex_Pair\tNum_Pairs\tPerc_of_Matched\tPerc_of_Total\n")
    match_dict_sorted = sorted(match_dict.items(), key=lambda item: item[1])
    for pair,count in match_dict_sorted:
        report.write(f"{pair}:\t{count}\t{round((count * 100 / match_counter),3)}\t{round((count * 100 / total),3)}\n")

    report.write("\n##### HOPPING RESULTS #####\nIndex_Pair\tNum_Pairs\tPerc_of_hoppinged\tPerc_of_Total\n")
    hopping_dict_sorted = sorted(hopping_dict.items(), key=lambda item: item[1])
    for pair,count in hopping_dict_sorted:
        report.write(f"{pair}:\t{count}\t{round((count * 100 / hopping_counter),3)}\t{round((count * 100 / total),3)}\n")

    ##### CLOSE ALL FILES #####
    unknown_R1_fastq.close()
    unknown_R2_fastq.close()
    hopping_R1_fastq.close()
    hopping_R2_fastq.close()

    for r1file, r2file in barcode_bank.values():
        r1file.close()
        r2file.close()

    # print(f"Demultiplexing Complete. Output contained in {outname}") # Commented out temporarily to avoid writing over output

if __name__ == "__main__":
    main()