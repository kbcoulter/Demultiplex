# Demultiplexing

## Overview

Demultiplexing is an essential step for downstream NGS analyses. It separates reads based on their index sequences, ensuring accurate sample identification. This repository provides a **basic demultiplexing workflow** for paired-end sequencing with dual-indexed libraries and Phred+33 encoded FASTQ files. Slurm scripts (`.sh`) are included for running the demultiplexing pipeline on HPC systems. Please reference these scripts to execute the workflow efficiently.

 Given:
- Four input FASTQ files
    - 2 files with biological reads (R1, R2)
    - 2 files with index reads (I1, I2)
- An index file containing known index sequences

Reads are demultiplexed by index pair, producing:

- One R1 FASTQ file and one R2 FASTQ file per matching index pair
- Two FASTQ files for index-hopping (mismatched index pairs)
- Two FASTQ files for unknown or low-quality indexes
    - Includes indexes with Ns or those failing the quality score cutoff
    - Includes indexes not matching the 24 known sequences

For all categories, the index pair sequence is appended to the header of both reads in the FASTQ files.

Example:
```
@SEQ_ID AAAAAAAA-CCCCCCCC
```

This indicates index1 = AAAAAAAA and index2 = CCCCCCCC. If either index is unknown or invalid, the read pair is placed in the unknown category.

Additionally, a Demultiplex_Report.txt is generated (Example: demultiplex_samples/_Demultiplex_Report.txt), summarizing:

- Number of read pairs with properly matched indexes (per index pair)
- Number of read pairs with index-hopping
- Number of index-hopping occurances by index combination
- Number of read pairs with unknown index(es)

## Assumptions
- Paired-end reads
- Phred+33 encoding (+64 Support will be added a later date)
- Tab-separated index file (format shown below)
- Input files are gzipped (Support for uncompressed files will be added at a later date)
- Dual matched indexes (Support for other indexes will be added at a later date)

## Repository Structure

### Non-Code References
```
Pseudocode_Outline.txt
Demultiplex_Notebook.txt
```

### EDA
Scripts for generating per-base quality score distributions for read1, read2, index1, and index2.
Quality scores are averaged at each position across all reads, producing a per-nucleotide mean distribution.
```
eda/qual_by_nuc_slurm.sh
eda/qual_by_nuc.py
```

### Demultiplex
Scripts for demultiplexing FASTQ files and reporting index-hopping.
```
demultiplex_samples/demultiplex.sh
demultiplex_samples/demultiplex.py
```

### Test Demulitplexing Functions
Note: This script only tests functions. While a test input and expected output are provided, these are not used within this script. These will be added to the testing script at a later date. 
```
demultiplex_samples/test_demultiplex_functions.py
```

### My Data
Index File
Submitted 24 indexed libraries for sequencing. Index file location:
```
/projects/bgmp/shared/2017_sequencing/indexes.txt
```

Example:
```
$ head -n4 /projects/bgmp/shared/2017_sequencing/indexes.txt
sample    group    treatment    index    index_sequence
1         2A       control      B1       GTAGCGTA
2         2B       control      A5       CGATCGAT
3         2B       control      C1       GATCAAGG
```

FASTQ Files
```
/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz
/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz
/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz
```
