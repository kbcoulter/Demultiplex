#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --job-name=demultiplex.sh
#SBATCH --error=demultiplex_error.log
#SBATCH --output=demultiplex_output.log

read1=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz
read2=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz
read3=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
read4=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz
index=/projects/bgmp/shared/2017_sequencing/indexes.txt

mamba activate demultiplex
/usr/bin/time -v ./demultiplex.py -r1 $read1 -r2 $read2 -r3 $read3 -r4 $read4 -i $index -q 5 -o Demultiplexed