#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log

#fqR1="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
#zcat "$fqR1" > R1FQ
#/usr/bin/time -v ./qual_by_nuc.py -f R1FQ -k 101 -o R1_FASTQ
#rm R1FQ

#fqR2="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
#zcat "$fqR2" > R2FQ
#/usr/bin/time -v ./qual_by_nuc.py -f R2FQ -k 8 -o R2_FASTQ
#rm R2FQ

#fqR3="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
#zcat "$fqR3" > R3FQ
#/usr/bin/time -v ./qual_by_nuc.py -f R3FQ -k 8 -o R3_FASTQ
#rm R3FQ

fqR4="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"
zcat "$fqR4" > R4FQ
/usr/bin/time -v ./qual_by_nuc.py -f R4FQ -k 101 -o R4_FASTQ
rm R4FQ