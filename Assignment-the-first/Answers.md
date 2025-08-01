# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it [here](./qual_by_nuc.py).

    Also, here is [Slurm Script](./qual_by_nuc_slurm.sh) to Run the python (qual_by_nuc.py). 

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Bio Read 1 | 101 | 33 |
| 1294_S1_L008_R2_001.fastq.gz | Index Read 1 | 8 | 33 |
| 1294_S1_L008_R3_001.fastq.gz | Index Read 2 | 8 | 33 |
| 1294_S1_L008_R4_001.fastq.gz | Bio Read 2 | 101 | 33 |

2. Per-base NT distribution
    1. [Histogram 1](./R1_FASTQ.png)
    2. [Histogram 2](./R2_FASTQ.png)
    3. [Histogram 3](./R3_FASTQ.png)
    4. [Histogram 4](./R4_FASTQ.png)

3. What is a good quality score cutoff for index reads and biological read pairs to utilize for sample identification and downstream analysis, respectively? Justify your answer.

    It is essential that we maintain a sufficient number of data while ensuring the quality of our indexed reads. Initially, it may seem logical to set a relatively high threshold for the mean quality score in our index. However, this is unecessary.

    Index hopping or non-matching barcodes will not appear in the output barcode FASTQ files, so there’s no need to set a quality score threshold to filter them out. Instead, the quality score will be used to filter reads that do match barcodes in our dataset. If a read matches a barcode in our database, it’s reasonable to assume it’s correct, as the likelihood of both barcodes index hopping to both match, and exist in our index bank, is incredibly low. Therefore, the quality score cutoff should only be used to exclude reads with EXTREMELY poor average quality, as the definitive errors are being sorted out already. 

    I propose a quality score cutoff of 5, meaning that if the mean quality score of a index is below 5, the read isn't passed into its respective FASTQ file and goes into an unknown file instead. This threshold should only cut reads of the worst quality. 

4. How many indexes have undetermined (N) base calls? (Utilize your command line tool knowledge. Submit the command(s) you used. CHALLENGE: use a one-line command)

        FASTQ R2: 3976613 via $ zcat 1294_S1_L008_R2_001.fastq.gz | awk '/^@/ {getline; print}' | grep N | wc -l
        FASTQ R3: 3328051 via $ zcat 1294_S1_L008_R3_001.fastq.gz | awk '/^@/ {getline; print}' | grep N | wc -l

## Part 2
1. Define the problem

    We have lane sequencing data from an Illumina sequencer and need to demultiplex the data , meaning we need to separate it into distinct files for our different biological reads. Currently, we have four files: Bioread1, Bioread2, Index1, and Index2. When matching the indexes (or barcodes) to their respective reads, you end up with two records of the same read, one being the reverse complement, with both reads oriented 5’ to 3’. So, we expect a corresponding record between Read1 and Read2. To ensure no index hopping has occurred, we need to match the indexes across these FASTQ files, and then split these large files into smaller FASTQ files for each matched read. That is, barcode A matches would produce two FASTQ files (Read1 and Read2), barcode B matches would produce another two files, etc. To properly separate these reads, we must match the Read1 and Read2 indexes (using the reverse complement of the Read2 index) for each FASTQ record. If index hopping occured, the mean quality score is not sufficient, or if the barcode doesn't match our key, we have additional output files for these records and counters to track the number of times this occurs. 

2. Describe output

    For each matching index pair, there is one Read1 FASTQ file and one Read2 FASTQ file.
    For index-hopping events (non-matching index pairs within our barcode set), there are two additional FASTQ files, one for Read1 (hopping_R1_fastq.txt) and one for Read2 (hopping_R2_fastq.txt).
    Two more FASTQ files are generated when one or both index reads are unknown (contain Ns), fail the quality score threshold, or do not exist in our index set

    This results in the following output:

        <barcode>_R1_fastq.txt for each unique barcode/index
        <barcode>_R2_fastq.txt for each unique barcode/index 
        hopping_R1_fastq.txt
        hopping_R2_fastq.txt
        unknown_R1_fastq.txt
        unknown_R2_fastq.txt
        
    Also, We must output:
    The number of read-pairs with properly matched indexes (per index-pair), the number of read pairs with index-hopping observed, the number of index-hopping occurances by index combination, and the number of read-pairs with unknown index(es).

3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).

4. [Pseudocode](./Pseudocode_Outline.txt)
    ( Click the link ) 
    
5. High level functions. (Also at top of [Pseudocode](./Pseudocode_Outline.txt)) For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functionsxit
    4. Return statement

    ```
    def reverse_comp(sequence: str) -> str:
        '''Takes a sequence (string) with or without newline characters and returns the 
            reverse compliment of the string without newline characters'''
        return rev_cop
    Input: CATG
    Expected Output: GTAC

    def convert_phred(letter: str) -> int:
        '''Takes a single ASCII character (string) encoded in Phred+33 and returns the 
            quality score value as an integer.'''
        return qscore  
    Input: I
    Expected output: 40

    def mean_qual_score(sequence: str) -> float:
        '''Takes a quality score sequence (string), uses convert_phred function to convert each phred+33 
            encoding into a integer quality score, then averages the quality score and returns the mean values. '''
        return mean_score
    Input: IIIIIIIIII
    Expected output: 40
    ```



