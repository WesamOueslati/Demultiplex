## File lengths
``` bash
DATA=/projects/bgmp/shared/2017_sequencing/

R1=${DATA}1294_S1_L008_R1_001.fastq.gz
R2=${DATA}1294_S1_L008_R2_001.fastq.gz
R3=${DATA}1294_S1_L008_R3_001.fastq.gz
R4=${DATA}1294_S1_L008_R4_001.fastq.gz

zcat $R# | wc -l
```
- R1: 1452986940 (forward reads)
- R4: 1452986940 (reverse reads)

- R2: 1452986940 (Index 1)
- R3: 1452986940 (Index 2)


## Define the problem
We have 4 files: two containing biological sequences which are the actual reads and two others containing the indexes attached to each of the sequences. One of the reads is the forward read and the other is the reverse read.

So if we have two strands:
NATGTAC    CAACGNCTA.........NNNCATCAG    NATGTAC  
ATACATN    GTTGCNGAT.........NNNGTAGTC    ATACATN

NGCCATG    GCCAGGCTA.........TGGAAAATA    NGCCATG  
ACGGTAN    CGGTCCGAT.........ACCTTTTAT    ACGGTAN

NGGGGGG    TTTTNTTTT.........GGGGNGGGG    NGGGGGG  
CCCCCCN    AAAANAAAA.........CCCCNCCCC    CCCCCCN

Then
    R1 would contain : CAACGNCTA, GCCAGGCTA, TTTTNTTTT
    R2 would contain : NATGTAC, NGCCATG, NGGGGGG
    R3 would contain : NTACATA, NATGGCA, NCCCCCC
    R4 would contain : CTGATGNNN, TATTTTCCA, CCCCNCCCC

We want to separate records into files based on the indexes in order to identify which reads belong to which samples. For each record in the index files, we check if the index is present in our list of indexes. If it is, we check if the index in the first file is the reverse complement of the index in the second file. If that's the case, we add it to the file containing records with that index. Otherwise, if it is present in the list of indexes but is not a reverse complement, then we add both the forward and reverse sequences to the hopped files. If the indexes are not present in the list of sample indexes, then we place the reads in the unknown files.

We expect 2*n + 4 files as the output where n is the number of samples in the multiplexed sequencing run.

Testing:
    Test files should include:
        - one low quality index
        - one mismatching index pair
        - one matching index pair
        - one entry with indexes that are not one of the indexes in the provided file
        - 

For verification:
    We can add the number of lines in the output files for each of the forwards and reverse reads and compare them to the number of lines in the original file.


# Pseudo Code
r1, r2, r3, r4

count_match = 0
count_hop = 0
counts = {index1_index2: 0}
count_unknown = 0

for line in r1, r2, r3, r4:
    every four lines:
        store r1_rec, r2_rec, r3_rec, r4_rec
        store index1 from r2_rec and index2 from r3_rec
        if index1 not in barcodes file or index2 not inbarcodes file;
            add r1_rec and r4_rec to unknown
            count_unknown += 1
            skip to next iteration
        
        counts[ index1_index2 ] +=1
        if index1 is reverse complement of index 2 and is_high_qual(index1) and is_high_qual(index2);
            add r1_rec and r4_rec to index1_index2.fq
            count_match += 1
        else
            add r1_rec and r4_rec to hopped.fq
            count_hop += 1

            
            
            


    record = {header: [3 lines]}

    if idx2 not in indexes or idx2 is low_quality
        write header + index1-index2 to files_r1[unknown]
        write record[header] to files_r1[unknown]
        continue
    if idx1 != rev_comp(idx2)
        write record to files_r1[hopped]
    else:
        write record to files_r1[index]


# Function declarations:
```
def reverse_complement(seq: str) -> str:
    ''' Takes in a DNA or RNA sequence and returns the reverse complement of that sequence.
    '''
    return reverse_comp
Input: ATCG
Output: CGAT

Input: NAAGTC
Output: GACTTN    


def is_low_qual(qual_seq: str, thresh: int) -> bool:
    '''returns true if the average quality score for seq is lower than the threshold, otherwise return false.'''
    return flag
Input: IIIII, 30
Output: False

Input: AAAAA, 30
Output: True
```