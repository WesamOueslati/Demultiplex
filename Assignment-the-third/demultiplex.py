#!/usr/bin/env python

import bioinfo
import gzip
import json

def reverse_complement(seq: str):
    tr_table = str.maketrans('ATCGN', 'TAGCN')
    reverse = seq[::-1]
    reverse_comp = reverse.translate(tr_table)
    return reverse_comp

def is_high_qual(qual_seq: str ,seq , thresh: float):
    return all(([bioinfo.convert_phred(char) > thresh for char in qual_seq])) and ('N' not in seq)

def is_valid(barcode, barcodes):
    return (barcode in barcodes.values())

def open_files(barcode_pairs: list[str]) -> dict:
    handles = {}
    
    for pair in barcode_pairs:
        handles[f"{pair}_r1"] = open(f"output/{pair}_r1.fq" ,'a')
        handles[f"{pair}_r2"] = open(f"output/{pair}_r2.fq" ,'a')
        
    handles[f"unknown_r1"] = open(f"output/unknown_r1.fq", 'a')
    handles[f"unknown_r2"] = open(f"output/unknown_r2.fq", 'a')
    handles[f"hopped_r1"] = open(f"output/hopped_r1.fq", 'a')
    handles[f"hopped_r2"] = open(f"output/hopped_r2.fq", 'a')

    return handles

def close_files(handles) -> None:
    for handle in handles:
        handles[handle].close()
    return

def write_record(handle, record: list[str], idx_str: str) -> None:
    for i, line in enumerate(record): 
        if i%4 == 0:
            handle.write(f"{line} {idx_str}\n")
        else:
            handle.write(f"{line}\n")
    return
    
def barcodes_to_dict(barcodes_file, sep):
    barcodes = { }
    with open(barcodes_file, 'r') as fh:
        for line in fh:
            if line.startswith('sample'):
                continue
            line = line.rstrip('\n')
            ls = line.split(sep)
            label = ls[3]
            barcode = ls[4]
            barcodes[label] = barcode
    return barcodes

def store_barcode_pairs(barcodes, output):
    barcodes_ls = barcodes.values()

    counts = {}
    
    files = []

    prod = itertools.product(barcodes_ls, repeat=2)
    with open(output, 'w') as out:
        for pair in prod:
            if pair[0] == pair[1]:
                files.append(f"{pair[0]}_{pair[1]}" )
            pair_str = f"{pair[0]}_{pair[1]}" 
            out.write(pair_str + "\n")
            counts[pair_str] = 0
    return counts, files
    

def demultiplex(R1: str, R2: str, R3: str, R4: str, barcodes: dict, output) -> dict:
    counts, files = store_barcode_pairs(barcodes, output)
    handles = open_files(files)
    with gzip.open(R1, 'rt') as r1, gzip.open(R2, 'rt') as r2, gzip.open(R3, 'rt') as r3, gzip.open(R4, 'rt') as r4:
        count_matches = 0
        count_unknown = 0
        count_hopped = 0
        r1_rec = []
        r2_rec = []
        r3_rec = []
        r4_rec = []
        
        index_headers = []
        sequence_headers = []
        
        for line in r1:
            
            r1_line = line.strip('\n')
            r2_line = r2.readline().strip('\n')
            r3_line = r3.readline().strip('\n')
            r4_line = r4.readline().strip('\n')
            
            r1_rec.append(r1_line)
            r2_rec.append(r2_line)
            r3_rec.append(r3_line)
            r4_rec.append(r4_line)
            
            

            if len(r1_rec) == 4:
                head_1, head_2 = r1_rec[0], r2_rec[0]
                idx_1, idx_2 = r2_rec[1], r3_rec[1]
                qual_1, qual_2 = r2_rec[3], r3_rec[3] 
                             
                index_headers.append(len(head_1))
                sequence_headers.append(len(head_2))
                                
                idx_str = f"{idx_1}_{reverse_complement(idx_2)}"
                
                rv = (idx_1 == reverse_complement(idx_2))
                # print(f"is reversed?: {rv}")
                low_qual = not (is_high_qual(qual_1, idx_1, 0.0) and is_high_qual(qual_2, idx_2, 0.0))
                # print(f"is low_quality?: {low_qual}")
                valid = (is_valid(idx_1, barcodes) and is_valid(reverse_complement(idx_2), barcodes))
                # print(f"is valid?: {valid}")
                
                if not valid or low_qual:
                    fh_r1 = handles[f"unknown_r1"]
                    fh_r2 = handles[f"unknown_r2"]
                    
                    write_record(fh_r1, r1_rec, idx_str)
                    write_record(fh_r2, r4_rec, idx_str)
                    count_unknown += 1
                else:
                    if rv:
                        fh_r1 = handles[f"{idx_str}_r1"]
                        fh_r2 = handles[f"{idx_str}_r2"]
                        
                        write_record(fh_r1, r1_rec, idx_str)
                        write_record(fh_r2, r4_rec, idx_str)
                        counts[idx_str] += 1
                        count_matches += 1
                    else:
                        fh_r1 = handles["hopped_r1"]
                        fh_r2 = handles["hopped_r2"]
                        
                        write_record(fh_r1, r1_rec, idx_str)
                        write_record(fh_r2, r4_rec, idx_str)
                        counts[idx_str] += 1
                        count_hopped += 1
                
                r1_rec.clear()
                r2_rec.clear()
                r3_rec.clear()
                r4_rec.clear()
    close_files(handles)
   
    return counts, count_matches, count_unknown, count_hopped
if __name__ == "__main__":
    
    import itertools
    import argparse
    
    assert is_high_qual("IIII", "AAAN", 5.0) == False
    assert is_high_qual("IIII", "AAAA", 5.0) == True
    
    BARCODES_PATH="/projects/bgmp/shared/2017_sequencing/indexes.txt"
    BARCODE_PAIRS_PATH="/projects/bgmp/oueslati/bioinfo/Bi622/AS1/Demultiplex/barcode_pairs.txt"
    
    # R1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
    # R2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
    # R3 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
    # R4 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--forward_read', help="Forward read")
    parser.add_argument('-1', '--index_1', help="index 1")
    parser.add_argument('-2', '--index_2', help="index 2")
    parser.add_argument('-r', '--reverse_read', help="Forward read")
    # parser.add_argument('-b', '--barcodes', help="Barcodes file")
    
    args = parser.parse_args()
    
    R1 = args.forward_read
    R2 = args.index_2
    R3 = args.index_1
    R4 = args.reverse_read

    # BARCODES_PATH = "/projects/bgmp/oueslati/bioinfo/Bi622/AS1/Demultiplex/TEST-input_FASTQ/barcodes.txt"
    # BARCODE_PAIRS_PATH="/projects/bgmp/oueslati/bioinfo/Bi622/AS1/Demultiplex/barcode_pairs.txt"
    
    # r1 = "/projects/bgmp/oueslati/bioinfo/Bi622/AS1/Demultiplex/TEST-input_FASTQ/r1.fq"
    # r2 = "/projects/bgmp/oueslati/bioinfo/Bi622/AS1/Demultiplex/TEST-input_FASTQ/r2.fq"
    # r3 = "/projects/bgmp/oueslati/bioinfo/Bi622/AS1/Demultiplex/TEST-input_FASTQ/r3.fq"
    # r4 = "/projects/bgmp/oueslati/bioinfo/Bi622/AS1/Demultiplex/TEST-input_FASTQ/r4.fq"
    
    # barcodes = barcodes_to_dict(BARCODES_PATH, sep=' ')
    
    barcodes = barcodes_to_dict(BARCODES_PATH, sep='\t')
    
    # demultiplex(r1, r2, r3, r4, barcodes, counts)
    counts, count_matches, count_unknown, count_hopped = demultiplex(R1, R2, R3, R4, barcodes, BARCODE_PAIRS_PATH)
    
    # Saving output to a file so i can use to make heat map later
    with open("counts.json", "w") as file:
        json.dump(counts, file, indent=4)
        
    with open("counts.txt", 'w') as counts_fh:
        counts_fh.write(f"Number of barcodes that matched:\t{count_matches}\n")
        counts_fh.write(f"Number of barcodes that were different:\t{count_hopped}\n")
        counts_fh.write(f"Number of barcodes that were not valid barcodes or were low quality:\t{count_unknown}/n")
        
    