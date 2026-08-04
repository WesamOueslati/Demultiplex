#!/usr/bin/env python

import bioinfo
import gzip

def reverse_complement(seq: str):
    tr_table = {
        'T': 'A',
        'C': 'G',
        'A': 'T',
        'G': 'C',
        'N': 'N'
    } 
    reverse = seq[::-1]
    reverse_comp = reverse.translate(tr_table)
    return reverse_comp

def is_high_qual(qual_seq: str ,seq , thresh: float):
    return all(([bioinfo.convert_phred(char) > thresh for char in qual_seq])) and ('N' not in seq)

def is_valid(barcode, barcodes):
    return (barcode in barcodes)

def write_record(file_name: str, record, idx_1, idx_2):
    fh = open(file_name, 'a')
    for i, line in enumerate(record): 
        if i%4 == 0:
            fh.write(f"{line} {idx_1}_{idx_2}\n")
        else:
            fh.write(f"{line}\n")
    fh.close()
    
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
    
    prod = itertools.product(barcodes_ls,barcodes_ls)
    with open(output, 'w') as out:
        for pair in prod:
            pair_str = f"{pair[0]}_{pair[1]}" 
            out.write(pair_str + "\n")
            counts[pair_str] = 0
    return counts
    

def demultiplex(R1: str, R2: str, R3: str, R4: str, barcodes: dict, counts) -> None:
    
    with open(R1, 'r') as r1, open(R2, 'r') as r2, open(R3, 'r') as r3, open(R4, 'r') as r4:
        count_matches = 0
        count_unknown = 0
        count_hopped = 0
        r1_rec = []
        r2_rec = []
        r3_rec = []
        r4_rec = []
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
                seq_1, seq_2 = r1_rec[0], r4_rec[0]
                idx_1, idx_2 = r2_rec[1], r3_rec[1]
                qual_1, qual_2 = r2_rec[3], r3_rec[3] 
                
                idx_str = f"{idx_1}_{idx_2}"
                
                rv = (idx_1 == reverse_complement(idx_2))
                low_qual = not (is_high_qual(qual_1, seq_1, 5.0) and is_high_qual(qual_2, seq_2, 5.0))
                valid = (is_valid(idx_1, barcodes) and is_valid(idx_2, barcodes))
                
                if not valid or low_qual:
                    write_record("unknown_r1", r1_rec, idx_1, idx_2)
                    write_record("unknown_r2", r4_rec, idx_1, idx_2)
                    count_unknown += 1
                else:
                    if rv:
                        write_record(f"{idx_str}_r1.fq", r1_rec, idx_1, idx_2)
                        write_record(f"{idx_str}_r2.fq", r4_rec, idx_1, idx_2)
                        counts[idx_str] += 1
                        count_matches += 1
                    else:
                        write_record("hopped_r1", r1_rec, idx_1, idx_2)
                        write_record("hopped_r2", r4_rec, idx_1, idx_2)
                        counts[idx_str] += 1
                        count_hopped += 1
                
                r1_rec.clear()
                r2_rec.clear()
                r3_rec.clear()
                r4_rec.clear()

if __name__ == "__main__":
    
    import itertools
    
    assert is_high_qual("IIII", "AAAN", 5.0) == False
    assert is_high_qual("IIII", "AAAA", 5.0) == True
    
    # BARCODES_PATH="/projects/bgmp/shared/2017_sequencing/indexes.txt"
    # BARCODE_PAIRS_PATH="/projects/bgmp/oueslati/bioinfo/Bi622/AS1/Demultiplex/barcode_pairs.txt"
    
    # R1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
    # R2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
    # R3 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
    # R4 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"
    
    
    BARCODES_PATH = "/projects/bgmp/oueslati/bioinfo/Bi622/AS1/Demultiplex/TEST-input_FASTQ/barcodes.txt"
    BARCODE_PAIRS_PATH="/projects/bgmp/oueslati/bioinfo/Bi622/AS1/Demultiplex/barcode_pairs.txt"
    
    r1 = "/projects/bgmp/oueslati/bioinfo/Bi622/AS1/Demultiplex/TEST-input_FASTQ/r1.fq"
    r2 = "/projects/bgmp/oueslati/bioinfo/Bi622/AS1/Demultiplex/TEST-input_FASTQ/r2.fq"
    r3 = "/projects/bgmp/oueslati/bioinfo/Bi622/AS1/Demultiplex/TEST-input_FASTQ/r3.fq"
    r4 = "/projects/bgmp/oueslati/bioinfo/Bi622/AS1/Demultiplex/TEST-input_FASTQ/r4.fq"
    
    barcodes = barcodes_to_dict(BARCODES_PATH, sep=' ')
    counts = store_barcode_pairs(barcodes, output=BARCODE_PAIRS_PATH)
    
    demultiplex(r1, r2, r3, r4, barcodes, counts)
    
    