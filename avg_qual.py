#!/usr/bin/env python
import bioinfo
import gzip

def calc_avg_qual(file_path, seq_length, output):
    
    with gzip.open(file_path, 'rt', encoding="utf-8") as fh:
        seq_length = int(seq_length)
        i = 0
        qual_scores = [0.0] * seq_length
        
        for line in fh:
            if i % 4 == 3:
                line = line.strip('\n')
                for j, char in enumerate(line):
                    qual_val = bioinfo.convert_phred(char)
                    qual_scores[j] += qual_val
            i+=1
                        
        avg_quals = [(score / (i//4)) for score in qual_scores]
        with open(output, 'w') as out:
            for base in avg_quals:
                out.write(str(f"{base}\n"))
    return avg_quals


if __name__ == "__main__":
    
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help="input file")
    parser.add_argument('-r', '--read_length', help="Sequence length for the reads")
    parser.add_argument('-o', '--output', help="output file")

    args = parser.parse_args()
    input = args.input
    read_length = args.read_length
    output = args.output
    
    avgs = calc_avg_qual(input, read_length, output)
    