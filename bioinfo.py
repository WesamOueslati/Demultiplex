#!/usr/bin/env python

# Author: WesamOueslati <optional@email.address>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.3"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = 'ATCGNatcgn'
RNA_bases = 'AUCGNaucgn'

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter) - 33

def qual_score(phred_score: str) -> float:
    '''returns the average quality score from a sequence'''
    sum = 0
    for letter in phred_score:
        sum+=convert_phred(letter)
    avg = sum/len(phred_score) 
    return avg


def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    seq = seq.upper()
    return len(seq) == seq.count("A") + seq.count("N") + seq.count("U" if RNAflag else "T") + seq.count("G") + seq.count("C")

def gc_content(DNA):
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(DNA)
    DNA = DNA.upper()
    return (DNA.count("G")+DNA.count("C"))/len(DNA)

def calc_median(lst: list):
    '''Given a sorted list, returns the median value of the list'''
    size = len(lst)
    if size % 2 == 0:
        return (lst[size // 2] + lst[(size // 2) - 1]) / 2
    return lst[size // 2]

def oneline_fasta(fasta_file):
    '''docstring'''
    output_file = f"oneline.fa"

    with open(fasta_file, 'r') as input:
            with open(output_file, 'w') as out:
                    lst = ""
                    for line in input:
                            if line.startswith('>'):
                                    if len(lst) != 0:
                                            out.write(f"{lst}\n")
                                            lst = ""
                                    out.write(line)
                            else:
                                    lst+= line.strip('\n')
                    out.write(lst)
    return output_file

def calc_sequenced_bases_fq(fq_file):
    num_sequenced_bases_fq = 0
    lengths = []
    i = 0
    with open(fq_file, 'r') as fq:
        for line in fq:
            if i % 4 == 1:
                line = line.strip('\n')
                lengths.append(len(line))
                num_sequenced_bases_fq += len(line)
            i+=1
    return num_sequenced_bases_fq, lengths

def calc_b_coverage(num_sequenced_bases, genome_size):
    """
    Calculates base coverage base on the following equation:
    C = (total sequenced bases) / (genome size)
    [PS6]
    """
    return (num_sequenced_bases) / (genome_size)

def calc_kmer_coverage(mean_read_length, kmer_size, base_coverage):
    """
    Calculates kmer coverage base on the following equation:
    Ck = C * (L - K + 1) / L
    [PS6]
    """    
    return base_coverage * (mean_read_length - kmer_size + 1) / mean_read_length



if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")
    
    
    assert gc_content("ACTTGC") == 0.5
    assert gc_content("CGGC") == 1
    assert gc_content("ATATATATAT") == 0
    
    assert calc_median([1,2,100]) == 2
    assert calc_median([1,2]) == 1.5
    assert calc_median([1,1,1,1,1,1,1,1,1,5000]) == 1
    assert calc_median([1,2,3,4,5,6,7,13]) == 4.5
    
    assert qual_score('IC') == 37
    assert qual_score('II') == 40    

    assert validate_base_seq('aTNcgN')
    assert validate_base_seq('not a dna seq') == False