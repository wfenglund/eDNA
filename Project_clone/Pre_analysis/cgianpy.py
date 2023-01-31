import os
import gzip
from operator import itemgetter

def get_unique_dicts(file_folder = '../Raw_data', line_length = 20, sample_lines = 100000):
    unique_dict_forward = {}
    unique_dict_reverse = {}
    for file in os.listdir(file_folder):
        if file.endswith('.gz'):
            with gzip.open(file_folder + '/' + file) as current_fastq:
                counter = 0
                for line in current_fastq:
                    line = line.decode('utf-8')
                    line_start = line[0:line_length]
                    if (counter == 1) or ((counter - 1) % 4 == 0): # if the position is a sequence position
                        if file.endswith('_1.fq.gz') or file.endswith('R1_001.fastq.gz'):
                            if line_start not in unique_dict_forward.keys():
                                unique_dict_forward[line_start] = 0
                            unique_dict_forward[line_start] = unique_dict_forward[line_start] + 1
                        if file.endswith('_2.fq.gz') or file.endswith('R2_001.fastq.gz'):
                            if line_start not in unique_dict_reverse.keys():
                                unique_dict_reverse[line_start] = 0
                            unique_dict_reverse[line_start] = unique_dict_reverse[line_start] + 1
                    counter = counter + 1
                    if counter >= sample_lines * 4:
                        break
    return unique_dict_forward, unique_dict_reverse

def print_dict(dict, primer, line_length = 20, table_length = 10):
    print('Sequence' + ' '*(line_length - 6) + 'Count')
    counter = 0
    for key, value in sorted(dict.items(), key=itemgetter(1), reverse = True):
        overlap = ''
        for nchars in range(1, len(key) + 1):
            if primer[:nchars] in key:
                overlap = primer[:nchars]
        sequence = key.replace(overlap, '\033[4m' + overlap + '\033[0m')
        print(f'{sequence}: {value}')
        counter = counter + 1
        if counter >= table_length:
            break

def rev_comp(sequence): # Function that takes a sequence (a string), upper or lower case, and reverse-complements it
    sequence = sequence[::-1] # Reverse the sequence
    new_seq = "" # Initialize the new sequence
    for nucleotide in sequence: # For every nucleotide/character in the sequence
        if nucleotide == "A" or nucleotide == "a":
            nucleotide = nucleotide.replace("A","T")
            nucleotide = nucleotide.replace("a","t")
            new_seq = new_seq + nucleotide
            continue
        elif nucleotide == "C" or nucleotide == "c":
            nucleotide = nucleotide.replace("C","G")
            nucleotide = nucleotide.replace("c","g")
            new_seq = new_seq + nucleotide
            continue
        elif nucleotide == "G" or nucleotide == "g":
            nucleotide = nucleotide.replace("G","C")
            nucleotide = nucleotide.replace("g","c")
            new_seq = new_seq + nucleotide
            continue
        elif nucleotide == "T" or nucleotide == "t" or nucleotide == "U" or nucleotide == "u":
            nucleotide = nucleotide.replace("T","A")
            nucleotide = nucleotide.replace("t","a")
            nucleotide = nucleotide.replace("U","A")
            nucleotide = nucleotide.replace("u","a")
            new_seq = new_seq + nucleotide
            continue
        elif nucleotide == "M" or nucleotide == "m":
            nucleotide = nucleotide.replace("M","K")
            nucleotide = nucleotide.replace("m","k")
            new_seq = new_seq + nucleotide
            continue
        elif nucleotide == "R" or nucleotide == "r":
            nucleotide = nucleotide.replace("R","Y")
            nucleotide = nucleotide.replace("r","y")
            new_seq = new_seq + nucleotide
            continue
        elif nucleotide == "Y" or nucleotide == "y":
            nucleotide = nucleotide.replace("Y","R")
            nucleotide = nucleotide.replace("y","r")
            new_seq = new_seq + nucleotide
            continue
        elif nucleotide == "K" or nucleotide == "k":
            nucleotide = nucleotide.replace("K","M")
            nucleotide = nucleotide.replace("k","m")
            new_seq = new_seq + nucleotide
            continue
        elif nucleotide == "V" or nucleotide == "v":
            nucleotide = nucleotide.replace("V","B")
            nucleotide = nucleotide.replace("v","b")
            new_seq = new_seq + nucleotide
            continue
        elif nucleotide == "H" or nucleotide == "h":
            nucleotide = nucleotide.replace("H","D")
            nucleotide = nucleotide.replace("h","d")
            new_seq = new_seq + nucleotide
            continue
        elif nucleotide == "D" or nucleotide == "d":
            nucleotide = nucleotide.replace("D","H")
            nucleotide = nucleotide.replace("d","h")
            new_seq = new_seq + nucleotide
            continue
        elif nucleotide == "B" or nucleotide == "b":
            nucleotide = nucleotide.replace("B","V")
            nucleotide = nucleotide.replace("b","v")
            new_seq = new_seq + nucleotide
            continue
        else:
            new_seq = new_seq + nucleotide
            continue
    return new_seq # Return the new sequence
