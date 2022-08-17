import os
import gzip
import re
from operator import itemgetter

def get_unique_dicts(file_folder = '../Data', line_length = 20, sample_lines = 100000):
    unique_dict_forward = {}
    unique_dict_reverse = {}
    for file in os.listdir(file_folder):
        if file.endswith('.gz'):
            with gzip.open(file_folder + '/' + file) as current_fastq:
                counter = 0
                for line in current_fastq:
                    line = line.decode('utf-8')
                    line_start = line[0:line_length]
                    if re.match("^[ATGC]{5}", line_start):
                        if file.endswith('_1.fq.gz') or file.endswith('R1_001.fastq.gz'):
                            if line_start not in unique_dict_forward.keys():
                                unique_dict_forward[line_start] = 0
                            unique_dict_forward[line_start] = unique_dict_forward[line_start] + 1
                        if file.endswith('_2.fq.gz') or file.endswith('R2_001.fastq.gz'):
                            if line_start not in unique_dict_reverse.keys():
                                unique_dict_reverse[line_start] = 0
                            unique_dict_reverse[line_start] = unique_dict_reverse[line_start] + 1
                        counter = counter + 1
                    if counter >= sample_lines:
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
