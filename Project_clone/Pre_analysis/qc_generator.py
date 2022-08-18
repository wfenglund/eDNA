import cgianpy
import os
import subprocess
import multiqc

###Inputs:
forward_primer = 'AAACTCGTGCCAGCCACC'

###Main:

# tables containing the most frequent reads in forward and reverse reads are printed
print('Collecting reads...\n')
unique_dict_forward, unique_dict_reverse = cgianpy.get_unique_dicts()
cgianpy.print_dict(unique_dict_forward, forward_primer) # a table with forward-sequences is printed
print()
cgianpy.print_dict(unique_dict_reverse, forward_primer) # a table with reverse-sequences is printed

# fastqc is run on every .gz (fastq) file in the Raw_data directory
gz_files = [file for file in os.listdir('../Raw_data/') if file.endswith('.gz')]
print(f'\nRunning FastQC on {len(gz_files)} files:')
for file in gz_files:
    print(f'File {(gz_files.index(file) + 1)} being processed...')
    subprocess.run(['fastqc ../Raw_data/' + file + ' --outdir=./ --quiet'], shell = True)
print('Done.')

# Multiqc is run in current directory to summarize the fastqc-logs
print('\nRunning MultiQC...')
multiqc.run('.')
