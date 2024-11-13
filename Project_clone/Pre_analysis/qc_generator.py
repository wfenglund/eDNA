import argparse
import cgianpy
import os
import subprocess
import multiqc

###Inputs:
parser = argparse.ArgumentParser(prog = "QC_Generator", description='Generate QC files and statistics from fastq-files.')
parser.add_argument('-f', '--primer_f',
                    default = 'AAACTCGTGCCAGCCACC',
                    help = "Sequence of the forward primer used to create the amplicon")
parser.add_argument('-d', '--data_folder',
                    default = '../Raw_data',
                    help = "Sequence of the forward primer used to create the amplicon")
parser.add_argument('-t', '--only_tables',
                    action = 'store_true',
                    help = "Option to run QC_Generator and only create the sequence tables")
args = parser.parse_args()

forward_primer = args.primer_f
data_folder = args.data_folder
only_tables = args.only_tables

###Main:

# tables containing the most frequent reads in forward and reverse reads are printed
print('Collecting reads...\n')
unique_dict_forward, unique_dict_reverse = cgianpy.get_unique_dicts(file_folder = data_folder)
cgianpy.print_dict(unique_dict_forward, forward_primer) # a table with forward-sequences is printed
print()
cgianpy.print_dict(unique_dict_reverse, forward_primer) # a table with reverse-sequences is printed

if only_tables == False: # if only_tables flag has not been used
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

