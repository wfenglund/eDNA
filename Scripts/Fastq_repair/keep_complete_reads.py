# A python script that filters out intact reads from gzipped fastq files
import os
import gzip

# Edit which folder to look in or move script to desired folder:
file_folder = './'

# For every file in that folder:
for file in os.listdir(file_folder):
    if file.endswith('1.fq.gz'): # if it is a fastq forward file
        # Open input files:
        with gzip.open(file_folder + '/' + file) as forward_input:
            forward_text = forward_input.readlines()
        with gzip.open(file_folder + '/' + file.replace('1.fq.gz', '2.fq.gz')) as reverse_input:
            reverse_text = reverse_input.readlines()
        # Open output files:
        with gzip.open(file.replace('1.fq.gz', 'filtered_1.fq.gz'), 'w') as forward_output:
            with gzip.open(file.replace('1.fq.gz', 'filtered_2.fq.gz'), 'w') as reverse_output:
                # Go through every line in file:
                counter = 1 # start counter for line in read
                for line in forward_text:
                    if line.startswith(b'@') and counter == 1: # if a read starts
                        forward_output.write(line)
                        reverse_output.write(line)
                        counter = 2
                    elif counter == 2: # if it is the second line of a read
                        forward_output.write(line)
                        reverse_output.write(line)
                        counter = 3
                    elif line.startswith(b'+') and counter == 3: # if it is the third line of a read
                        forward_output.write(line)
                        reverse_output.write(line)
                        counter = 4
                    elif counter == 4: # if it is the fourth line of a read
                        forward_output.write(line)
                        reverse_output.write(line)
                        counter = 1 # reset line counter
