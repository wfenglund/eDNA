import cgianpy
import os
import subprocess
import argparse

parser = argparse.ArgumentParser(prog = "cutadapt_adapter", description='Filter sequence reads using cutadapt')

parser.add_argument('-p', '--project_name',
                    default = 'proj_test',
                    help = "Usually the journal number (diarienummer)")
parser.add_argument('-f', '--primer_f',
                    default = 'AAACTCGTGCCAGCCACC',
                    help = "Sequence of the forward primer used to create the amplicon")
parser.add_argument('-r', '--primer_r',
                    default = 'GGGTATCTAATCCCAGTTTG',
                    help = "Sequence of the reverse primer used to create the amplicon")
parser.add_argument('-s', '--sequenced_by',
                    default = 'bmk', choices = ['bmk', 'novo', 'scilife', 'minion'],
                    help = "Sequence center that generated the data")
parser.add_argument('-a', '--linked',
                    default = 'no', choices = ['yes', 'no'],
                    help = "Is the amplicon shorter than the read length")
parser.add_argument('-b', '--reverse_flag',
                    default = 'no', choices = ['yes', 'no'],
                    help = "Is the sequence library generated directional")
parser.add_argument('-R', '--reads_folder',
                    default = '../Raw_data', choices = ['../Raw_data', '../Unfiltered_data'],
                    help = "Path with the read data to be filtered")
parser.add_argument('-W', '--write_folder',
                    default = '../Filtered_data',
                    help = "Folder to write the filtered data to (must exist)")
parser.add_argument('-O', '--out_prefix',
                    default = '',
                    help = "Prefix to be added before every output file name")
args = parser.parse_args()

###Inputs:
project_name = args.project_name
primer_f = args.primer_f
primer_r = args.primer_r
f_format = args.sequenced_by
linked_flag = args.linked
reverse_flag = args.reverse_flag
reads_folder = args.reads_folder
write_folder = args.write_folder
out_prefix = args.out_prefix

###Main:

# parameters are set
test1 = 0
test2 = 0
if f_format == "scilife":
    file_end = "_R1_001.fastq.gz"
    repl_from = "R1"
    repl_to = "R2"
    paired_flag = 'yes'
    test2 = 1
elif f_format == "novo" or f_format == "bmk":
    file_end = "_1.fq.gz"
    repl_from = "_1.fq"
    repl_to = "_2.fq"
    paired_flag = 'yes'
    test2 = 1
elif f_format == 'minion':
    file_end = '.fastq.gz'
    repl_from = ''
    repl_to = ''
    linked_flag = 'yes' # since minion reads are always assumed to have both primers
    paired_flag = 'no'
    test2 = 1
if linked_flag == "no":
    primer_f_rc = ""
    primer_r_rc = ""
    regular_primer = primer_f
    reverse_primer = primer_r
    test1 = 1
elif linked_flag == "yes":
    primer_f_rc = cgianpy.rev_comp(primer_f)
    primer_r_rc = cgianpy.rev_comp(primer_r)
    regular_primer = primer_f + "..." + primer_r_rc
    reverse_primer = primer_r + "..." + primer_f_rc
    test1 = 1
if test1 + test2 != 2:
    print("Please revisit your input variables. Program terminates.")
    quit()

# cutadapt is run on each fastq-file in the data-folder
for file in os.listdir(reads_folder): #For every file in folder containing reads
    if file.endswith(file_end): #Check if file is a .fastq-read-file
        sample = file.replace(file_end, "") #Extract sample name
        current_list = list([write_folder + "/" + out_prefix + sample + "_outFwd_1.fastq.gz", write_folder + "/" + out_prefix + sample + "_outFwd_2.fastq.gz", write_folder + "/" + out_prefix + sample + "_outRev_1.fastq.gz", write_folder + "/" + out_prefix + sample + "_outRev_2.fastq.gz"]) #Create filtering names and store in a list
        
        if paired_flag == "yes":
            subprocess.run([f'cd {reads_folder} ; cutadapt -j 0 --max-n=0 --discard-untrimmed -g {regular_primer} -G {reverse_primer} -o {current_list[0]} -p {current_list[1]} {file} {file.replace(repl_from, repl_to)}'], shell=True) #Run cutadapt for regular direction
            if reverse_flag == "yes":
                subprocess.run([f'cd {reads_folder} ; cutadapt -j 0 --max-n=0 --discard-untrimmed -g {regular_primer} -G {reverse_primer} -o {current_list[2]} -p {current_list[3]} {file.replace(repl_from, repl_to)} {file}'], shell=True) #Run cutadapt in reverse direction
        else: # if input is single ended
            subprocess.run([f'cd {reads_folder} ; cutadapt -j 0 --max-n=0 --revcomp --discard-untrimmed -g {regular_primer} -o {current_list[0]} {file}'], shell=True) #Run cutadapt single ended

