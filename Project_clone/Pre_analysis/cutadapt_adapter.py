import cgianpy
import os
import subprocess

###Inputs:
project_name = "proj_test"
primer_f = "AAACTCGTGCCAGCCACC" # Forward primer (default = MiFish)
primer_r = "GGGTATCTAATCCCAGTTTG" # Reverse primer (default = MiFish)
###
f_format = "novo" #scilife or novo
anchored_flag = "yes" #yes or no
reverse_flag = "yes" #yes or no
reads_folder = "../Raw_data" # Location of raw fastq-file
write_folder = "../Filtered_data" # Location where trimmed reads will be written

###Main:

# parameters are set
test1 = 0
test2 = 0
if anchored_flag == "no":
    primer_f_rc = ""
    primer_r_rc = ""
    regular_primer = primer_f
    reverse_primer = primer_r
    test1 = 1
elif anchored_flag == "yes":
    primer_f_rc = cgianpy.rev_comp(primer_f)
    primer_r_rc = cgianpy.rev_comp(primer_r)
    regular_primer = primer_f + "..." + primer_r_rc
    reverse_primer = primer_r + "..." + primer_f_rc
    test1 = 1
if f_format == "scilife":
    file_end = "_R1_001.fastq.gz"
    repl_from = "R1"
    repl_to = "R2"
    test2 = 1
elif f_format == "novo":
    file_end = "_1.fq.gz"
    repl_from = "_1.fq"
    repl_to = "_2.fq"
    test2 = 1
if test1 + test2 != 2:
    print("Please revisit your input variables. Program terminates.")
    quit()

# cutadapt is run on each fastq-file in the data-folder
for file in os.listdir(reads_folder): #For every file in folder containing reads
    if file.endswith(file_end): #Check if file is a .fastq-read-file
        sample = file.replace(file_end, "") #Extract sample name
        current_list = list([write_folder + "/" + sample + "_outFwd_1.fastq.gz", write_folder + "/" + sample + "_outFwd_2.fastq.gz", write_folder + "/" + sample + "_outRev_1.fastq.gz", write_folder + "/" + sample + "_outRev_2.fastq.gz"]) #Create filtering names and store in a list
        
        subprocess.run([f'cd {reads_folder} ; cutadapt -j 0 --max-n=0 --discard-untrimmed -g {regular_primer} -G {reverse_primer} -o {current_list[0]} -p {current_list[1]} {file} {file.replace(repl_from, repl_to)}'], shell=True) #Run cutadapt for regular direction
        if reverse_flag == "yes":
            subprocess.run([f'cd {reads_folder} ; cutadapt -j 0 --max-n=0 --discard-untrimmed -g {regular_primer} -G {reverse_primer} -o {current_list[2]} -p {current_list[3]} {file.replace(repl_from, repl_to)} {file}'], shell=True) #Run cutadapt in reverse direction
