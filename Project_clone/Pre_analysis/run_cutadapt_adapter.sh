# if unsure about settings, run: python cutadapt_adapter.py --help
python3 cutadapt_adapter.py --primer_f AAACTCGTGCCAGCCACC --primer_r GGGTATCTAATCCCAGTTTG --sequenced_by bmk --linked no --reverse_flag no | tee cutadapt_log.out 
