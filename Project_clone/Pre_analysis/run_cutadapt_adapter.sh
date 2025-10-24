# if unsure about settings, run: python cutadapt_adapter.py --help
python3 cutadapt_adapter.py --primer_f AAACTCGTGCCAGCCACC --primer_r GGGTATCTAATCCCAGTTTG --sequenced_by bmk --linked no --reverse_flag no --out_prefix mifish | tee cutadapt_log.out 

# python3 cutadapt_adapter.py --primer_f AAACTCGTGCCAGCCACC --primer_r GGGTATCTAATCCCAGTTTG --sequenced_by bmk --linked no --reverse_flag no --out_prefix mifish | tee cutadapt_log_mifish.out 
# python3 cutadapt_adapter.py --primer_f ACGAGAAGACCCYRYGRARCTT --primer_r TCTHRRRANAGGATTGCGCTGTTA --sequenced_by bmk --linked no --reverse_flag no --out_prefix v16s | tee cutadapt_log_v16s.out 
