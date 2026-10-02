# if unsure about settings, run: python cutadapt_adapter.py --help

# - Only mifish:
python3 cutadapt_adapter.py --primer_f AAACTCGTGCCAGCCACC --primer_r GGGTATCTAATCCCAGTTTG --sequenced_by bmk --linked no --reverse_flag no --out_prefix mifish | tee cutadapt_log.out 

# - mifish and v16s:
# python3 cutadapt_adapter.py --primer_f AAACTCGTGCCAGCCACC --primer_r GGGTATCTAATCCCAGTTTG --sequenced_by bmk --linked no --reverse_flag no --out_prefix mifish | tee cutadapt_log_mifish.out 
# python3 cutadapt_adapter.py --primer_f ACGAGAAGACCCYRYGRARCTT --primer_r TCTHRRRANAGGATTGCGCTGTTA --sequenced_by bmk --linked no --reverse_flag no --out_prefix v16s | tee cutadapt_log_v16s.out 

# - IBA insect primers (should just be concatenated in dada2):
# python3 cutadapt_adapter.py --primer_f CCHGAYATRGCHTTYCCHCG --primer_r TCDGGRTGNCCRAARAAYCA --sequenced_by bmk --linked no --reverse_flag no --out_prefix IBA | tee cutadapt_log_IBA.out 

# - EPT insect primers:
# python3 cutadapt_adapter.py --primer_f GGDACWGGWTGAACWGTWTAYCCHCC --primer_r CAAACAAATARDGGTATTCGDTY --sequenced_by bmk --linked no --reverse_flag yes --out_prefix EPT | tee cutadapt_log_EPT.out 

# - batra amphibian primers:
# python3 cutadapt_adapter.py --primer_f ACACCGCCCGTCACCCT --primer_r GTAYACTTACCATGTTACGACTT --sequenced_by bmk --linked no --reverse_flag no --out_prefix batra | tee cutadapt_log_batra.out

# - mussels, ITS and 16S (ITS should just be concatenated in dada2):
# python3 cutadapt_adapter.py --primer_f AGACTGGGTTGCGGAGGT --primer_r CGAGTGATCCACCGCTTAGA --sequenced_by bmk --linked no --reverse_flag no --out_prefix batra | tee cutadapt_log_mITS.out
# python3 cutadapt_adapter.py --primer_f GCTGTTATCCCCGGGGTAR --primer_r AAGACGAAAAGACCCCGC --sequenced_by bmk --linked no --reverse_flag no --out_prefix batra | tee cutadapt_log_m16S.out
