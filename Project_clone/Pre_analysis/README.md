This folder contains pre-analysis tools.

###### QC
Run:
```
$ python qc_generator.py
```
To generate a QC-report based on fastq-files in the "../Raw_data"-folder.
**Requirements:** FastQC, MultiQC (installed via pip or conda).

###### Trimming
Run:
```
$ python cutadapt_adapter.py
```
To trim reads in fastq-files in the "../Data"-folder. Edit script in a text-editor to specify settings (primers, file format, etc.).
**Requirements:** Cutadapt.
