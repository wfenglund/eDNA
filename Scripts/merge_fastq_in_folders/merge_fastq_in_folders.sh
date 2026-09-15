# place script in folder with folders of fastas to be combined, run with:
# $ bash merge_fastq_in_folders.sh
for i in barcode*
do
  zcat ./$i/* > $i.fastq
  gzip $i.fastq
done
