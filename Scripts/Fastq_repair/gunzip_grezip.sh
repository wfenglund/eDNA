# A bash script that gunzips and re-gzips fastq-files in current folder
# (use if there is something wrong with the gzipping)

for zipped_file in ./*.gz
do
	echo "Unzipping $zipped_file"
	gunzip $zipped_file
	unzipped_file=${zipped_file/".gz"/""}
	echo "Rezipping $unzipped_file"
	gzip $unzipped_file
done
