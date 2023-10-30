#!/bin/bash

# parse input file with new credentials:
awsinfo=(`awk '{print $4}' input.txt`)

# File that stores aws credentials:
credentials_file=~/.aws/credentials

# Save current credentials to backup if present:
if test -f "$credentials_file"; then
    cat ~/.aws/credentials >> ~/.aws/credentials.bak
    #echo "vimrules"
fi

# Write new credentials:
echo "# "${awsinfo[0]} > $credentials_file
echo "[default]" >> $credentials_file
echo "aws_access_key_id = "${awsinfo[1]} >> $credentials_file
echo "aws_secret_access_key = "${awsinfo[2]} >> $credentials_file

# Download data from amazon S3:
aws s3 cp s3://${awsinfo[0]} . --recursive

# Check md5sums:
md5sum */*/Data/*.fq.gz > localmd5sum.txt
sort -k2 localmd5sum.txt > localmd5sum.txt
sort -k2 */*/Data/data_md5.txt > bmkmd5sum.txt
echo "Differences between md5sums (if any):"
diff localmd5sum.txt bmkmd5sum.txt
