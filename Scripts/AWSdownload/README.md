## Download and check files from BMK
We get the sequence files delivered in Amazon S3 buckets and these
needs to be installed using
[aws-cli](https://aws.amazon.com/cli/). This tool is a bit cumbersome
to use and here we supply a shell script that automate downloading the
data and run checksums to make sure all files are intact after download.

The script assumes that aws-cli is installed and available in your
executable path. In addition you need to be on an operative system
that have `md5sum` available (most unix based system do have this
tool, but on bsd systems, including Mac os X, the tool is called `md5`
and the functionality might differ slightly).

## Howto use
Note: if you use this script several times, make sure to execute it
in different folders, otherwise the script will mix the md5sums from
both your current and previous download(s) (this would impact step 5).

1. Clone this directory including the scripts and the input.txt file
   to your local machine.
2. Copy the Index, access key id, and secret access key from the BMK
   delivery mail and save these three rows in the input.txt file.
   After this you should have a input.txt file looking as follows.
   
   ```
   Index in AWS: bmkdatarelease-5/delivery_xxxxxxxxxxxxxxxxx
                                            
   Access key ID: XXXXXXXXXXXXXXXXXXXX
                   
   Secret Access key: xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   ```
3. With the awsdownload.sh script and the input.txt in folder where
   you want data to be saved. Run the script like this (note that this
   will create a folder named "~/.aws" on your system containing your
   aws credentials):
   
   `bash awsdownload.sh`
   
5. Once the data has been downloaded you should see a list of files
   followed by ": OK" if everything is okay.
   ```
   Unknown_v1BMK230919-BP392-ZX01-060001-01_good_1.fq.gz: OK
   Unknown_v1BMK230919-BP392-ZX01-060001-01_good_2.fq.gz: OK
   Unknown_v2BMK230919-BP392-ZX01-060001-01_good_1.fq.gz: OK
   Unknown_v2BMK230919-BP392-ZX01-060001-01_good_2.fq.gz: OK
   ```

Any other type of outcome indicates issues with the download so make
sure to check downloaded files if you do not get an OK at the end.


