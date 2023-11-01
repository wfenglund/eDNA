## Download and check files from BMK
We get the sequence files delivered in Amazon S3 buckets. The easiest way to download data is to use [aws cli](https://aws.amazon.com/cli/).

In order to use the script downloadscript available here you need to have aws cli installed and available in your path. In addition you need to be on an operative system that have `md5sum` available (most unix based system do have this tool, but on bsd systems, including Mac os X, the tool is called `md5`.

