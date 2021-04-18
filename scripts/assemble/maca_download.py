import os
import boto3
import sys

base, sample, read = sys.argv[1], sys.argv[2], sys.argv[3]

s3 = boto3.resource('s3')
bucket = s3.Bucket('czbiohub-maca')

fld = os.path.join(base, sample)

if not os.path.exists(base):
    os.mkdir(base)

if not os.path.exists(fld):
    os.mkdir(fld)


bucket_file = '%s_%s_001.fastq.gz' % (sample, read)
print(bucket_file)

bucket.download_file(bucket_file, os.path.join(fld, '%s.fastq.gz' % read))
