import sys
import glob
import os
import pandas as pd
import uuid
""" Makes a samplesheet (and debugging seedfile with only 4 directories) for use with Derek's Snakemake pipeline
    must supply full directory path for optimal use(rather than relative paths)"""



base = sys.argv[1] # full directory path to the raw data directory
#name of samplesheet
file_name = sys.argv[2]
#gather all sub directories of the provided directory (non-recursivley)
directories = glob.glob(base +'*/')




# lot stands for list of tuples
lot = []

#generate lists that ultimately become columns of the dataframe

for path in directories:
    split = path.split('/')
    s = '/'
    tup = s.join(split[:-2]), split[-2], split[0], str(uuid.uuid4())
    lot.append(tup)              

#create dataframe
df = pd.DataFrame(lot, columns = ['base', 'samplename', 'group', 'unique'])

#write dataframe to .tsv file 
df.to_csv(file_name, sep='\t')

df.iloc[:8,:].to_csv(file_name+'.debug', sep = '\t')

