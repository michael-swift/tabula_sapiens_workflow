from pprint import pprint
import os
import pandas as pd
import glob
from collections import defaultdict

test = False

shell.prefix("set +euo pipefail;")

configfile: "config/config.yaml"

outdir = config['out_dir']

os.makedirs(outdir, exist_ok=True)

# Resources
SIF="/oak/stanford/groups/quake/mswift/resources/immcantation_suite-4.1.0.sif"



# Data and Params
DIST="auto"
# Extract Sample Names from file paths 
FASTQ_DIR=config['raw_data']
fq_paths = glob.glob(FASTQ_DIR + '*5prime*')

samples = []
for path in fq_paths:
    fq = path.split('/')[-1]
    sample = fq.rsplit('_', 4)[0]
    samples.append(sample)
samples = list(set(samples))


# Testing
if test == True:
    samples = samples[:2]

rule all:
     input: expand("{outdir}/{sample}/outs/metrics_summary.csv", sample = samples, outdir = outdir)
     params: name="all", partition="normal"
     threads: 1

#### CellRanger VDJ #####
rule cellranger_vdj:
     input: config['raw_data']
     output: summary="{outdir}/{sample}/outs/metrics_summary.csv",
             fasta="{outdir}/{sample}/outs/filtered_contig.fasta", 
             anno="{outdir}/{sample}/outs/filtered_contig_annotations.csv"
     log: "{outdir}/log_{sample}/log.log"
     params: name="vdj_cellranger", 
             partition="quake",
             outdir=config['out_dir'],
             cell_ranger=config['cell_ranger'],
             ref=config['cell_ranger_ref']
     shell: "cd {params.outdir} && rm -rf {wildcards.sample} && {params.cell_ranger} vdj --id={wildcards.sample} --fastqs={input} --reference={params.ref} --sample={wildcards.sample}"
