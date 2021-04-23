from pprint import pprint
import os
import pandas as pd
import glob
from collections import defaultdict

shell.prefix("set +euo pipefail;")


configfile: "config/config.yaml"


outdir = config["out_dir"]
# should be loaded as config file
donors = [
    "pilot_10",
    "pilot_13",
    "pilot_4",
    "pilot_8",
    "pilot_9",
    "pilot_6",
    "pilot_3",
    "pilot_12",
]


os.makedirs(outdir, exist_ok=True)

# Extract Sample Names from file paths
FASTQ_DIR = config["raw_data"]
fq_paths = glob.glob(FASTQ_DIR + "*5prime*")

libs = []

for path in fq_paths:
    fq = path.split("/")[-1]
    lib = fq.rsplit("_", 4)[0]
    libs.append(lib)
    libs = list(set(libs))


def get_changeos(wildcards):
    parentdir = wildcards["outdir"]
    donor = wildcards["donor"]
    changeo = "{}/SS2/{}/bracer/filtered_BCR_summary/changeodb.tab".format(
        parentdir, donor
    )
    return [changeo]


test = False
# Testing
if test == True:
    libs = libs[:2]


include: "rules/get_containers.smk"
include: "rules/vdj.smk"

rule all:
    input:
        "{}/combined_igblast.airr.tsv".format(outdir),
    params:
        name="all",
        partition="normal",
    threads: 1
