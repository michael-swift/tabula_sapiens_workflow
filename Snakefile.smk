import os
import pandas as pd
import glob
from collections import defaultdict

shell.prefix("set +euo pipefail;")


configfile: "config/config.yaml"


base = config["base"]
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


os.makedirs(base, exist_ok=True)

# Extract Sample Names from file paths
def get_10X_fastqs(FASTQ_DIR):

    fq_paths = glob.glob(FASTQ_DIR + "*5prime*")
    libs = []
    for path in fq_paths:
        fq = path.split("/")[-1]
        lib = fq.rsplit("_", 4)[0]
        libs.append(lib)
    libs = list(set(libs))
    libs.sort()
    return libs


libs = get_10X_fastqs(config["raw_data"])
# Extract changeos from angelas data


def get_changeos(wildcards):
    parentdir = wildcards["base"]
    donor = wildcards["donor"]
    changeo = "{}/SS2/{}/bracer/filtered_BCR_summary/changeodb.tab".format(
        parentdir, donor
    )
    return [changeo]


# Sense the 10X library type
def sense_lib_type(wildcards):
    """defaults to Ig"""

    lib = wildcards["lib"]
    if "TCR" in lib:
        return "TCR"
    elif "BCR" in lib:
        return "Ig"
    else:
        "break IgBlast"


## Debug mode
test = False
# Testing
if test == True:
    libs = libs[:2]


include: "rules/get_resources.smk"
include: "rules/vdj.smk"


localrules:
    get_tracer_contigs,


rule all:
    input:
        "{}/vdjc/combined_vdjc.tsv.gz".format(base),
    params:
        name="all",
        partition="normal",
    threads: 1
