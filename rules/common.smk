"""
Snakefile for basic zcat operations on gzipped fastqs and generation of directories.
Also includes paths to common scripts
"""
from glob import glob
import subprocess
import socket  # used instead of os for hostname


current_host = socket.gethostbyaddr(socket.gethostname())[0]

# select appropriate path based on current cluster
if "sh" in current_host:
    SCRIPTS = "/oak/stanford/groups/quake/mswift/B_cells/singleBcell_RNA_seq/scripts"
    REFS = "/oak/stanford/groups/quake/mswift/resources/"
    scratch_storage = "LOCAL_SCRATCH"

    add_col_htseq_script = SCRIPTS + "/alignment/add_col_htseq.py"


# -----------------------------------
# Function for loading seeds


def load_seeds(infile, seeds=[]):
    with open(infile, "rU") as f:
        for line in f:
            seeds.append(line.rstrip())
    return seeds


def load_seeds_with_basenames(infile, seeds=[], basenames=[]):
    with open(infile, "rU") as f:
        for line in f:
            vals = line.rstrip().split("\t")
            seeds.append(vals[0])
            basenames.append(vals[1])
    return seeds, basenames


def load_seeds_with_or_without_basenames(infile, directory=""):
    """ Detects whether the seedfile has 1 or 2 columns:
        For 2 columns:
            'base' is the first column (library e.g. PASL_1)
        For 1 column:
    """
    samples = []
    base = None
    with open(infile) as f:
        for line in f:
            sp = line.strip().split("\t")
            if len(sp) == 2:
                # if 2 columns, define base as list of directories
                if not base:
                    base = []
                base.append(os.path.join(directory, sp[0]))
                samples.append(sp[1])
            else:
                samples.append(sp[0])
    if base:
        pass
    else:
        # need to pad to len(samples) because using 'zip' in Snakefile expansion
        base = [directory] * len(samples)

    return base, samples


# -----------------------------------
# Functions for manipulating multiple fastq files


def get_all_files(d):
    return [d + "/" + f for f in os.listdir(d) if os.path.isfile(os.path.join(d, f))]


def get_fastq_gzs(wildcards):
    if "base" in wildcards.keys() and "sample" in wildcards.keys():
        prefix = "{}/{}".format(wildcards["base"], wildcards["sample"])
    print(prefix)
    r1 = glob("{}/*R1*.fastq.gz".format(prefix), recursive=True)
    r2 = glob("{}/*R2*.fastq.gz".format(prefix), recursive=True)
    try:
        assert len(r1) == 1
        assert len(r2) == 1
    except AssertionError:
        print("more than one read (maybe dealing with hiseq data)")
        print(r1)
        print(r2)
        pass
    return [r1[0], r2[0]]


def get_R1_R2_assembly(wildcards):
    parent_dir = wildcards.base
    sample_dir = wildcards.sample
    R1_R2_list = [
        glob(parent_dir + "/" + sample_dir + "/*R1*.fastq.gz")[0],
        glob(parent_dir + "/" + sample_dir + "/*R2*.fastq.gz")[0],
    ]
    return R1_R2_list

def get_all_fastq_gzs_R1_v3(wildcards):

    all_files = get_all_files(wildcards.base + "/" + wildcards.sample)
    fastq_gzs = []
    for f in all_files:
        if (
            ("00" in f)
            and ("_R1_" in f)
            and (".fastq" in f)
            and (os.path.splitext(f)[1] == ".gz")
        ):
            fastq_gzs.append(f)
    return sorted(fastq_gzs)


def get_all_fastq_gzs_R2_v3(wildcards):

    all_files = get_all_files(wildcards.base + "/" + wildcards.sample)
    fastq_gzs = []
    for f in all_files:
        if (
            ("00" in f)
            and ("_R2_" in f)
            and (".fastq" in f)
            and (os.path.splitext(f)[1] == ".gz")
        ):
            fastq_gzs.append(f)
    return sorted(fastq_gzs)


def unzip_fastq(f):
    cmd = "gunzip " + f
    p = subprocess.Popen(cmd, shell=True)
    return None


# -----------------------------------
# Common rules


rule sort_and_index_bam:
    input:
        "{infile}.{bam}",
    output:
        bam=temp("{infile}.sorted.{bam}"),
        ind=temp("{infile}.sorted.{bam}.bai"),
    threads: 4
    params:
        name="sort_bam",
        partition=config["default_partition"],
    resources:
        mem_mb=5000,
    conda:
        os.path.join(workflow.basedir, "envs/singlecell.yaml")
    shell:
        "samtools sort --threads {threads} -o {output.bam} {input} && "
        "samtools index {output.bam}"


rule cat_R1_v3:
    """ Concatenate fastq files """
    input:
        get_fastq_gzs,
    output:
        "{base}/{sample}/R1.fastq.gz",
    params:
        name="zcat",
        partition="normal",
        time="0:1",
    resources:
        mem_mb=5300,
    # wildcard_constraints: sample="^\d+"  # avoid interference with subsampling
    shell:
        "mv {input[0]} {output}"


rule cat_R2_v3:
    input:
        get_fastq_gzs,
    output:
        "{base}/{sample}/R2.fastq.gz",
    params:
        name="zcat",
        partition="normal",
        time="0:1",
    resources:
        mem_mb=5300,
    # wildcard_constraints: sample="^\d+"  # avoid interference with subsampling
    shell:
        "mv {input[1]} {output}"
