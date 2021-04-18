import os
import pandas as pd


include: 'common.smk'

rule fastp:
    """ fastp does what fastqc does but better afaik"""
    input:
        #get_r1_r2_fastq_gz
        #get_fastq_gz_MS #change this for R1 R2 vs. older hiseq output
        rules.cat_R1_v3.output, rules.cat_R2_v3.output
    output:
          fastp_output1='{base}/{sample}/qc/fastp/read1.fastq.gz',
          fastp_output2='{base}/{sample}/qc/fastp/read2.fastq.gz',
          failed_out='{base}/{sample}/qc/fastp/failed.out'
    params:
          name='fastp',
          partition="normal,owners,quake"
    resources:
          mem_mb=5000
    conda:
          os.path.join(workflow.basedir, 'envs/singlecell.yaml')
    shell:
        "celldir=$(dirname {output[0]}) &&"
        " cd $celldir &&"
        " fastp --in1 {input[0]} --in2 {input[1]}"
        " --out1 {output[0]} --out2 {output[1]}"
        " --failed_out {output[2]}"
        " --disable_adapter_trimming"

rule fastqc:
    """ fastqc does not let you specify an output name thus
        moving files to standardize
    """ 
    input:
        rules.cat_R1_v3.output
    output:
        fastqc_html='{base}/{sample}/qc/{R}_fastqc.html',
        fastqc_zip='{base}/{sample}/qc/{R}_fastqc.zip'
    threads: 2
    params:
        name="fastqc",
        partition="normal,owners",
        no_ext=lambda wildcards, input: str(input).strip('.gz').strip('.fastq')
    resources:
        mem_mb=5000
    conda:
        os.path.join(workflow.basedir, 'envs/singlecell.yaml')
    shell:
        "mkdir -p $(dirname {output.fastqc_html}) &&"
        " fastqc {input} &&"
        " mv {params.no_ext}_fastqc.html {output.fastqc_html} &&"
        " mv {params.no_ext}_fastqc.zip {output.fastqc_zip}"
    

rule multiqc:
    input:
        expand('{base}/{sample}/qc/{R}_fastqc.html',
               zip,
               base=samplesheet.base.values.tolist(),
               sample=samplesheet.samplename.values.tolist(),
               R=['R1']*len(samplesheet.base.values.tolist()) + ['R2']*len(samplesheet.base.values.tolist())
              )
    output:
        '{processed_data}/multiqc_report.html'
    params:
        name="multiqc",
        partition="quake,owners",
        time="4-0",
        # multiqc searches one directory, so find common prefix to all base directories
        commonbasedir=os.path.dirname(
            os.path.commonprefix(
                (pd.Series(samplesheet.base.values.tolist())+'/').unique().tolist()))
    resources:
        mem_mb=20000
    conda:
        os.path.join(workflow.basedir, 'envs/singlecell.yaml')
    shell:
        "multiqc -d -o $(dirname {output}) {params.commonbasedir}"

