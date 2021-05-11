import sys
import os

# -----------------------------
# Configuration

species = config["species"]
# Cell Ranger
#### CellRanger VDJ #####


rule cellranger_vdj:
    input:
        fastqs=config["raw_data"],
        ref=rules.get_imgt_db.output.ref,
    output:
        "{base}/10X/{lib}/outs/metrics_summary.csv",
        "{base}/10X/{lib}/outs/filtered_contig.fasta",
        "{base}/10X/{lib}/outs/filtered_contig_annotations.csv",
        "{base}/10X/{lib}/outs/airr_rearrangement.tsv",
    log:
        "{base}/logs/cellranger_{lib}.log",
    params:
        name="vdj_cellranger",
        base=config["base"],
        cell_ranger=config["cell_ranger"],
    shell:
        "cd {params.base}/10X && rm -rf {wildcards.lib} && {params.cell_ranger}/cellranger vdj --id={wildcards.lib} --fastqs={input.fastqs} --reference={input.ref} --sample={wildcards.lib}"

rule igblast_10X:
    input:
        "{base}/10X/{lib}/outs/filtered_contig.fasta",
    output:
        "{base}/10X/igblast/{lib}_igblast.airr.tsv",
    conda:
        os.path.join(workflow.basedir, "envs/vdj.yaml")
    params:
        organism="human",
        IGDBDIR=config["IGDBDIR"],
        seqtype =sense_lib_type
    shell:
        """
        ml load system libuv
        wdir=$(dirname {input[0]})
        export IGDATA={params.IGDBDIR}
        igblastn \
        -germline_db_V {params.IGDBDIR}/database/imgt_{params.organism}_{params.seqtype[1]}_v \
        -germline_db_D {params.IGDBDIR}/database/imgt_{params.organism}_{params.seqtype[1]}_d \
        -germline_db_J {params.IGDBDIR}/database/imgt_{params.organism}_{params.seqtype[1]}_j \
        -auxiliary_data {params.IGDBDIR}/optional_file/{params.organism}_gl.aux \
        -domain_system imgt \
        -ig_seqtype {params.seqtype[0]} \
        -organism {params.organism} \
        -outfmt 19 \
        -query {input} \
        -out {output}
        """

rule edit_10X_igblast:
    input:
        tsv="{base}/10X/igblast/{lib}_igblast.airr.tsv",
        airr_10X="{base}/10X/{lib}/outs/airr_rearrangement.tsv",
    output:
        tsv="{base}/10X/igblast/{lib}_igblast_edit.airr.tsv",
    run:
        df = pd.read_table(input.tsv, sep="\t")
        df.loc[:,"sequence_id"] = df.sequence_id + "_" + wildcards.lib
        df.loc[:,"library"] = "10X_vdj"
        df.to_csv(output.tsv, sep="\t", index=False, header=True)


### get_bracer_contigs:

rule get_bracer_contigs:
    input:
        get_changeos,
    output:
        fasta="{base}/SS2/{donor}/contigs.fasta",
    conda:
        os.path.join(workflow.basedir, "envs/bcr.yaml")
    log:
        "{base}/logs/{donor}/get_fastas.log",
    params:
        name="extract_seqs",
        partition="quake,owners",
        scripts=os.path.join(workflow.basedir, "scripts"),
    shell:
        "python {params.scripts}/get_bracer_contigs.py {input} {output}"

# TODO run this with biopython
rule combine_bracer_contigs:
    input:
        expand("{base}/SS2/{donor}/contigs.fasta", base=base, donor=donors),
    output:
        "{base}/SS2/BCR_combined_contigs.fasta",
    conda:
        os.path.join(workflow.basedir, "envs/bcr.yaml")
    log:
        "{base}/logs/concat_fastas.log",
    params:
        name="combine_seqs",
        partition="quake,owners",
    shell:
        "cat {input} > {output}"


rule igblast_bracer:
    input:
        rules.combine_bracer_contigs.output,
    output:
        "{base}/SS2/igblast/bracer.airr.tsv",
    conda:
        os.path.join(workflow.basedir, "envs/vdj.yaml")
    params:
        organism="human",
        IGDBDIR=config["IGDBDIR"],
        seqtype="Ig",
    shell:
        """
        ml load system libuv
        wdir=$(dirname {input[0]})
        export IGDATA={params.IGDBDIR}
        igblastn \
        -germline_db_V {params.IGDBDIR}/database/imgt_{params.organism}_ig_v \
        -germline_db_D {params.IGDBDIR}/database/imgt_{params.organism}_ig_d \
        -germline_db_J {params.IGDBDIR}/database/imgt_{params.organism}_ig_j \
        -auxiliary_data {params.IGDBDIR}/optional_file/{params.organism}_gl.aux \
        -domain_system imgt \
        -ig_seqtype {params.seqtype} \
        -organism {params.organism} \
        -outfmt 19 \
        -query {input} \
        -out {output}
        """


rule edit_bracer_igblast:
    input:
        tsv="{base}/SS2/igblast/bracer.airr.tsv",
    output:
        tsv="{base}/SS2/igblast/bracer_lib.airr.tsv",
    run:
        df = pd.read_table(input.tsv, sep="\t")
        df['library'] = 'bracer'
        df.to_csv(output.tsv, sep="\t", index=False, header=True)


# This is currently so fucked. should specify files explicity
rule get_tracer_contigs:
    input:
        di="{base}/SS2/",
    output:
        fasta="{base}/SS2/tracer_contigs.fasta",
    log:
        "{base}/logs/get_tracer_fastas.log",
    params:
        name="extract_tcr_seqs",
        partition="quake,owners",
    run:
        from Bio import SeqIO
        fastas = glob.glob(input.di + "/*/tracer/assembled/*/filtered*/*.fa*")
        records = []
        for fasta in fastas:
            cellname = fasta.split("/")[-3]
            donor = fasta.split("/")[-6]
            for record in SeqIO.parse(fasta, "fasta"):
                record.id = "{}|{}|{}".format(
                    cellname, record.description, donor)
                records.append(record)
        SeqIO.write(records, output.fasta, "fasta")


rule igblast_tracer:
    input:
        fasta="{base}/SS2/tracer_contigs.fasta",
        sif="{}/resources/repertoire/{}.sif".format(workflow.basedir, "immcantation"),
    output:
        "{base}/SS2/igblast/tracer.airr.tsv",
    conda:
        os.path.join(workflow.basedir, "envs/vdj.yaml")
    params:
        organism="human",
        IGDBDIR=config["IGDBDIR"],
        seqtype="TCR",
        dbtype = "tr",
    shell:
        """
        ml load system libuv
        wdir=$(dirname {input[0]})
        export IGDATA={params.IGDBDIR}
        igblastn \
        -germline_db_V {params.IGDBDIR}/database/imgt_{params.organism}_{params.dbtype}_v \
        -germline_db_D {params.IGDBDIR}/database/imgt_{params.organism}_{params.dbtype}_d \
        -germline_db_J {params.IGDBDIR}/database/imgt_{params.organism}_{params.dbtype}_j \
        -auxiliary_data {params.IGDBDIR}/optional_file/{params.organism}_gl.aux \
        -domain_system imgt \
        -ig_seqtype {params.seqtype} \
        -organism {params.organism} \
        -outfmt 19 \
        -query {input.fasta} \
        -out {output}
        """
rule edit_tracer_igblast:
    input:
        tsv="{base}/SS2/igblast/tracer.airr.tsv",
    output:
        tsv="{base}/SS2/igblast/tracer_lib.airr.tsv",
    run:
        df = pd.read_table(input.tsv, sep="\t")
        df['library'] = 'tracer'
        df.to_csv(output.tsv, sep="\t", index=False, header=True)

## Combined Outputs
rule combine_igblast:
    input:
        TenXs=expand(
            "{base}/10X/igblast/{lib}_igblast_edit.airr.tsv", lib=libs, base=base
        ),
        bracer="{base}/SS2/igblast/bracer_lib.airr.tsv",
        tracer="{base}/SS2/igblast/tracer_lib.airr.tsv",
    output:
        tsv="{base}/vdj/combined_igblast.airr.tsv",
    log:
        "{base}/logs/combineigblast.log",
    run:
        dfs = []
        infiles = input.TenXs
        infiles.append(input.bracer)
        infiles.append(input.tracer)
        for i in infiles:
            df = pd.read_table(i, sep="\t")
            dfs.append(df)

        combined = pd.concat(dfs)
        combined.to_csv(output.tsv, sep="\t", index=False, header=True)

rule split_loci:
    input:db="{base}/vdj/combined_igblast.airr.tsv"
    output:bcr="{base}/vdj/ig_airr.tsv",tcr="{base}/vdj/tr_airr.tsv"
    log: "{base}/logs/split.log"
    run:
        df = pd.read_table(input.db, sep = "\t")
        df.dropna(subset=['locus'], inplace = True)
        df_out = df[df.locus.str.contains('IG')]
        df_out.to_csv(output.bcr, index=False, header = True, sep = "\t")
        df_out = df[~df.locus.str.contains('IG')]
        df_out.to_csv(output.tcr, index=False, header = True, sep = "\t")

rule changeo_clone:
    input:
        bcr_db="{base}/vdj/ig_airr.tsv", sif=rules.get_immcantation_image.output
    output:
        "{base}/vdj/changeo/combined_germ-pass.tsv"
    conda:
        "../envs/vdj.yaml"
    params:
        dist="0.15",
        sample_name="combined",
        nproc= '2'
    log:
        "{base}/logs/changeo_clone.log",
    shell:
        "singularity exec -B {wildcards.base}:/data {input.sif} changeo-clone -x {params.dist} -d {input.bcr_db} -n {params.sample_name} -o /data/vdj/changeo -p {params.nproc}" 

rule annotate_constant_region:
    input:
        "{base}/vdj/changeo/combined_germ-pass.tsv"
    output:
        "{base}/vdjc/combined_vdjc.tsv.gz",
    conda:
        "../envs/vdj.yaml"
    params:
        ighc_db=config["ighc_db"]["human"],
        scripts=config["scripts"],
    log:
        "{base}/logs/annotate_constant_region.log",
    shell:
        "cat {input} > {output}"
"""
        "python {params.scripts}/blast_constant_region.py "
        "{input} "
        "--min_j_sequence_length 15 "
        "--allow_missing_cdr3 True "
        "--allow_ns_in_sequence True "
        "-ighc_db {params.ighc_db} "
        "-outdir {wildcards.base}/vdjc "
        "> {log}"
"""
