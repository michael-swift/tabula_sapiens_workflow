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
        seqtype=sense_lib_type,
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


rule edit_10X_igblast:
    input:
        tsv="{base}/10X/igblast/{lib}_igblast.airr.tsv",
        airr_10X="{base}/10X/{lib}/outs/airr_rearrangement.tsv",
    output:
        tsv="{base}/10X/igblast/{lib}_igblast_edit.airr.tsv",
    params:
        organism="human",
    run:
        df = pd.read_table(input.tsv, sep="\t")
        df["library"] = wildcards.lib
        df.loc[:, "cell_id"] = wildcards.lib + df["sequence_id"]
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


rule combine_bracer_contigs:
    input:
        expand("{base}/SS2/{donor}/contigs.fasta", base=base, donor=donors),
    output:
        "{base}/SS2/combined_contigs.fasta",
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
        tsv="{base}/SS2/igblast/igblast_cell.airr.tsv",
    params:
        organism="human",
    run:
        df = pd.read_table(input.tsv, sep="\t")
        df["cell_id"] = df["sequence_id"].str.split("|", expand=True)[-1]
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

        print(input.di)
        fastas = glob.glob(input.di + "/*/tracer/assembled/*/filtered*/*.fa*")
        records = []
        print(fastas)
        for fasta in fastas:
            cellname = fasta.split("/")[-3]
            donor = fasta.split("/")[-6]
            for record in SeqIO.parse(fasta, "fasta"):
                record.description = "{}|{}|{}".format(
                    record.description, donor, cellname
                )
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
        -query {input.fasta} \
        -out {output}
        """


## Combined Outputs
rule combine_igblast:
    input:
        TenXs=expand(
            "{base}/10X/igblast/{lib}_igblast_edit.airr.tsv", lib=libs, base=base
        ),
        bracer="{base}/SS2/igblast/bracer.airr.tsv",
        tracer="{base}/SS2/igblast/tracer.airr.tsv",
    output:
        tsv="{base}/vdj/combined_igblast.airr.tsv.gz",
    log:
        "{base}/log/combineigblast.log",
    run:
        dfs = []
        infiles = input.TenXs
        infiles.append(input.bracer)
        for i in infiles:
            df = pd.read_table(i, sep="\t")
            df["sample_id"] = i.split("/")[-1].split("_")[0]
            dfs.append(df)

        combined = pd.concat(dfs)
        combined.to_csv(output.tsv, sep="\t", index=False, header=True)


rule annotate_constant_region:
    input:
        "{base}/vdj/combined_igblast.airr.tsv.gz",
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
        "python {params.scripts}/blast_constant_region.py "
        "{input} "
        "--min_j_sequence_length 15 "
        "--allow_missing_cdr3 True "
        "--allow_ns_in_sequence True "
        "-ighc_db {params.ighc_db} "
        "-outdir {wildcards.base}/vdjc "
        "> {log}"
