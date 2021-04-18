# snakemake rules for pulling the immcantation image depending on
# whether docker or singularity is being used
# put this in config
immcantation_docker = "docker://kleinstein/immcantation:4.1.0"
baldr_address = "docker://bosingerlab/baldr"

rule get_BALDR_image:
        """ Pull the baldr image
        """
        output:
            "{}/resources/repertoire/{}.sif".format(workflow.basedir, "baldr"),
        threads: 1
        params:
            name="container_pull"
        resources:
            mem_mb=10000,
        shell:
            "module load system && "
            "mkdir -p $(dirname {output}) && cd $(dirname {output}) && "
            "img=$(basename {output}) && "
            "singularity pull --name $img "
            "{baldr_address}"
 
rule get_immcantation_image:
        """ pull the immcantation image using singularity

        """
        output:"{}/resources/repertoire/{}.sif".format(workflow.basedir, "immcantation"),
        threads: 1
        params:
            name="container_pull",
        resources:
            mem_mb=10000,
        shell:
            "module load system && "
            "mkdir -p $(dirname {output}) && cd $(dirname {output}) && "
            "img=$(basename {output}) && "
            "singularity pull --name $img "
            "{immcantation_address}"
            
