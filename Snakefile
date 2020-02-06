configfile: "config.json"

rule all:
    input:
        "stuff"


rule download_annotation:
    """ Download the different annotations as .GTF.GZ """
    output:
        annotation = "data/references/{annotation}.gtf.gz"
    params:
        link = lambda wildcards: config["download"][wildcards.annotation]
    shell:
        "wget --quiet -O {output.annotation} {params.link}"


rule download_genome:
    """ Download the genome as .FA """
    output:
        genome = "data/references/{genome}.fa"
    params:
        link = lambda wildcards: config["download"][wildcards.genome]
    shell:
        "wget --quiet -O {output.genome}.gz {params.link} && "
        "gunzip {output}.gz"


rule create_gene_bed12:
    input:
        annotation = download_annotation.output.annotation
    output:
        bed12 = "data/BED12/genes_{annotation}.bed"
    shell:
        "zcat {input.annotation} | "
        "awk '{{if ($3 ==\"gene\") {{print $0}}}}'"
#grep -oP 'gene_id "\K.*?(?=")'
