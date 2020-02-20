configfile: "config.json"

import glob

nsplit = 20

wildcard_constraints:
    id = "[0-9]+"


def get_parsed(wildcards):
    """ Returns the parsed BLAST files for the combine function """
    file = "data/{annotation}/genes/seq.{{id}}.fa".format(**wildcards)
    parsed = "results/{annotation}/blast/parsed_{{id}}.tsv".format(**wildcards)
    ids = glob_wildcards(file)
    return expand(parsed, id=ids.id)


rule all:
    """ Results files """
    input:
        "results/ensembl99/blast.tsv",
        "results/ensembl99/blast_stuff.tsv"


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
        genome = "data/references/genome.fa"
    params:
        link = lambda wildcards: config["download"]["genome"]
    shell:
        "wget --quiet -O {output.genome}.gz {params.link} && "
        "gunzip {output}.gz"


rule create_gene_bed12:
    """ Create a BED12 for all genes and for the longuest spliced transcript """
    input:
        gtf = rules.download_annotation.output.annotation
    output:
        bed12 = "data/{annotation}/BED12/genes.bed"
    conda:
        "envs/python.yaml"
    script:
        "scripts/create_gene_bed12.py"


rule extract_gene_sequences:
    """ From the BED12 file, extract all sequence in a FASTA file """
    input:
        bed12 = rules.create_gene_bed12.output.bed12,
        genome = rules.download_genome.output.genome
    output:
        fasta = "data/{annotation}/genes/seq.fa.temp"
    conda:
        "envs/bedtools.yaml"
    shell:
        "bedtools getfasta -split -name -s -fi {input.genome} -bed {input.bed12} -fo {output.fasta}"


rule truncate_fasta_gene_names:
    """ Remove the coordinates from the identifiers in the sequence FASTA """
    input:
        fasta = rules.extract_gene_sequences.output.fasta
    output:
        fasta = "data/{annotation}/genes/seq.fa"
    shell:
        "sed -r 's/\:.+//' {input.fasta} > {output.fasta}"


rule extract_seq_len_from_fasta:
    """
        From the sequence FASTA, extract the nucleotide lenght of all sequences
    """
    input:
        fasta = rules.truncate_fasta_gene_names.output.fasta
    output:
        seq_len = "data/{annotation}/seq_len.tsv"
    conda:
        "envs/bioawk.yaml"
    shell:
        "cat {input.fasta} | "
        "bioawk -t -c fastx '{{ print $name, length($seq) }}' > {output.seq_len}"


rule creating_blast_db:
    """ Indexing the sequence FASTA for the BLAST software """
    input:
        fasta = rules.truncate_fasta_gene_names.output.fasta
    output:
        db = "data/{annotation}/blast.db.nhr"
    params:
        db_name = "data/{annotation}/blast.db"
    conda:
        "envs/blast.yaml"
    shell:
        "makeblastdb -in {input.fasta} -dbtype nucl -parse_seqids -out {params.db_name}"


rule splitting_fasta:
    """ Splitting the sequence FASTA for parallelization """
    input:
        fasta = rules.truncate_fasta_gene_names.output.fasta
    output:
        tkn = "data/{annotation}/genes/.tkn"
    conda:
        "envs/pyfasta.yaml"
    shell:
        "pyfasta split -n {nsplit} {input.fasta} && touch {output.tkn}"


rule blasting_fasta:
    """  """
    input:
        db = rules.creating_blast_db.output.db,
        tkn = rules.splitting_fasta.output.tkn
    output:
        results = "results/{annotation}/blast/{id}.tsv"
    params:
        db_name = rules.creating_blast_db.params.db_name,
        query = "data/{annotation}/genes/seq.{id}.fa"
    threads:
        32
    conda:
        "envs/blast.yaml"
    shell:
        "blastn -db {params.db_name} -query {params.query} -outfmt 6 -num_threads {threads} -max_hsps 1 -evalue 1e-20 -out {output.results}"


rule parsing_blast_results:
    """ Parsing the BLAST results. """
    input:
        blast = rules.blasting_fasta.output.results,
        seq_len = rules.extract_seq_len_from_fasta.output.seq_len
    output:
        parsed = "results/{annotation}/blast/parsed_{id}.tsv"
    conda:
        "envs/python.yaml"
    script:
        "scripts/parsing_blast_results.py"


rule combine_parsed_results:
    """ """
    input:
        gtf = rules.download_annotation.output.annotation,
        files = get_parsed,
        tkn = "data/{annotation}/genes/.tkn"
    output:
        results = "results/{annotation}/blast.tsv"
    conda:
        "envs/python.yaml"
    script:
        "scripts/combine_parsed_results.py"


rule explore_results:
    input:
        data = rules.combine_parsed_results.output.results
    output:
        results = "results/{annotation}/blast_stuff.tsv"
    conda:
        "envs/python.yaml"
    script:
        "scripts/explore_results.py"


rule prepare_cytoscape_data:
    input:
        gtf = rules.download_annotation.output.annotation,
        data = rules.combine_parsed_results.output.results
    output:
        node = "results/{annotation}/cytoscape/node_score{score}.tsv",
        edge = "results/{annotation}/cytoscape/edge_score{score}.tsv"
    conda:
        "envs/python.yaml"
    script:
        "scripts/prepare_cytoscape_data.py"


rule scoring_edges_data:
    input:
        gtf = rules.download_annotation.output.annotation,
        data = rules.combine_parsed_results.output.results
    output:
        "results/{annotation}/scoring/score{score}.tsv"
    conda:
        "envs/python.yaml"
    script:
        "scripts/scoring_edges_data.py"
