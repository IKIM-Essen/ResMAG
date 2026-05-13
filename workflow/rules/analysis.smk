# gene identification
rule gene_identification:
    input:
        contigs=get_assembly,
    output:
        faa="results/{project}/output/proteins/{sample}/{sample}_proteins.faa",
        fna="results/{project}/output/proteins/{sample}/{sample}_nucleotides.fna",
        gff="results/{project}/output/proteins/{sample}/{sample}_annotations.gff",
    log:
        "logs/{project}/proteins/{sample}.log",
    conda:
        "../envs/pprodigal.yaml"
    threads: 32
    shell:
        "pprodigal -i {input.contigs} -o {output.gff} -a {output.faa} "
        "-d {output.fna} -p meta --tasks {threads} > {log} 2>&1"


rule gzip_proteins:
    input:
        faa=rules.gene_identification.output.faa,
        fna=rules.gene_identification.output.fna,
        gff=rules.gene_identification.output.gff,
    output:
        faa="results/{project}/output/proteins/{sample}/{sample}_proteins.faa.gz",
        fna="results/{project}/output/proteins/{sample}/{sample}_nucleotides.fna.gz",
        gff="results/{project}/output/proteins/{sample}/{sample}_annotations.gff.gz",
    log:
        "logs/{project}/proteins/{sample}_gzip.log",
    priority: 1
    conda:
        "../envs/unix.yaml"
    threads: 20
    shell:
        "gzip {input.faa} {input.fna} {input.gff} > {log} 2>&1"


# Plasmid analysis

if not config["genomad"]["use-shared"]:

    rule load_genomad_DB:
        output:
            folder=get_genomad_DB_folder(),
            file=get_genomad_DB_file(),
        log:
            "logs/load_genomad_DB.log",
        conda:
            "../envs/genomad.yaml"
        params:
            res_folder=lambda wildcards, output: Path(output.folder).parent,
        shell:
            "genomad download-database {params.res_folder}/ > {log} 2>&1"


rule genomad_run:
    input:
        db=get_genomad_DB_folder(),
        asmbl=rules.gzip_assembly.output,
    output:
        plasmid_tsv=temp(
            "results/{project}/genomad/{sample}/{sample}_summary/{sample}_plasmid_summary.tsv"
        ),
        virus_tsv=temp(
            "results/{project}/genomad/{sample}/{sample}_summary/{sample}_virus_summary.tsv"
        ),
    log:
        "logs/{project}/plasmids/{sample}_run.log",
    conda:
        "../envs/genomad.yaml"
    threads: 64
    params:
        outdir=lambda wildcards, output: Path(output.plasmid_tsv).parent.parent,
    shell:
        "genomad end-to-end --cleanup -t {threads} "
        "{input.asmbl} {params.outdir}/ "
        "{input.db}/ > {log} 2>&1"


rule move_genomad_output:
    input:
        plasmid_tsv=rules.genomad_run.output.plasmid_tsv,
        virus_tsv=rules.genomad_run.output.virus_tsv,
    output:
        plasmid_tsv="results/{project}/output/plasmids/{sample}/{sample}_plasmid_summary.tsv",
        virus_tsv="results/{project}/output/plasmids/{sample}/{sample}_virus_summary.tsv",
    log:
        "logs/{project}/plasmids/{sample}_move_output.log",
    conda:
        "../envs/unix.yaml"
    threads: 32
    params:
        outdir=lambda wildcards, output: Path(output.plasmid_tsv).parent,
        sum_folder=lambda wildcards, input: Path(input.plasmid_tsv).parent,
    shell:
        "cp {params.sum_folder}/*.tsv {params.sum_folder}/*.json {params.outdir}/ > {log} 2>&1"
