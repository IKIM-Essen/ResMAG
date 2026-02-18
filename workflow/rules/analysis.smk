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
    threads: 32
    conda:
        "../envs/pprodigal.yaml"
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
    threads: 20
    log:
        "logs/{project}/proteins/{sample}_gzip.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "gzip {input.faa} {input.fna} {input.gff} > {log} 2>&1"


# Plasmid analysis
rule load_genomad_DB:
    output:
        folder=get_genomad_DB_folder(),
        file=get_genomad_DB_file(),
    params:
        res_folder=lambda wildcards, output: Path(output.folder).parent,
    log:
        "logs/load_genomad_DB.log",
    conda:
        "../envs/genomad.yaml"
    shell:
        "genomad download-database {params.res_folder}/ > {log} 2>&1"


rule genomad_run:
    input:
        db=rules.load_genomad_DB.output.folder,
        asmbl=rules.gzip_assembly.output,
    output:
        #outdir=temp(directory("results/{project}/genomad/{sample}/")),
        plasmid_tsv=temp(
            "results/{project}/genomad/{sample}/{sample}_summary/{sample}_plasmid_summary.tsv"
        ),
        virus_tsv=temp(
            "results/{project}/genomad/{sample}/{sample}_summary/{sample}_virus_summary.tsv"
        ),
    params:
        outdir=lambda wildcards, output: Path(output.plasmid_tsv).parent.parent,
    log:
        "logs/{project}/plasmids/{sample}_run.log",
    threads: 64
    conda:
        "../envs/genomad.yaml"
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
    params:
        outdir=lambda wildcards, output: Path(output.plasmid_tsv).parent,
        sum_folder=lambda wildcards, input: Path(input.plasmid_tsv).parent,
    log:
        "logs/{project}/plasmids/{sample}_move_output.log",
    threads: 32
    conda:
        "../envs/unix.yaml"
    shell:
        "scp {params.sum_folder}/*.tsv {params.sum_folder}/*.json {params.outdir}/ > {log} 2>&1"


# Resistance analysis
rule download_CARD_data:
    output:
        json=get_card_db_file(),
    params:
        download=config["card"]["url"],
        folder=lambda wildcards, output: Path(output.json).parent,
        filename=lambda wildcards, output: Path(output.json).name,
    log:
        "logs/CARD_data_download.log",
    threads: 30
    conda:
        "../envs/unix.yaml"
    shell:
        "(cd {params.folder} && "
        "wget {params.download} && "
        "tar -xvf data ./{params.filename}) > {log} 2>&1"


rule CARD_load_DB_for_reads:
    input:
        get_card_db_file(),
    output:
        temp(touch("results/CARD_load_DB_for_reads.done")),
    log:
        "logs/CARD_load_DB_for_reads.log",
    conda:
        "../envs/card.yaml"
    shell:
        "rgi clean --local && "
        "rgi load --card_json {input} --local > {log} 2>&1"


rule CARD_annotation:
    input:
        json=get_card_db_file(),
        load=rules.CARD_load_DB_for_reads.output,
    output:
        done=temp(touch("results/CARD_annotation.done")),
        ann=get_card_annotation_file(),
    params:
        folder=lambda wildcards, input: Path(input.json).parent,
        file=lambda wildcards, input: Path(input.json).name,
    log:
        "logs/CARD_annotation.log",
    conda:
        "../envs/card.yaml"
    shell:
        "(cd {params.folder}/ && "
        "rgi card_annotation -i {params.file}) && "
        "rgi load -i {input.json} --card_annotation {output.ann} --local > {log} 2>&1"


rule CARD_read_run:
    input:
        fa=get_filtered_gz_fastqs,
        db=rules.CARD_annotation.output,
    output:
        txt="results/{project}/output/resistance/CARD/reads/{sample}/{sample}.gene_mapping_data.txt",
    params:
        folder=lambda wildcards, output: Path(output.txt).parent,
    log:
        "logs/{project}/ARGs/reads/{sample}.log",
    threads: 64
    conda:
        "../envs/card.yaml"
    shell:
        "rgi bwt -1 {input.fa[0]} -2 {input.fa[1]} "
        "-o {params.folder}/{wildcards.sample} --local "
        "--clean -n {threads} > {log} 2>&1"


rule CARD_read_sample_summary:
    input:
        txt=rules.CARD_read_run.output.txt,
    output:
        csv="results/{project}/output/resistance/CARD/reads/{sample}/{sample}_read_ARGs.csv",
    params:
        case="reads",
    log:
        "logs/{project}/ARGs/reads/{sample}.log",
    threads: 2
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/arg_summary_sample.py"


# updates CARD database to use for contigs instead of reads
# read based classification must be finished before
rule CARD_load_DB:
    input:
        db=get_card_db_file(),
        read_args=expand(
            "results/{project}/output/resistance/CARD/reads/{sample}/{sample}.gene_mapping_data.txt",
            sample=get_samples(),
            project=get_project(),
        ),
    output:
        touch("results/CARD_load_DB.done"),
    log:
        "logs/CARD_load_DB.log",
    conda:
        "../envs/card.yaml"
    shell:
        "rgi clean --local && "
        "rgi load --card_json {input.db} --local > {log} 2>&1"


rule CARD_assembly_run:
    input:
        faa=rules.gzip_proteins.output.faa,
        db=rules.CARD_load_DB.output,
    output:
        txt="results/{project}/output/resistance/CARD/assembly/{sample}/{sample}.txt",
        json="results/{project}/output/resistance/CARD/assembly/{sample}/{sample}.json",
    params:
        path_wo_ext=lambda wildcards, output: Path(output.txt).with_suffix(""),
    log:
        "logs/{project}/ARGs/assembly/{sample}.log",
    threads: 64
    conda:
        "../envs/card.yaml"
    shell:
        "rgi main -i {input.faa} -o {params.path_wo_ext} "
        "-t protein -a DIAMOND --local "
        "-n {threads} --clean > {log} 2>&1"


use rule CARD_read_sample_summary as CARD_assembly_sample_summary with:
    input:
        txt=rules.CARD_assembly_run.output.txt,
        json=rules.CARD_assembly_run.output.json,
    output:
        csv="results/{project}/output/resistance/CARD/assembly/{sample}/{sample}_assembly_ARGs.csv",
    params:
        case="assembly",
    log:
        "logs/{project}/ARGs/assembly/{sample}.log",


use rule CARD_assembly_run as CARD_mag_run with:
    input:
        fa="results/{project}/output/fastas/{sample}/mags/{binID}.fa.gz",
        db=rules.CARD_load_DB.output,
    output:
        txt="results/{project}/output/resistance/CARD/mags/{sample}/{binID}.txt",
    params:
        path_wo_ext=lambda wildcards, output: Path(output.txt).with_suffix(""),
    log:
        "logs/{project}/ARGs/mags/{sample}/{binID}.log",


rule wrap_mag_ARGs:
    input:
        get_mag_ARGs,
    output:
        touch("results/{project}/output/ARGs/mags/{sample}/all_mags.done"),
    log:
        "logs/{project}/ARGs/mags/{sample}/all_mags.log",
