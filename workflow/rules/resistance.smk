# Resistance analysis
if not config["card"]["use-shared"]:

    rule download_CARD_data:
        output:
            json=get_card_db_file(),
        log:
            "logs/CARD_data_download.log",
        conda:
            "../envs/unix.yaml"
        threads: 30
        params:
            download=config["card"]["url"],
            folder=lambda wildcards, output: Path(output.json).parent,
            filename=lambda wildcards, output: Path(output.json).name,
        shell:
            "(cd {params.folder} && "
            "wget {params.download} && "
            "tar -xvf data ./{params.filename}) > {log} 2>&1"


rule CARD_load_DB_for_reads:
    input:
        get_card_db_file(),
    output:
        touch("results/CARD_load_DB_for_reads.done"),
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
        ann=get_card_annotation_file(),
    log:
        "logs/CARD_annotation.log",
    conda:
        "../envs/card.yaml"
    params:
        folder=lambda wildcards, input: Path(input.json).parent,
        file=lambda wildcards, input: Path(input.json).name,
    shell:
        "(cd {params.folder}/ && "
        "rgi card_annotation -i {params.file}) && "
        "rgi load -i {input.json} --card_annotation {output.ann} --local > {log} 2>&1"


rule CARD_read_run:
    input:
        fa=get_filtered_gz_fastqs,
        db=rules.CARD_annotation.output,
        load=rules.CARD_load_DB_for_reads.output,
    output:
        txt="results/{project}/output/resistance/CARD/reads/{sample}/{sample}.gene_mapping_data.txt",
    log:
        "logs/{project}/ARGs/reads/{sample}.log",
    conda:
        "../envs/card.yaml"
    threads: 60
    params:
        folder=lambda wildcards, output: Path(output.txt).parent,
    shell:
        "rgi bwt -1 {input.fa[0]} -2 {input.fa[1]} "
        "-o {params.folder}/{wildcards.sample} --local "
        "--clean -n {threads} > {log} 2>&1"


rule uniCARD_makeDB:
    input:
        db=get_uniCARD_db(),
    output:
        dmdb=get_unicard_dmnd(),
    log:
        "logs/unicard_makeDB.log",
    conda:
        "../envs/diamond.yaml"
    threads: 60
    params:
        filename=lambda wildcards, input: Path(input.db).name,
        filename_wo_ext=Path(get_uniCARD_db_wo_ext()).name,
        folder=lambda wildcards, input: Path(input.db).parent,
    shell:
        "(cd {params.folder} && diamond makedb --in {params.filename} -d {params.filename_wo_ext}) > {log} 2>&1"


rule uniCARD_run:
    input:
        faa=get_proteins,
        dmdb=get_unicard_dmnd(),
    output:
        tsv=temp(
            "results/{project}/output/resistance/uniCARD/assembly/{sample}_out.tsv"
        ),
    log:
        "logs/{project}/uniCARD/{sample}.log",
    # ensure priority over potential contig classification
    priority: 1
    conda:
        "../envs/diamond.yaml"
    threads: 60
    # to ensure only one of these commands are run at the same time
    resources:
        heavy=1,
    params:
        db_wo_ext=lambda w, input: os.path.splitext(input.dmdb)[0],
        max_targets=15,
    shell:
        "diamond blastp -q {input.faa} -d {params.db_wo_ext} "
        "-p {threads} --max-target-seqs {params.max_targets} "
        "-o {output.tsv} > {log} 2>&1"


rule gz_uniCARD_tsv:
    input:
        unicard_in=rules.uniCARD_run.output.tsv,
    output:
        tsv=temp(
            "results/{project}/output/resistance/uniCARD/assembly/{sample}_out.tsv.gz"
        ),
    log:
        "logs/{project}/uniCARD/{sample}_gz.log",
    conda:
        "../envs/unix.yaml"
    threads: 16
    shell:
        "(gzip -k {input.unicard_in}) > {log} 2>&1"


rule filter_uniCARD:
    input:
        tsv=rules.gz_uniCARD_tsv.output.tsv,
        json=get_CARD_hierarchy(),
    output:
        csv="results/{project}/output/resistance/uniCARD/assembly/{sample}.csv",
    log:
        "logs/{project}/uniCARD/{sample}_filter.log",
    conda:
        "../envs/python.yaml"
    threads: 16
    params:
        # filter window for filtering UniCARD results, default: 1%
        filter_window=0.01,
    shell:
        "python workflow/scripts/uniCARD_filtering.py --infile {input.tsv} "
        "--card_hierarchy {input.json} --outfile {output.csv} "
        "--filter_window {params.filter_window} > {log} 2>&1"


rule resistance_abundance_per_sample:
    input:
        csv=rules.filter_uniCARD.output.csv,
        json=get_CARD_hierarchy(),
        contig_file="results/{project}/output/classification/assembly/{sample}/{sample}_classified_contigs.csv",
        contig_len=rules.extract_contig_headers.output.contig_len,
    output:
        abd_class="results/{project}/output/resistance/uniCARD/per_sample_resistance_abundance/{sample}_drug_class_abundance.csv",
        abd_ARGs="results/{project}/output/resistance/uniCARD/per_sample_resistance_abundance/{sample}_ARGs_abundance.csv",
    log:
        "logs/{project}/uniCARD/{sample}_abundance.log",
    conda:
        "../envs/python.yaml"
    threads: 5
    script:
        "../scripts/uniCARD_abundance.py"


"""
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
"""
