# if there is no local database version to use, it is downloaded
if not config["kaiju"]["use-local"]:

    rule download_kaiju:
        output:
            db_files=get_kaiju_files(),
        params:
            download=config["kaiju"]["download"],
            db_folder=lambda wildcards, output: Path(output.db_files[0]).parent,
        log:
            "logs/kaiju_DB_download.log",
        conda:
            "../envs/unix.yaml"
        shell:
            "(mkdir -p {params.db_folder} && "
            "wget -c {params.download} -O - | "
            "tar -zxv -C {params.db_folder}) > {log} 2>&1"


rule run_kaiju:
    input:
        db_files=get_kaiju_files(),
        fastqs=get_filtered_gz_fastqs,
    output:
        report=temp("results/{project}/output/classification/reads/{sample}/kaiju.out"),
    params:
        evalue=0.00001,
    threads: 64
    log:
        "logs/{project}/kaiju/{sample}_run.log",
    conda:
        "../envs/kaiju.yaml"
    shell:
        "kaiju -z {threads} -t {input.db_files[0]} -f {input.db_files[1]} "
        "-v -E {params.evalue} -i {input.fastqs[0]} "
        "-j {input.fastqs[1]} -o {output.report} > {log} 2>&1"


rule kaiju_table:
    input:
        report=rules.run_kaiju.output.report,
        db_files=get_kaiju_files(),
    output:
        summary="results/{project}/output/classification/reads/{sample}/{sample}_{level}_kaiju_summary.tsv",
    threads: 2
    log:
        "logs/{project}/kaiju/{sample}_{level}_table.log",
    conda:
        "../envs/kaiju.yaml"
    shell:
        "(kaiju2table -t {input.db_files[0]} -n {input.db_files[2]} "
        "-r {wildcards.level} -o {output.summary} "
        "{input.report}) > {log} 2>&1"
        # -p to print full taxonomy


rule kaiju2krona:
    input:
        db_files=get_kaiju_files(),
        report=rules.run_kaiju.output.report,
    output:
        krona=temp(
            "results/{project}/output/classification/reads/{sample}/kaiju.out.krona"
        ),
        html=report(
            "results/{project}/output/report/{sample}/{sample}_kaiju.out.html",
            htmlindex="index.html",
            category="5. Taxonomic classification",
            subcategory="5.2 Read classification",
            labels={"sample": "{sample}"},
        ),
    log:
        "logs/{project}/kaiju/{sample}_2krona.log",
    conda:
        "../envs/kaiju.yaml"
    shell:
        "(kaiju2krona -t {input.db_files[0]} -n {input.db_files[2]} "
        "-i {input.report} -o {output.krona} && "
        "ktImportText -o {output.html} {output.krona}) > {log} 2>&1"


rule prepare_gtdb:
    output:
        done=temp(touch("results/GTDB_prep.done")),
    params:
        db_folder=get_gtdb_folder(),
    threads: 1
    log:
        "logs/GTDB_prep.log",
    conda:
        "../envs/gtdbtk.yaml"
    shell:
        "(conda env config vars set "
        "GTDBTK_DATA_PATH='{params.db_folder}') > {log} 2>&1"


rule gtdbtk_classify_wf:
    input:
        bins=rules.gzip_bins.output.bins,
        db_prep=rules.prepare_gtdb.output.done,
    output:
        json="results/{project}/output/classification/bins/{sample}/gtdbtk.json",
        outdir=temp(
            directory("results/{project}/output/classification/bins/{sample}/gtdbtk/")
        ),
    params:
        clf_outdir=lambda wildcards, output: Path(output.json).parent,
        json=lambda wildcards, output: Path(output.json).name,
    threads: 40
    # to ensure only one of these commands are run at the same time
    resources:
        heavy=1,
    log:
        "logs/{project}/gtdbtk/{sample}_classify.log",
    conda:
        "../envs/gtdbtk.yaml"
    shell:
        "(gtdbtk classify_wf --prefix {wildcards.sample} -x fa.gz "
        "--cpus {threads} --pplacer_cpus {threads} "
        "--genome_dir {input.bins}/ --out_dir {output.outdir}/ && "
        "cp {output.outdir}/{params.json} {params.clf_outdir}/) > {log} 2>&1"


# combines bacterial and archaeal
rule gtdb_summary:
    input:
        json=rules.gtdbtk_classify_wf.output.json,
        gtdb_dir=rules.gtdbtk_classify_wf.output.outdir,
    output:
        summary="results/{project}/output/classification/bins/{sample}/{sample}.summary.tsv",
    params:
        outdir=lambda wildcards, input: Path(input.json).parent,
    threads: 10
    log:
        "logs/{project}/gtdbtk/{sample}_summary.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/gtdb_summary.py"
