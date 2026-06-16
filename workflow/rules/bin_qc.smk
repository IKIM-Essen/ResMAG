## bin QC
if not config["checkm2"]["use-shared"]:

    rule checkm2_DB_download:
        output:
            dbfile=get_checkm2_db(),
        log:
            "logs/checkm2_DB_download.log",
        conda:
            "../envs/checkm2.yaml"
        params:
            direct=lambda wildcards, output: Path(output.dbfile).parent.parent,
        shell:
            "checkm2 database --download --path {params.direct} > {log} 2>&1"


rule checkm2_run:
    input:
        bins="results/{project}/output/fastas/{sample}/bins/",
        dbfile=get_checkm2_db(),
    output:
        stats="results/{project}/output/report/prerequisites/binning/{sample}/checkm2_quality_report.tsv",
    log:
        "logs/{project}/checkm2/{sample}.log",
    conda:
        "../envs/checkm2.yaml"
    threads: 24
    params:
        outdir="results/{project}/qc/checkm2/{sample}/",
    shell:
        "(checkm2 predict -x fa.gz --threads {threads} --force "
        "--input {input.bins}/ --output-directory {params.outdir}/ "
        "--database_path {input.dbfile} && "
        "cp {params.outdir}/quality_report.tsv {output.stats}) > {log} 2>&1"


rule bin_summary_sample:
    input:
        tool="results/{project}/output/report/prerequisites/binning/{sample}/{sample}_DASTool_summary.tsv",
        checkm="results/{project}/output/report/prerequisites/binning/{sample}/checkm2_quality_report.tsv",
        gtdb="results/{project}/output/classification/bins/{sample}/{sample}.summary.tsv",
    output:
        csv_bins="results/{project}/output/report/{sample}/{sample}_summary_bins.csv",
        csv_mags="results/{project}/output/report/{sample}/{sample}_summary_mags.csv",
    log:
        "logs/{project}/bin_summary/{sample}.log",
    conda:
        "../envs/python.yaml"
    threads: 4
    params:
        max_cont=config["MAG-criteria"]["max-contamination"],
        min_comp=config["MAG-criteria"]["min-completeness"],
    script:
        "../scripts/bin_summary.py"
