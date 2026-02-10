## read QC
# fastp in paired-end mode for Illumina paired-end data
# version in this wrapper: fastp=1.0.1
rule fastp:
    input:
        sample=get_local_fastqs,
    output:
        trimmed=temp(
            [
                "results/{project}/trimmed/fastp/{sample}.1.fastq.gz",
                "results/{project}/trimmed/fastp/{sample}.2.fastq.gz",
            ]
        ),
        html=temp("results/{project}/trimmed/fastp/{sample}.html"),
        json="results/{project}/output/report/prerequisites/qc/{sample}.fastp.json",
    params:
        adapters=get_adapters,
        extra="--qualified_quality_phred {phred} --length_required {minlen}".format(
            phred=(config["quality-criteria"]["min-PHRED"]),
            minlen=(config["quality-criteria"]["min-length-reads"]),
        ),
    log:
        "logs/{project}/fastp/{sample}.log",
    threads: 16
    wrapper:
        "v7.1.0/bio/fastp"


"""# version in this wrapper: fastqc=0.12.1
rule fastqc:
    input:
        rules.fastp.output.trimmed,
        #get_trimmed_fastqs,
    output:
        html=temp("results/{project}/qc/fastqc/{sample}_trimmed.html"),
        zip=temp("results/{project}/qc/fastqc/{sample}_trimmed_fastqc.zip"),
    threads: 4
    resources:
        mem_mb=1024,
    log:
        "logs/{project}/fastqc/{sample}.log",
    wrapper:
        "v7.6.0/bio/fastqc"
"""


rule fastqc:
    input:
        rules.fastp.output.trimmed,
    output:
        zip=temp(
            expand(
                "results/{{project}}/trimmed/fastp/{{sample}}.{read}_fastqc.zip",
                read=["1", "2"],
            )
        ),
        html=temp(
            expand(
                "results/{{project}}/trimmed/fastp/{{sample}}.{read}_fastqc.html",
                read=["1", "2"],
            )
        ),
    threads: 4
    resources:
        mem_mb=1024,
    log:
        "logs/{project}/fastqc/{sample}.log",
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc --memory {resources.mem_mb} --threads {threads} "
        "--format fastq --quiet {input} > {log} 2>&1"


# version in this wrapper: multiqc=1.33
rule multiqc:
    input:
        expand(
            [
                "results/{{project}}/trimmed/fastp/{sample}.{read}_fastqc.zip",
                "results/{{project}}/output/report/prerequisites/qc/{sample}.fastp.json",
            ],
            sample=get_samples(),
            read=["1", "2"],
        ),
        config="config/multiqc_config.yaml",
    output:
        report(
            "results/{project}/output/report/all/multiqc_{project}.html",
            htmlindex="multiqc_{project}.html",
            category="1. Quality control",
            labels={"sample": "all samples"},
        ),
        "results/{project}/output/report/all/multiqc_{project}.zip",
    params:
        extra=("--title 'Results for data from {project} project'"),
        use_input_files_only=True,
    log:
        "logs/{project}/multiqc.log",
    wrapper:
        "v8.1.1/bio/multiqc"


rule qc_summary:
    input:
        # move these to report_prerequistes
        jsons=expand(
            "results/{{project}}/output/report/prerequisites/qc/{sample}.fastp.json",
            sample=get_samples(),
        ),
        human_logs=expand(
            "results/{{project}}/output/report/prerequisites/qc/{sample}_filter_human.log",
            sample=get_samples(),
        ),
        host_logs=get_host_map_statistics,
    output:
        csv="results/{project}/output/report/all/quality_summary.csv",
        vis_csv=temp("results/{project}/output/report/all/quality_summary_visual.csv"),
    params:
        other_host=config["host-filtering"]["do-host-filtering"],
        hostname=config["host-filtering"]["host-name"],
    log:
        "logs/{project}/report/qc_summary.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/quality_summary.py"


rule qc_summary_report:
    input:
        rules.qc_summary.output.vis_csv,
    output:
        temp(
            report(
                directory("results/{project}/output/report/all/quality_summary/"),
                htmlindex="index.html",
                category="1. Quality control",
                labels={
                    "sample": "all samples",
                },
            )
        ),
    params:
        pin_until="sample",
        styles="resources/report/tables/",
        name="quality_summary",
        header="Quality summary based on fastp report and mapping to host genome(s)",
        pattern=config["tablular-config"],
    log:
        "logs/{project}/report/qc_summary_rbt_csv.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt csv-report {input} --pin-until {params.pin_until} {output} && "
        "(sed -i '{params.pattern} {params.header}</a>' "
        "{output}/indexes/index1.html && "
        "sed -i 's/report.xlsx/{params.name}_report.xlsx/g' {output}/indexes/index1.html) && "
        "mv {output}/report.xlsx {output}/{params.name}_report.xlsx && "
        "cp {params.styles}* {output}/css/ > {log} 2>&1"
