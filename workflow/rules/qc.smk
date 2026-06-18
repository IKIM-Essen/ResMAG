## read QC
# fastp in paired-end mode for Illumina paired-end data
# version in this wrapper: fastp=1.0.1
rule fastp:
    input:
        sample=get_fastqs,
    output:
        trimmed=temp(
            [
                "results/{project}/trimmed/fastp/{sample}.1.fastq.gz",
                "results/{project}/trimmed/fastp/{sample}.2.fastq.gz",
            ]
        ),
        html=temp("results/{project}/trimmed/fastp/{sample}.html"),
        json="results/{project}/output/report/prerequisites/qc/{sample}.fastp.json",
    log:
        "logs/{project}/fastp/{sample}.log",
    threads: 16
    params:
        adapters=get_adapters,
        extra="--qualified_quality_phred {phred} --length_required {minlen}".format(
            phred=(config["quality-criteria"]["min-PHRED"]),
            minlen=(config["quality-criteria"]["min-length-reads"]),
        ),
    wrapper:
        "v7.1.0/bio/fastp"


rule fastqc:
    input:
        rules.fastp.output.trimmed,
    output:
        zip1=temp(
            "results/{project}/trimmed/fastp/{sample}.1_fastqc.zip",
        ),
        zip2=temp(
            "results/{project}/trimmed/fastp/{sample}.2_fastqc.zip",
        ),
    log:
        "logs/{project}/fastqc/{sample}.log",
    conda:
        "../envs/fastqc.yaml"
    threads: 2
    resources:
        mem_mb=1024,
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
        "results/{project}/output/report/multiqc_{project}.html",
        "results/{project}/output/report/multiqc_{project}.zip",
    log:
        "logs/{project}/multiqc.log",
    params:
        use_input_files_only=True,
    wrapper:
        "v8.1.1/bio/multiqc"


rule qc_summary:
    input:
        # move these to report_prerequistes
        json=rules.fastp.output.json,#"results/{project}/output/report/prerequisites/qc/{sample}.fastp.json",
        human_log="results/{project}/output/report/prerequisites/qc/{sample}_filter_human.log",
    output:
        csv="results/{project}/output/report/{sample}/{sample}_summary_quality.csv",
    log:
        "logs/{project}/report/{sample}_qc_summary.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/quality_summary.py"
