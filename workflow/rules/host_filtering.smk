from pathlib import Path

# Download human reference genome if not local file is given
if not config["human-filtering"]["use-shared"]:

    rule download_human_ref:
        output:
            fasta=get_human_ref(),
        log:
            "logs/human_ref_download.log",
        group:
            "refGenome_depended"
        conda:
            "../envs/unix.yaml"
        params:
            download=get_human_ref_download(),
            folder=lambda wildcards, output: Path(output.fasta).parent,
        shell:
            "(mkdir -p {params.folder} && "
            "cd {params.folder} && "
            "wget {params.download}) > {log} 2>&1"


rule map_to_human:
    input:
        fastqs=get_trimmed_fastqs,
        ref=get_human_ref(),
    output:
        bam=temp("results/{project}/human_filtering/alignments/{sample}.bam"),
    log:
        "logs/{project}/human_filtering/map_to_human_{sample}.log",
    conda:
        "../envs/minimap2.yaml"
    threads: 64
    shell:
        "(minimap2 -a -xsr -t {threads} {input.ref} {input.fastqs} | "
        "samtools view -bh | "
        "samtools sort --threads {threads} -o {output.bam}) > {log} 2>&1"


rule index_human_alignment:
    input:
        rules.map_to_human.output.bam,
    output:
        bai=temp("results/{project}/human_filtering/alignments/{sample}.bam.bai"),
    log:
        "logs/{project}/human_filtering/index_human_alignment_{sample}.log",
    conda:
        "../envs/minimap2.yaml"
    threads: 20
    shell:
        "samtools index {input} > {log} 2>&1"


rule filter_human:
    input:
        bam=rules.map_to_human.output.bam,
        bai=rules.index_human_alignment.output.bai,
    output:
        filtered=temp(
            expand(
                "results/{{project}}/output/filtered_reads/{{sample}}_{read}.fastq",
                read=["R1", "R2"],
            )
        ),
    log:
        "results/{project}/output/report/prerequisites/qc/{sample}_filter_human.log",
    conda:
        "../envs/minimap2.yaml"
    threads: 64
    shell:
        "(samtools fastq --threads {threads} -F 3584 -f 77 "
        "-o {output.filtered[0]} {input.bam} && "
        "samtools fastq --threads {threads} -F 3584 -f 141 "
        "-o {output.filtered[1]} {input.bam}) > {log} 2>&1"


rule gzip_filtered_reads:
    input:
        "results/{project}/output/filtered_reads/{sample}_{read}.fastq",
    output:
        "results/{project}/output/filtered_reads/{sample}_{read}.fastq.gz",
    log:
        "logs/{project}/human_filtering/gzip_{sample}_{read}.log",
    priority: 1
    conda:
        "../envs/unix.yaml"
    threads: 20
    shell:
        "gzip -k {input} > {log} 2>&1"
