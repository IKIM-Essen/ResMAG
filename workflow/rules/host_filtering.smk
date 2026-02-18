from pathlib import Path


if config["human-filtering"]["use-local"]:

    rule copy_local_human_ref:
        output:
            fasta=get_human_ref(),
        params:
            local=get_human_local_folder(),
            folder=lambda wildcards, output: Path(output.fasta).parent,
            file=lambda wildcards, output: Path(output.fasta).name,
        log:
            "logs/human_ref_local_copy.log",
        group:
            "refGenome_depended"
        conda:
            "../envs/unix.yaml"
        shell:
            "(mkdir -p {params.folder} && "
            "tar cpfz - -C {params.local} {params.file} | "
            "(cd {params.folder} ; tar xpfz -)) > {log} 2>&1"

else:

    rule download_human_ref:
        output:
            fasta=get_human_ref(),
        params:
            download=config["human-filtering"]["download-path"],
            folder=lambda wildcards, output: Path(output.fasta).parent,
        log:
            "logs/human_ref_download.log",
        group:
            "refGenome_depended"
        conda:
            "../envs/unix.yaml"
        shell:
            "(mkdir -p {params.folder} && "
            "cd {params.folder} && "
            "wget {params.download}) > {log} 2>&1"


rule map_to_human:
    input:
        fastqs=get_prefiltered_fastqs,
        ref=get_human_ref(),
    output:
        bam=temp("results/{project}/human_filtering/alignments/{sample}.bam"),
    threads: 64
    log:
        "logs/{project}/human_filtering/map_to_human_{sample}.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "(minimap2 -a -xsr -t {threads} {input.ref} {input.fastqs} | "
        "samtools view -bh | "
        "samtools sort --threads {threads} -o {output.bam}) > {log} 2>&1"


rule index_human_alignment:
    input:
        rules.map_to_human.output.bam,
    output:
        bai=temp("results/{project}/human_filtering/alignments/{sample}.bam.bai"),
    threads: 20
    log:
        "logs/{project}/human_filtering/index_human_alignment_{sample}.log",
    conda:
        "../envs/minimap2.yaml"
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
    threads: 64
    log:
        "results/{project}/output/report/prerequisites/qc/{sample}_filter_human.log",
    conda:
        "../envs/minimap2.yaml"
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
    threads: 20
    conda:
        "../envs/unix.yaml"
    shell:
        "gzip -k {input} > {log} 2>&1"


if config["host-filtering"]["do-host-filtering"]:
    # if there is a different host than human,
    # this is run before filtering human reads

    use rule map_to_human as map_to_host with:
        input:
            fastqs=get_trimmed_fastqs,
            ref=config["host-filtering"]["ref-genome"],
        output:
            bam=temp("results/{project}/host_filtering/alignments/{sample}.bam"),
        params:
            ref=config["host-filtering"]["ref-genome"],
        threads: 20
        log:
            "logs/{project}/host_filtering/map_to_host_{sample}.log",

    use rule index_human_alignment as index_host_alignment with:
        input:
            rules.map_to_host.output.bam,
        output:
            bai=temp("results/{project}/host_filtering/alignments/{sample}.bam.bai"),
        threads: 3
        log:
            "logs/{project}/host_filtering/index_host_alignment_{sample}.log",

    rule filter_host:
        input:
            bam=rules.map_to_host.output.bam,
            bai=rules.index_host_alignment.output.bai,
        output:
            filtered=temp(
                expand(
                    "results/{{project}}/host_filtering/non_host/{{sample}}_{read}.fastq.gz",
                    read=["R1", "R2"],
                )
            ),
        threads: 64
        log:
            "results/{project}/output/report/prerequisites/qc/{sample}_filter_host.log",
        conda:
            "../envs/minimap2.yaml"
        shell:
            "(samtools fastq -F 3584 -f 77 {input.bam} | "
            "gzip -c > {output.filtered[0]} && "
            "samtools fastq -F 3584 -f 141 {input.bam} | "
            "gzip -c > {output.filtered[1]}) > {log} 2>&1"
