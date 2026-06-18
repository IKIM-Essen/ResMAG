from pathlib import Path


rule megahit:
    input:
        fastqs=get_filtered_gz_fastqs,
    output:
        contigs=temp("results/{project}/megahit/{sample}/final.contigs.fa"),
        outdir=temp(directory("results/{project}/megahit/{sample}/")),
        log="results/{project}/output/report/prerequisites/assembly/{sample}_megahit.log",
    params:
        threshold=get_contig_length_threshold(),
    threads: 64
    priority: 2
    resources:
        heavy=2,
    log:
        "logs/{project}/assembly/{sample}_megahit.log",
    conda:
        "../envs/megahit.yaml"
    shell:
        "(megahit -1 {input.fastqs[0]} -2 {input.fastqs[1]} "
        "--min-contig-len {params.threshold} -t {threads} "
        "--out-dir {output.outdir} -f > {log} 2>&1) && "
        "cp {log} {output.log}"


rule map_to_assembly:
    input:
        contigs=get_assembly,
        fastqs=get_filtered_gz_fastqs,
    output:
        bam=temp(
            "results/{project}/output/report/prerequisites/assembly/{sample}_reads_mapped.bam"
        ),
    threads: 60
    log:
        "logs/{project}/assembly/{sample}_mapping_reads.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "(minimap2 -a -xsr -t {threads} {input.contigs} {input.fastqs} | "
        "samtools view -bh | "
        "samtools sort --threads {threads} -o {output.bam}) > {log} 2>&1"


rule index_assembly_alignment:
    input:
        rules.map_to_assembly.output.bam,
    output:
        bai=temp(
            "results/{project}/output/report/prerequisites/assembly/{sample}_reads_mapped.bam.bai"
        ),
    threads: 2
    log:
        "logs/{project}/assembly/{sample}_mapping_reads_index.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "samtools index {input} > {log} 2>&1"


rule reads_mapped_assembly:
    input:
        bam=rules.map_to_assembly.output.bam,
        bai=rules.index_assembly_alignment.output.bai,
    output:
        txt="results/{project}/output/report/prerequisites/assembly/{sample}_reads_mapped.txt",
    threads: 16
    log:
        "logs/{project}/assembly/{sample}_mapping_reads.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "samtools view -c -F 4 --threads {threads} "
        "-o {output.txt} {input.bam} > {log} 2>&1"


rule gzip_assembly:
    input:
        contigs=get_assembly,
    output:
        "results/{project}/output/fastas/{sample}/{sample}.fa.gz",
    threads: 16
    priority: 1
    log:
        "logs/{project}/assembly/{sample}_gzip.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "gzip -c {input.contigs} > {output} 2> {log}"


rule cleanup_assembly:
    input:
        # dependency only
        txt=rules.reads_mapped_assembly.output.txt,

        # actually delete these
        to_delete=[
            rules.map_to_assembly.output.bam,
            rules.index_assembly_alignment.output.bai,
        ] 
    output:
        touch("results/{project}/cleanup/{sample}_assembly.done"),
    log:
        "logs/{project}/cleanup_assembly/{sample}.log",
    conda:
        "../envs/unix.yaml"
    threads: 1
    priority: 1
    script:
        "../scripts/cleanup_files.py"


rule assembly_summary:
    input:
        qc_csv=rules.qc_summary.output.csv,
        asbl="results/{project}/output/report/prerequisites/assembly/{sample}_megahit.log",
        mapped="results/{project}/output/report/prerequisites/assembly/{sample}_reads_mapped.txt",
        csv_bins="results/{project}/output/report/{sample}/{sample}_summary_bins.csv",
        csv_mags="results/{project}/output/report/{sample}/{sample}_summary_mags.csv",
    output:
        csv="results/{project}/output/report/{sample}/{sample}_summary_assembly.csv",
    log:
        "logs/{project}/report/{sample}_summary_assembly.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/assembly_summary.py"
