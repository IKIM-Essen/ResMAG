from pathlib import Path


rule postprocess_metabat:
    input:
        outdir=rules.metabat.output.outdir,
    output:
        tsv="results/{project}/output/contig2bins/{sample}/metabat_contig2bin.tsv",
    params:
        binner="metabat",
        prefix="bin",
    threads: 12
    log:
        "logs/{project}/contig2bins/{sample}/postprocess_metabat.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/metabat_postprocess.py"


rule postprocess_metacoag:
    input:
        c2bin="results/{project}/binning/metacoag/{sample}/contig_to_bin.tsv",
    output:
        tsv="results/{project}/output/contig2bins/{sample}/metacoag_contig2bin.tsv",
    threads: 12
    log:
        "logs/{project}/contig2bins/{sample}/postprocess_metacoag.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(awk '{{print $1 \"\t\" $NF}}' {input.c2bin} | "
        "sed 's/len=[0-9]*,/metacoag_/g' > {output.tsv}) > {log} 2>&1"


rule dastool_run:
    input:
        metacoag=rules.postprocess_metacoag.output.tsv,
        metabat=rules.postprocess_metabat.output.tsv,
        contigs=get_assembly,
    output:
        summary=temp(
            "results/{project}/binning/das_tool/{sample}/{sample}_DASTool_summary.tsv"
        ),
        contig2bin=temp(
            "results/{project}/binning/das_tool/{sample}/{sample}_DASTool_contig2bin.tsv"
        ),
        bins=temp(
            directory(
                "results/{project}/binning/das_tool/{sample}/{sample}_DASTool_bins/"
            )
        ),
    params:
        outdir=lambda wildcards, output: Path(output.bins).parent,
        threshold=0.001,
    threads: 64
    log:
        "logs/{project}/das_tool/{sample}/das_tool_run.log",
    conda:
        "../envs/das_tool.yaml"
    shell:
        "DAS_Tool --debug --write_bins "
        "-i {input.metabat},{input.metacoag} "
        "-l metabat,metacoag "
        "-c {input.contigs} "
        "-o {params.outdir}/{wildcards.sample} "
        "--score_threshold {params.threshold} "
        "--threads={threads} "
        "> {log} 2>&1 "


if bins_for_sample:

    rule move_dastool_output:
        input:
            contig2bin=rules.dastool_run.output.contig2bin,
            summary=rules.dastool_run.output.summary,
        output:
            contig2bin="results/{project}/output/contig2bins/{sample}/DASTool_contig2bin.tsv",
            summary="results/{project}/output/report/{sample}/{sample}_DASTool_summary.tsv",
            done=touch("results/{project}/binning/das_tool/{sample}_move.done"),
        threads: 20
        log:
            "logs/{project}/das_tool/{sample}/move_output.log",
        conda:
            "../envs/unix.yaml"
        shell:
            "(cp {input.contig2bin} {output.contig2bin} && "
            "cp {input.summary} {output.summary}) > {log} 2>&1 "

    rule gzip_bins:
        input:
            bins=rules.dastool_run.output.bins,
        output:
            bins=directory("results/{project}/output/fastas/{sample}/bins/"),
            done=touch("results/{project}/binning/das_tool/{sample}_bins.done"),
        threads: 15
        log:
            "logs/{project}/bins/{sample}/gz_bins.log",
        conda:
            "../envs/unix.yaml"
        shell:
            "(gzip -k {input.bins}/*.fa && "
            "mkdir -p {output.bins}/ && "
            "mv {input.bins}/*.fa.gz {output.bins}/ ) > {log} 2>&1"

    rule move_MAGs:
        input:
            bin_dir=rules.gzip_bins.output.bins,
            csv="results/{project}/output/report/{sample}/{sample}_mags_summary.csv",
        output:
            outdir=directory("results/{project}/output/fastas/{sample}/mags/"),
        threads: 20
        log:
            "logs/{project}/bins/{sample}/move_MAGs.log",
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/move_MAGs.py"
