
rule coverm_metabat:
    input:
        bact_reads=get_filtered_gz_fastqs,
        contigs=get_gz_assembly,
    output:
        abd=temp("results/{project}/binning_prep/{sample}/abundance_metabat.tsv"),
    threads: 30
    log:
        "logs/{project}/coverm_metabat/{sample}.log",
    conda:
        "../envs/coverm.yaml"
    shell:
        "coverm contig -m metabat -1 {input.bact_reads[0]} "
        "-2 {input.bact_reads[1]} -r {input.contigs} "
        "-t {threads} -o {output.abd} > {log} 2>&1"


rule metabat:
    input:
        contigs=get_gz_assembly,
        abd=rules.coverm_metabat.output.abd,
    output:
        outdir=temp(directory("results/{project}/binning/metabat/{sample}/")),
    threads: 60
    params:
        prefix="bin",
    log:
        "logs/{project}/metabat/{sample}.log",
    conda:
        "../envs/metabat.yaml"
    shell:
        "metabat1 -i {input.contigs} -a {input.abd} "
        "--seed 1 --p1 95 --p2 90 "
        "-B 20 --pB 50 --minProb 80 --minBinned 20 "
        "-m 1500 -s 100000 "
        "--minCorr 95 --minContigByCorr 300 --minSamples 1 "
        "-t {threads} -v -l -o {output}/{params.prefix} > {log} 2>&1"
