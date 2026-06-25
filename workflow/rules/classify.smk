# if there is no local database version to use, it is downloaded
if not config["kaiju"]["use-shared"]:

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


# classification of reads
rule run_kaiju_reads:
    input:
        db_files=get_kaiju_files(),
        fastqs=get_filtered_gz_fastqs,
    output:
        report=temp("results/{project}/output/classification/reads/{sample}/kaiju.out"),
    params:
        evalue=0.00001,
    threads: 60
    log:
        "logs/{project}/kaiju/{sample}_reads_run.log",
    conda:
        "../envs/kaiju.yaml"
    shell:
        "kaiju -z {threads} -t {input.db_files[0]} -f {input.db_files[1]} "
        "-v -E {params.evalue} -i {input.fastqs[0]} "
        "-j {input.fastqs[1]} -o {output.report} > {log} 2>&1"


rule kaiju_table_reads:
    input:
        report=rules.run_kaiju_reads.output.report,
        db_files=get_kaiju_files(),
    output:
        summary="results/{project}/output/classification/reads/{sample}/{sample}_{level}_abundance.tsv",
    threads: 2
    log:
        "logs/{project}/kaiju/{sample}_reads_{level}_table.log",
    conda:
        "../envs/kaiju.yaml"
    shell:
        "(kaiju2table -t {input.db_files[0]} -n {input.db_files[2]} "
        "-r {wildcards.level} -o {output.summary} "
        "{input.report}) > {log} 2>&1"
        # -p to print full taxonomy


rule kaiju2krona_reads:
    input:
        db_files=get_kaiju_files(),
        report=rules.run_kaiju_reads.output.report,
    output:
        krona=temp(
            "results/{project}/output/classification/reads/{sample}/kaiju.out.krona"
        ),
        html="results/{project}/output/report/{sample}/{sample}_kaiju_krona_reads.html",
    log:
        "logs/{project}/kaiju/{sample}_reads_2krona.log",
    conda:
        "../envs/kaiju.yaml"
    shell:
        "(kaiju2krona -t {input.db_files[0]} -n {input.db_files[2]} "
        "-i {input.report} -o {output.krona} && "
        "ktImportText -o {output.html} {output.krona}) > {log} 2>&1"


# classification of contigs
rule run_kaiju_contigs:
    input:
        db_files=get_kaiju_files(),
        assembly=get_gz_assembly,
    output:
        report=temp(
            "results/{project}/output/classification/assembly/{sample}/kaiju.out"
        ),
    threads: 60
    log:
        "logs/{project}/contig_classification/{sample}_kaiju_run.log",
    conda:
        "../envs/kaiju.yaml"
    shell:
        "kaiju -z {threads} -t {input.db_files[0]} -f {input.db_files[1]} "
        "-i {input.assembly} -o {output.report} > {log} 2>&1"


rule gzip_kaiju_contigs:
    input:
        report=rules.run_kaiju_contigs.output.report,
    output:
        report="results/{project}/output/classification/assembly/{sample}/kaiju.out.gz",
    threads: 16
    priority: 1
    log:
        "logs/{project}/contig_classification/{sample}_kaiju_gzip.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "gzip -k {input.report} > {log} 2>&1"


use rule kaiju_table_reads as kaiju_table_contigs with:
    input:
        report=rules.run_kaiju_contigs.output.report,
        db_files=get_kaiju_files(),
    output:
        summary="results/{project}/output/classification/assembly/{sample}/{sample}_{level}_abundance.tsv",
    threads: 2
    log:
        "logs/{project}/contig_classification/{sample}_{level}_table.log",


use rule kaiju2krona_reads as kaiju2krona_contigs with:
    input:
        report=rules.run_kaiju_contigs.output.report,
        db_files=get_kaiju_files(),
    output:
        krona=temp(
            "results/{project}/output/classification/assembly/{sample}/kaiju.out.krona"
        ),
        html="results/{project}/output/report/{sample}/{sample}_kaiju_krona_contigs.html",
    log:
        "logs/{project}/contig_classification/{sample}_2krona.log",


rule extract_contig_headers:
    input:
        assembly=get_gz_assembly,
    output:
        contig_len="results/{project}/output/classification/assembly/{sample}/{sample}_contigs_length.csv",
    threads: 4
    log:
        "logs/{project}/contig_classification/{sample}_contig_length.log",
    conda:
        "../envs/unix.yaml"
    shell:
        r"""
        (
            gzip -cd {input.assembly} | \
            awk '
                BEGIN{{print "contig,contig_len"}}
                /^>/ {{
                    sub(/^>/,"",$1);
                    for(i=1;i<=NF;i++)
                        if($i~/^len=/){{
                            split($i,a,"=");
                            print $1","a[2]
                        }}
                }}
            '
        ) > {output.contig_len} 2> {log}
        """


rule contig_classification:
    input:
        out=rules.run_kaiju_contigs.output.report,
        genus="results/{project}/output/classification/assembly/{sample}/{sample}_genus_abundance.tsv",
        species="results/{project}/output/classification/assembly/{sample}/{sample}_species_abundance.tsv",
        contig_len=rules.extract_contig_headers.output.contig_len,
    output:
        csv="results/{project}/output/classification/assembly/{sample}/{sample}_classified_contigs.csv",
    log:
        "logs/{project}/contig_classification/{sample}_classified_contigs.log",
    threads: 4
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/contig_classification.py"


rule gtdbtk_classify_wf:
    input:
        bins=rules.gzip_bins.output.bins,
        #db_prep=rules.prepare_gtdb.output.done,
    output:
        json="results/{project}/output/classification/bins/{sample}/gtdbtk.json",
        outdir=temp(
            directory("results/{project}/output/classification/bins/{sample}/gtdbtk/")
        ),
    params:
        clf_outdir=lambda wildcards, output: Path(output.json).parent,
        json=lambda wildcards, output: Path(output.json).name,
        db_folder=get_gtdb_folder(),
    threads: 32
    # to ensure only one of these commands are run at the same time
    resources:
        heavy=1,
    log:
        "logs/{project}/gtdbtk/{sample}_classify.log",
    conda:
        "../envs/gtdbtk.yaml"
    shell:
        r"""
        (
            export GTDBTK_DATA_PATH="{params.db_folder}" &&

            gtdbtk classify_wf \
                --prefix {wildcards.sample} \
                -x fa.gz \
                --cpus {threads} \
                --pplacer_cpus {threads} \
                --genome_dir {input.bins}/ \
                --out_dir {output.outdir}/ &&
            
            cp {output.outdir}/{params.json} {params.clf_outdir}/
        ) > {log} 2>&1
        """


# combines bacterial and archaeal
rule gtdb_summary:
    input:
        json=rules.gtdbtk_classify_wf.output.json,
        gtdb_dir=rules.gtdbtk_classify_wf.output.outdir,
    output:
        summary="results/{project}/output/classification/bins/{sample}/{sample}.summary.tsv",
    params:
        outdir=lambda wildcards, input: Path(input.json).parent,
    threads: 4
    log:
        "logs/{project}/gtdbtk/{sample}_summary.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/gtdb_summary.py"


rule prep_ncbi_database:
    output:
        temp(touch("results/{project}/ncbi_database_prep.done")),
    log:
        "logs/{project}/tax_abundance/ncbi_prep.log",
    params: db_prep=True,
    conda:
        "../envs/python.yaml"
    threads: 5
    script:
        "../scripts/ncbi_abundance.py"


rule prep_gtdb_database:
    output:
        temp(touch("results/{project}/gtdb_database_prep.done")),
    log:
        "logs/{project}/tax_abundance/gtdb_prep.log",
    params: db_prep=True,
    conda:
        "../envs/python.yaml"
    threads: 5
    script:
        "../scripts/gtdb_abundance.py"

        
rule tax_abundance_ncbi:
    input:
        tsv="results/{project}/output/classification/{stage}/{sample}/{sample}_{level}_abundance.tsv",
        db=rules.prep_ncbi_database.output,
    output:
        csv="results/{project}/output/classification/{stage}/{sample}/{sample}_{stage}_{level}_abundance.csv",
    log:
        "logs/{project}/tax_abundance/{sample}/{stage}_{level}.log",
    params: db_prep=False,
    conda:
        "../envs/python.yaml"
    threads: 5
    script:
        "../scripts/ncbi_abundance.py"


rule tax_abundance_gtdb:
    input:
        csv="results/{project}/output/report/{sample}/{sample}_summary_{stage}.csv",
        db=rules.prep_gtdb_database.output,
    output:
        csv="results/{project}/output/classification/bins/{sample}/{sample}_{stage}_abundance.csv",
    log:
        "logs/{project}/tax_abundance/{sample}/{stage}.log",
    params: db_prep=False,
    conda:
        "../envs/python.yaml"
    threads: 5
    script:
        "../scripts/gtdb_abundance.py"
