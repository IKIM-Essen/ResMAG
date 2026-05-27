import os


configfile: "config/config.yaml"


def get_project():
    return config["project-name"]


def get_samples():
    return list(pep.sample_table["sample_name"].values)


def get_fastqs(wildcards):
    file_r1 = pep.sample_table.loc[wildcards.sample]["fq1"]
    file_r2 = pep.sample_table.loc[wildcards.sample]["fq2"]
    return [
        file_r1,
        file_r2,
    ]


def get_adapters(wildcards):
    return config["adapter-seqs"]


def get_trimmed_fastqs(wildcards):
    return [
        "results/{project}/trimmed/fastp/{sample}.1.fastq.gz",
        "results/{project}/trimmed/fastp/{sample}.2.fastq.gz",
    ]


def get_prefiltered_fastqs(wildcards):
    if config["host-filtering"]["do-host-filtering"]:
        return [
            "results/{project}/host_filtering/non_host/{sample}_R1.fastq.gz",
            "results/{project}/host_filtering/non_host/{sample}_R2.fastq.gz",
        ]
    else:
        return get_trimmed_fastqs(wildcards)


def get_host_map_statistics(wildcards):
    if config["host-filtering"]["do-host-filtering"]:
        logs = expand(
            "results/{{project}}/output/report/prerequisites/qc/filter_host_{sample}.log",
            sample=get_samples(),
        )
        return logs
    else:
        return []


def get_filtered_gz_fastqs(wildcards):
    return [
        "results/{project}/output/filtered_reads/{sample}_R1.fastq.gz",
        "results/{project}/output/filtered_reads/{sample}_R2.fastq.gz",
    ]


def get_assembly(wildcards):
    return "results/{project}/megahit/{sample}/final.contigs.fa"


def get_gz_assembly(wildcards):
    return "results/{project}/output/fastas/{sample}/{sample}.fa.gz"


def get_proteins(wildcards):
    return ("results/{project}/output/proteins/{sample}/{sample}_proteins.faa.gz",)


## binning parameters
def get_contig_length_threshold():
    return config["min-contig-length"]


def bins_for_sample(wildcards):
    binfolder = "results/{project}/binning/das_tool/{sample}/{sample}_DASTool_bins/"
    bin_list = [b for b in os.listdir(binfolder) if b.endswith(".fa")]
    if len(bin_list) > 0:
        return True
    else:
        return False


def get_human_ref_download():
    return config["human-filtering"]["download-path"]


def get_human_ref():
    if config["human-filtering"]["use-shared"]:
        path = config["human-filtering"]["local-path"]
        return path
    else:
        path = get_human_ref_download()
        filename = path.split("/")[-1]
        local_ref = "".join(["resources/ref_genome/", filename])
        return local_ref


def get_host_ref():
    return config["host-filtering"]["ref-genome"]


def get_kaiju_files():
    file = config["kaiju"]["fmi-file"]
    path = str(Path(file).parent)
    fmi = Path(file).name
    names = ["nodes.dmp", fmi, "names.dmp"]
    files = ["/".join([path, name]) for name in names]
    return files


def get_checkm2_db():
    path = "".join(
        [config["checkm2"]["db-folder"], "CheckM2_database/uniref100.KO.1.dmnd"]
    )
    return path


def get_gtdb_folder():
    return config["gtdb"]


def get_genomad_DB_folder():
    path = "".join([config["genomad"]["db-folder"], "genomad_db/"])
    return path


def get_genomad_DB_file():
    path = "".join([get_genomad_DB_folder(), "names.dmp"])
    return path


def get_card_db_file():
    return config["card"]["db-file"]


def get_card_annotation_file():
    version = config["card"]["version"]
    filename = "".join(["card_database_", version, ".fasta"])
    folder = str(Path(get_card_db_file()).parent)

    path = "/".join([folder, filename])
    return path


def get_uniCARD_db():
    return config["uniCARD"]["db-file"]


def get_CARD_hierarchy():
    return config["uniCARD"]["hierarchy-json"]


def get_uniCARD_db_wo_ext():
    db_wo_ext = str(Path(get_uniCARD_db()).with_suffix("").with_suffix(""))
    return db_wo_ext


def get_unicard_dmnd():
    return "".join([get_uniCARD_db_wo_ext(), ".dmnd"])


"""

def get_mag_fa(wildcards):
    folder = "results/{}/output/fastas/{}/mags/".format(
        wildcards.project, wildcards.sample
    )
    files = [
        os.path.join(folder, binID)
        for binID in os.listdir(folder)
        if binID.endswith("fa.gz")
    ]
    return files


def get_mag_ARGs(wildcards):
    bin_fastas = (get_mag_fa(wildcards),)
    bin_fastas = list(bin_fastas)[0]
    binIDs = [os.path.basename(binID) for binID in bin_fastas]
    folder = "results/{}/output/ARGs/mags/{}/".format(
        wildcards.project, wildcards.sample
    )
    arg_files = [
        os.path.join(folder, binID.replace(".fa.gz", ".txt")) for binID in binIDs
    ]
    return arg_files
"""
