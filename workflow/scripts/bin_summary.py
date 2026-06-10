import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

## input files
in_dastool = snakemake.input.tool
in_checkm = snakemake.input.checkm
in_gtdb = snakemake.input.gtdb
sample = snakemake.wildcards.sample

## output files
### csv
csv_path_mags = snakemake.output.csv_mags
csv_path_bins = snakemake.output.csv_bins

## params
max_cont = snakemake.params.max_cont
min_comp = snakemake.params.min_comp

final_columns = [
    "completeness",
    "contamination",
    "bin_score",
    "contigs",
    "genome_size",
    "contig_N50",
    "GC_content",
    "classification",
    "closest_genome_reference",
]

# read in binning information
tool_df = pd.read_table(
    in_dastool,
    usecols=[
        "bin",
        "bin_score",
        "contigs",
    ],
)
tool_df.set_index("bin", inplace=True)


# read in quality information
checkm_df = pd.read_table(
    in_checkm,
    usecols=[
        "Name",
        "Completeness",
        "Contamination",
        "Genome_Size",
        "GC_Content",
        "Contig_N50",
    ],
)

checkm_df.columns = checkm_df.columns.str.lower()
checkm_df.rename(
    {"name": "bin", "gc_content": "GC_content", "contig_n50": "contig_N50"},
    axis=1,
    inplace=True,
)
checkm_df["bin"] = checkm_df["bin"].str.replace(".fa", "", 1, case=True)
checkm_df.set_index("bin", inplace=True)


# read in classification information
gtdb_df = pd.read_table(
    in_gtdb,
    usecols=[
        "user_genome",
        "classification",
        "closest_genome_reference",
    ],
)
gtdb_df.set_index("user_genome", inplace=True)
gtdb_df.index.name = "bin"

# combine information for a bin summary
bins_df = pd.concat(
    [
        tool_df,
        checkm_df,
        gtdb_df,
    ],
    axis=1,
)
# reorder columns
bins_df = bins_df[final_columns]


# add sample information column
bins_df.index.name = "bin"
bins_df = bins_df.reset_index()
bins_df.insert(0, "sample", sample)

# save to csv
bins_df.to_csv(csv_path_bins, index=False)

# filter with MAG criteria
mags_df = bins_df.loc[
    (bins_df["completeness"] >= min_comp) & (bins_df["contamination"] <= max_cont), :
]
mags_df.rename(
    {"bin": "MAG"},
    axis=1,
    inplace=True,
)

# save to csv
mags_df.to_csv(csv_path_mags, index=False)
