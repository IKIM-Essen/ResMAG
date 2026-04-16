import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

## input files
in_dastool = snakemake.input.tool
in_checkm = snakemake.input.checkm
in_gtdb = snakemake.input.gtdb

## output files
### csv
csv_path_mags = snakemake.output.csv_mags
csv_path_bins = snakemake.output.csv_bins
csv_path_tax = snakemake.output.csv_tax

## params
max_cont = snakemake.params.max_cont
min_comp = snakemake.params.min_comp


def save_csv_table(csv_path, summary_df):
    summary_df.to_csv(csv_path)

    summary_df.columns.name = "bin"
    summary_df.index.name = None


tool_df = pd.read_table(in_dastool)
tool_df.drop(["bin_set"], axis=1, inplace=True)

tool_df.set_index("bin", inplace=True)
col_order = [
    "bin_score",
    "contigs",
    "size",
    "N50",
    "SCG_set",
    "unique_SCGs",
    "SCG_completeness",
    "redundant_SCGs",
    "SCG_redundancy",
]
tool_df = tool_df[col_order]

del col_order[:2]
tool_red_df = tool_df.drop(col_order, axis=1)

checkm_df = pd.read_table(in_checkm)
rm_list = ["Translation_Table_Used", "Additional_Notes"]
checkm_df.drop(rm_list, axis=1, inplace=True)
checkm_df.columns = checkm_df.columns.str.lower()
checkm_df.rename(
    {"name": "bin", "gc_content": "GC_content", "contig_n50": "contig_N50"},
    axis=1,
    inplace=True,
)
checkm_df["bin"] = checkm_df["bin"].str.replace(".fa", "", 1, case=True)
checkm_df.set_index("bin", inplace=True)
col_order = [
    "completeness",
    "contamination",
    "genome_size",
    "GC_content",
    "contig_N50",
    "total_coding_sequences",
    "coding_density",
    "average_gene_length",
    "completeness_model_used",
]
checkm_df = checkm_df[col_order]

del col_order[:5]
checkm_red_df = checkm_df.drop(col_order, axis=1)

gtdb_df = pd.read_table(in_gtdb)
gtdb_df.rename({"user_genome": "bin"}, axis=1, inplace=True)
gtdb_df.set_index("bin", inplace=True)
cols = [
    "classification",
    "closest_genome_reference",
    "closest_genome_reference_radius",
    "closest_genome_ani",
    "closest_genome_af",
    "classification_method",
]
gtdb_red_df = gtdb_df[cols]

save_csv_table(csv_path_tax, gtdb_red_df)


bins_df = pd.concat(
    [
        tool_red_df,
        checkm_red_df,
        gtdb_red_df[
            ["classification", "closest_genome_reference", "closest_genome_ani"]
        ],
    ],
    axis=1,
)
col_order = [
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
bins_df = bins_df[col_order]
bins_df.index.name = "bin"

save_csv_table(csv_path_bins, bins_df)


mags_df = bins_df.loc[
    (bins_df["completeness"] >= min_comp) & (bins_df["contamination"] <= max_cont), :
]
mags_df.index.name = "MAG"
save_csv_table(csv_path_mags, mags_df)
