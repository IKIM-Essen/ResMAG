import pandas as pd
import json
import sys
import re

sys.stderr = open(snakemake.log[0], "w")

card_drugs = snakemake.input.json
args_file = snakemake.input.csv
contig_file = snakemake.input.contig_file
contig_len_file = snakemake.input.contig_len
sample_id = snakemake.wildcards.sample

# abundance outfiles
abd_class = snakemake.output.abd_class
abd_ARGs = snakemake.output.abd_ARGs

# Load CARD database drug class information
with open(card_drugs) as f:
    data = json.load(f)

# Convert to DataFrame
df_card = pd.DataFrame.from_dict(data, orient="index")

# Split semicolon-separated values into lists
df_card["classes"] = df_card["classes"].str.split("; ")
df_card["antibiotics"] = df_card["antibiotics"].str.split("; ")

# Flatten and get unique values for drug classes and antibiotics
all_classes = pd.Series(
    [item for sublist in df_card["classes"].dropna() for item in sublist]
).unique()
all_classes = sorted(all_classes)

all_antibiotics = pd.Series(
    [item for sublist in df_card["antibiotics"].dropna() for item in sublist if item]
).unique()
all_aros = df_card["ARO_name"].unique()
all_aros = sorted(all_aros)

# ---------------------------------------------------------
# Create abundance table for ARG drug classes per contig
# ---------------------------------------------------------

# Read contig lengths
contig_len_df = pd.read_csv(contig_len_file)

# read in contig classification file
classification_df = pd.read_csv(contig_file)
classification_df.drop(["contig_len"], axis=1, inplace=True)

# Read in resistance information
df_args = pd.read_csv(args_file)

# If no ARGs were detected, write empty output files with correct headers
if df_args.empty:
    # ----- drug class output -----
    class_cols = [
        "sample",
        "contig",
        "contig_len",
        "level",
        "classification",
        "#ARGs",
    ] + all_classes

    pd.DataFrame(columns=class_cols).to_csv(abd_class, index=False)

    # ----- ARG output -----
    aro_cols = [
        "sample",
        "contig",
        "contig_len",
        "level",
        "classification",
        "#ARGs",
    ] + all_aros

    pd.DataFrame(columns=aro_cols).to_csv(abd_ARGs, index=False)

    sys.exit(0)

# Make a copy for drug classes to avoid modifying original df
df_args_for_classes = df_args.copy()

# Convert semicolon-separated class strings into lists
# Example: "aminoglycoside antibiotic;carbapenem" -> ["aminoglycoside antibiotic", "carbapenem"]
df_args_for_classes["class_list"] = (
    df_args_for_classes["classes"]
    .fillna("")
    .apply(lambda x: [c.strip() for c in x.split(";")])
)

# Create one-hot encoded columns for each drug class
# (1 = class present for this ARG, 0 = absent)
for c in all_classes:
    df_args_for_classes[c] = df_args_for_classes["class_list"].apply(
        lambda lst: int(c in lst)
    )

# Aggregate counts per contig
# "#ARGs" = number of ARGs (ARO_IDs) on the contig
# Drug class columns are summed to give abundance per class
agg_dict = {"ARO_ID": "count"}
agg_dict.update(
    {c: "sum" for c in all_classes}  # "sum" counts how many times the class appears
)

class_abundance_df = (
    df_args_for_classes.groupby("contig")
    .agg(agg_dict)
    .reset_index()
    .rename(columns={"ARO_ID": "#ARGs"})
)

# Merge with taxonomic classification
# Use LEFT JOIN to keep all contigs containing ARGs
merged_classes_df = class_abundance_df.merge(classification_df, on="contig", how="left")

# Add contig lengths
merged_classes_df = merged_classes_df.merge(contig_len_df, on="contig", how="left")

# Fill missing taxonomy information
# Contigs without classification are labelled as unclassified
merged_classes_df["level"] = merged_classes_df["level"].fillna("unclassified")

merged_classes_df["classification"] = merged_classes_df["classification"].fillna(
    "unclassified"
)

# Reorder columns
# Sort drug classes
class_columns = sorted([a for a in all_classes if a in merged_classes_df.columns])
cols = [
    "contig",
    "contig_len",
    "level",
    "classification",
    "#ARGs",
] + class_columns
merged_classes_df = merged_classes_df[cols]
merged_classes_df.insert(0, "sample", sample_id)

# Write out class abundance file
merged_classes_df.to_csv(abd_class, index=False)


# ---------------------------------------------------------
# Create ARG abundance table per contig
# ---------------------------------------------------------

"""# Add a column indicating presence of an ARG hit
df_args["present"] = 1

# Create a contig × ARG matrix
# Values correspond to the number of occurrences of each ARG
aro_matrix = df_args.pivot_table(
    index="contig",
    columns="ARO_name",
    values="present",
    aggfunc="sum",  # changed from "max" to "sum"
    fill_value=0,
)"""
# Create a contig × ARG matrix
# Values correspond to the number of occurrences of each ARG
aro_matrix = df_args.groupby(["contig", "ARO_name"]).size().unstack(fill_value=0)


# Count the total number of ARGs per contig
gene_counts = df_args.groupby("contig").size().rename("#ARGs")

# Combine total ARG counts and ARG abundance matrix
names_abundance_df = pd.concat([gene_counts, aro_matrix], axis=1).reset_index()

# Add columns for ARGs that are absent in this sample
# so all output files have identical columns
missing_aros = [aro for aro in all_aros if aro not in names_abundance_df.columns]

if missing_aros:
    missing_df = pd.DataFrame(0, index=names_abundance_df.index, columns=missing_aros)

    names_abundance_df = pd.concat([names_abundance_df, missing_df], axis=1)

# Merge with taxonomic classification
# Use LEFT JOIN to keep all contigs containing ARGs
merged_names_df = names_abundance_df.merge(classification_df, on="contig", how="left")

# Add contig lengths
merged_names_df = merged_names_df.merge(contig_len_df, on="contig", how="left")

# Fill missing taxonomy information
# Contigs without classification are labelled as unclassified
merged_names_df["level"] = merged_names_df["level"].fillna("unclassified")

merged_names_df["classification"] = merged_names_df["classification"].fillna(
    "unclassified"
)

# Reorder columns
# Sort ARGs
aro_columns = sorted([a for a in all_aros if a in merged_names_df.columns])
cols = [
    "contig",
    "contig_len",
    "level",
    "classification",
    "#ARGs",
] + aro_columns
merged_names_df = merged_names_df[cols]
merged_names_df.insert(0, "sample", sample_id)

# Write out ARG abundance file
merged_names_df.to_csv(abd_ARGs, index=False)
