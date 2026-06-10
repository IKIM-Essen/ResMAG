import pandas as pd
import json
import sys

sys.stderr = open(snakemake.log[0], "w")

card_drugs = snakemake.input.json
args_file = snakemake.input.csv
contig_file = snakemake.input.contig_file
sample_id = snakemake.wildcards.sample

# abundance outfiles
abd_class = snakemake.output.abd_class
abd_ARGs = snakemake.output.abd_ARGs

## Write drug classes and ARO names abundance:
with open(card_drugs) as f:
    data = json.load(f)

# Convert to DataFrame
df_card = pd.DataFrame.from_dict(data, orient="index")

# Split semicolon-separated values into lists
df_card["classes"] = df_card["classes"].str.split("; ")
df_card["antibiotics"] = df_card["antibiotics"].str.split("; ")

# Flatten and get unique values
all_classes = pd.Series(
    [item for sublist in df_card["classes"].dropna() for item in sublist]
).unique()
all_classes.sort()
all_antibiotics = pd.Series(
    [item for sublist in df_card["antibiotics"].dropna() for item in sublist if item]
).unique()
all_aros = df_card["ARO_name"].unique()
all_aros.sort()

# read in contig classification
classification_df = pd.read_csv(contig_file)

##Create classes and ARG names abundance files
df_args = pd.read_csv(args_file)
df_args_for_classes = df_args.copy()
# Step 1: Preprocess the classes column into lists
df_args_for_classes["class_list"] = (
    df_args_for_classes["classes"]
    .fillna("")
    .apply(lambda x: [c.strip() for c in x.split(";")])
)

# Step 2: Create one-hot columns per class
for c in all_classes:
    df_args_for_classes[c] = df_args_for_classes["class_list"].apply(
        lambda lst: int(c in lst)
    )

# Step 3: Group by contig
agg_dict = {"ARO_ID": "count"}
agg_dict.update(
    {c: "sum" for c in all_classes}  # was max before
)  # "sum" counts how many times the class appears


class_abundance_df = df_args_for_classes.groupby("contig").agg(agg_dict).reset_index()
# Write out abundance file
class_abundance_df = class_abundance_df.rename(columns={"ARO_ID": "#ARGs"})

merged_classes_df = classification_df.merge(
    class_abundance_df, on="contig", how="inner"
)
merged_classes_df.insert(0, "sampleID", sample_id)

merged_classes_df.to_csv(abd_class, index=False)


df_args["present"] = 1

# Step 2: Pivot to get binary presence/absence for all ARO_name per contig
aro_matrix = df_args.pivot_table(
    index="contig",
    columns="ARO_name",
    values="present",
    aggfunc="sum",  # <-- changed from "max" to "sum"
    fill_value=0,
)

# Step 3: Count ARGs (genes) per contig
gene_counts = df_args.groupby("contig").size().rename("#ARGs")

# Step 4: Join gene counts and presence matrix
names_abundance_df = pd.concat([gene_counts, aro_matrix], axis=1).reset_index()

# Find missing AROs
missing_aros = [aro for aro in all_aros if aro not in names_abundance_df.columns]

# Create a DataFrame with 0s for missing AROs (same index as result)
missing_df = pd.DataFrame(0, index=names_abundance_df.index, columns=missing_aros)

# Concatenate in one go — avoids fragmentation
names_abundance_df = pd.concat([names_abundance_df, missing_df], axis=1)

# Optional: reorder columns — num_genes first, then sorted AROs
aro_columns = sorted([a for a in all_aros if a in names_abundance_df.columns])
names_abundance_df = names_abundance_df[["contig", "#ARGs"] + aro_columns]

# Write out abundance file
merged_names_df = classification_df.merge(names_abundance_df, on="contig", how="inner")
merged_names_df.insert(0, "sampleID", sample_id)

merged_names_df.to_csv(abd_ARGs, index=False)
