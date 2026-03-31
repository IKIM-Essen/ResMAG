import pandas as pd
import re
import sys

sys.stderr = open(snakemake.log[0], "w")

kaiju_out = snakemake.input.out
kaiju_species = snakemake.input.species
kaiju_genus = snakemake.input.genus
headers = snakemake.input.headers

classification = snakemake.output.csv


## Prepare classification and length information for all contigs classified to at least genus level ##

df_kaiju_raw = pd.read_table(kaiju_out, names=["classified", "contig", "taxid"])
kaiju_df = df_kaiju_raw[df_kaiju_raw["classified"] != "U"].copy()

species_df = pd.read_csv(kaiju_species, sep="\t")
# remove rows for unclassified and virus
species_df = species_df.iloc[:-3].copy()

genus_df = pd.read_csv(kaiju_genus, sep="\t")
# remove rows for unclassified and virus
genus_df = genus_df.iloc[:-3].copy()

# make taxon name to taxon id dictionary for all occurring species and genera
species_map = dict(zip(species_df["taxon_id"], species_df["taxon_name"]))
genus_map = dict(zip(genus_df["taxon_id"], genus_df["taxon_name"]))

# kaiju_df is the original DataFrame with only classified contigs
kaiju_df["species_name"] = kaiju_df["taxid"].map(species_map)
kaiju_df["genus_name"] = kaiju_df["taxid"].map(genus_map)

# removing all not classified to species or genus level
kaiju_df_red = kaiju_df[
    ~(kaiju_df["species_name"].isna() & kaiju_df["genus_name"].isna())
].copy()
# adding info about the level of classification
kaiju_df_red["level"] = (
    kaiju_df_red["species_name"].notna().map({True: "species", False: "genus"})
)
kaiju_df_red["classification"] = kaiju_df_red["species_name"].combine_first(
    kaiju_df_red["genus_name"]
)
kaiju_df_red = kaiju_df_red.drop(
    columns=["species_name", "genus_name", "classified", "taxid"]
)

# get contig length
records = []

with open(headers) as f:
    for line in f:
        # Example header: >k141_570066 flag=0 multi=1.0000 len=450
        contig_match = re.search(r"^>(\S+)", line)
        length_match = re.search(r"len=(\d+)", line)

        if contig_match and length_match:
            contig = contig_match.group(1)
            length = int(length_match.group(1))
            records.append((contig, length))

# Convert to DataFrame
length_df = pd.DataFrame(records, columns=["contig", "contig_len"])

# merge contig classification with length information
merged_df = kaiju_df_red.merge(length_df, on="contig", how="left")
cols = merged_df.columns.tolist()
cols.insert(1, cols.pop(cols.index("contig_len")))
merged_df = merged_df[cols]

# Write out contig classification file
merged_df.to_csv(classification, index=False)
