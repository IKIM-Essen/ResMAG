import os
import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")


contig_files = snakemake.input.contig_files
abd_class_files = snakemake.input.abd_class_files
abd_ARGs_files = snakemake.input.abd_ARGs_files


## Combining contig and ARG information for all samples
classes_all_samples = []
names_all_samples = []
for contig_file in contig_files:
    sample_id = os.path.basename(os.path.dirname(contig_file))

    classification_df = pd.read_csv(contig_file)
    abd_class_file = [file for file in abd_class_files if file.find(sample_id) >= 0][0]
    abd_classes_df = pd.read_csv(abd_class_file)
    abd_ARGs_file = [file for file in abd_ARGs_files if file.find(sample_id) >= 0][0]
    abd_names_df = pd.read_csv(abd_ARGs_file)

    merged_classes_df = classification_df.merge(
        abd_classes_df, on="contig", how="inner"
    )
    merged_classes_df.insert(0, "sampleID", sample_id)

    classes_all_samples.append(merged_classes_df)

    merged_names_df = classification_df.merge(abd_names_df, on="contig", how="inner")
    merged_names_df.insert(0, "sampleID", sample_id)

    names_all_samples.append(merged_names_df)


final_classes_df = pd.concat(classes_all_samples, ignore_index=True)
final_classes_df.to_csv(snakemake.output.abd_class_possible, index=False)

# in antibiotic classes per sample file all antibiotic classes in CARD are included
# All that do not appear in at least one sample are removed:
final_classes_df_red = final_classes_df.loc[:, (final_classes_df != 0).any(axis=0)]
final_classes_df_red.to_csv(snakemake.output.abd_class_present, index=False)

final_names_df = pd.concat(names_all_samples, ignore_index=True)
final_names_df.to_csv(snakemake.output.abd_ARGs_possible, index=False)

# in ARG per sample file all ARGs in CARD are included
# All that do not appear in at least one sample are removed:
final_names_df_red = final_names_df.loc[:, (final_names_df != 0).any(axis=0)]
final_names_df_red.to_csv(snakemake.output.abd_ARGs_present, index=False)
