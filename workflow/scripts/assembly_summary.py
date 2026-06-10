import pandas as pd
import sys
import os

sys.stderr = open(snakemake.log[0], "w")

asbl_log = snakemake.input.asbl
txt = snakemake.input.mapped
qc_csv = snakemake.input.qc_csv
mag_file = snakemake.input.csv_mags
bin_file = snakemake.input.csv_bins
sample = snakemake.wildcards.sample

df = pd.read_csv(qc_csv, index_col="sample")

summary_dict = {}

summary_sample_dict = {}
summary_sample_dict["#reads_after_filtering"] = df.loc[sample, "#reads_after_filtering"]

with open(asbl_log, "r") as a_log:
    for line in a_log:
        if line.find("total") >= 0:
            info = line.split(" - ")[-1]
            seps = info.strip().split(", ")
            for sep in seps:
                sep = sep.split()
                if "contigs" in sep:
                    colname = "#contigs"
                    value = int(sep[0])
                elif "total" in sep:
                    colname = "assembly_size(bp)"
                    value = int(sep[1])
                elif "max" in sep:
                    colname = "longest_contig(bp)"
                    value = int(sep[1])
                elif "avg" in sep:
                    colname = "average_contig_size(bp)"
                    value = int(sep[1])
                elif "N50" in sep:
                    colname = "N50(bp)"
                    value = int(sep[1])

                summary_sample_dict[colname] = value

with open(txt) as t:
    mapped = int(t.readline())
    summary_sample_dict["#assembled_reads"] = mapped
    summary_sample_dict["%assembled_reads"] = round(
        ((mapped / summary_sample_dict["#reads_after_filtering"]) * 100), 2
    )

if bin_file.rfind(sample) >= 0:

    bin_df = pd.read_csv(bin_file)
    summary_sample_dict["#bins"] = len(bin_df)

    mag_df = pd.read_csv(mag_file)
    summary_sample_dict["#MAGs"] = len(mag_df)

summary_dict[sample] = summary_sample_dict

summary_df = pd.DataFrame.from_dict(summary_dict, orient="index")
summary_df.index.name = "sample"

# reordering columns of output table
summary_df.insert(3, "#assembled_reads", summary_df.pop("#assembled_reads"))
summary_df.insert(4, "%assembled_reads", summary_df.pop("%assembled_reads"))

summary_df.to_csv(snakemake.output.csv)
