import pandas as pd
import sys
import os

sys.stderr = open(snakemake.log[0], "w")

asbl_logs = snakemake.input.asbl
txts = snakemake.input.mapped
qc_csv = snakemake.input.qc_csv
csv_mags = snakemake.input.csv_mags
csv_bins = snakemake.input.csv_bins

df=pd.read_csv(qc_csv,index_col="sample")
samples=df.index.to_list()

summary_dict = {}
for sample in samples:
    summary_sample_dict = {}
    summary_sample_dict["#reads_after_filtering"] = df.loc[sample,"#reads_after_filtering"]

    for asbl_log in asbl_logs:
        if asbl_log.rfind(f'{sample}_megahit') >= 0:
            with open(asbl_log, "r") as a_log:
                for line in a_log:
                    if line.find("total") >= 0:
                        info = line.split(" - ")[-1]
                        seps = info.strip().split(", ")
                        for sep in seps:
                            sep = sep.split()
                            if "contigs" in sep:
                                colname = "#contigs"
                                value=int(sep[0])
                            elif "total" in sep:
                                colname = "assembly_size(bp)"
                                value=int(sep[1])
                            elif "max" in sep:
                                colname = "longest_contig(bp)"
                                value=int(sep[1])
                            elif "avg" in sep:
                                colname = "average_contig_size(bp)"
                                value=int(sep[1])
                            elif "N50" in sep:
                                colname = "N50(bp)"
                                value=int(sep[1])
                            
                            summary_sample_dict[colname] = value

    for txt in txts:
        if txt.rfind(f'{sample}_reads_mapped') >= 0:
            with open(txt) as t:
                mapped=int(t.readline())
                summary_sample_dict["#assembled_reads"] = mapped
                summary_sample_dict["%assembled_reads"] = round(((mapped/summary_sample_dict["#reads_after_filtering"]) * 100),2)
            break
    
    for bin_file in csv_bins:
        if bin_file.rfind(sample) >= 0:

            bin_df=pd.read_csv(bin_file)
            summary_sample_dict["#bins"] = len(bin_df)

            mag_file=[file for file in csv_mags if os.path.basename(os.path.dirname(file)) == sample][0]
            mag_df=pd.read_csv(mag_file)
            summary_sample_dict["#MAGs"] = len(mag_df)

    summary_dict[sample] = summary_sample_dict

summary_df = pd.DataFrame.from_dict(summary_dict, orient="index")
summary_df.index.name = "sample"

# reordering columns of output table
summary_df.insert(3, "#assembled_reads", summary_df.pop("#assembled_reads"))
summary_df.insert(4, "%assembled_reads", summary_df.pop("%assembled_reads"))

summary_df.to_csv(snakemake.output.csv)

## edit decimal separator for visual of table
def format_int_with_commas(x):
    #Formats an integer with commas as thousand separators
    return f"{x:,}"

def format_float_with_percent(x):
    #Formats a float with percent sign
    return f"{x:.2f}%"

cols=[col for col in summary_df.columns.to_list() if col.find("%")>=0]
summary_df[cols]=summary_df[cols].map(format_float_with_percent)

cols=[col for col in summary_df.columns.to_list() if col.find("%")<0]
summary_df[cols]=summary_df[cols].map(format_int_with_commas)

summary_df.to_csv(snakemake.output.vis_csv)
