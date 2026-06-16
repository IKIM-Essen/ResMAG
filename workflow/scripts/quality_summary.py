import pandas as pd
import json
import sys

sys.stderr = open(snakemake.log[0], "w")

json_file = snakemake.input.json
human_log = snakemake.input.human_log

outfile = snakemake.output.csv

results_dict = {}

# for json_file in jsons:
sample = json_file[json_file.rfind("/") + 1 : json_file.rfind(".fastp.json")]
sample_results_dict = {}

with open(json_file, "r") as read_file:
    jdata = json.load(read_file)
    # fastp json shows reads, not read pairs
    reads_before = int(jdata["summary"]["before_filtering"]["total_reads"])
    sample_results_dict["#reads"] = reads_before

    total_pre_host_filt = int(jdata["summary"]["after_filtering"]["total_reads"])
    sample_results_dict["#reads_afterQC"] = total_pre_host_filt

    sample_results_dict["%reads_filtered_by_quality"] = round(
        ((1 - (total_pre_host_filt / reads_before)) * 100), 2
    )

    sample_results_dict["#bp_afterQC"] = int(
        jdata["summary"]["after_filtering"]["total_bases"]
    )

    sample_results_dict["%Q30_reads_afterQC"] = round(
        (float(jdata["summary"]["after_filtering"]["q30_rate"]) * 100), 2
    )

with open(human_log, "r") as file:
    for line in file.readlines():
        if line.find("processed") >= 0:
            non_human_reads = int(line.split()[2]) * 2
            no_reads = total_pre_host_filt - non_human_reads
            prct_reads = round(((int(no_reads) / int(total_pre_host_filt)) * 100), 2)
            break

sample_results_dict[f"#human_reads"] = no_reads
sample_results_dict[f"%human"] = prct_reads

sample_results_dict["#reads_after_filtering"] = non_human_reads
sample_results_dict["%reads_filtered_total"] = round(
    ((1 - (non_human_reads / reads_before)) * 100), 2
)

results_dict[sample] = sample_results_dict

results_df = pd.DataFrame.from_dict(results_dict, orient="index")
results_df.index.name = "sample"

results_df = results_df.sort_index()
results_df.to_csv(outfile)
