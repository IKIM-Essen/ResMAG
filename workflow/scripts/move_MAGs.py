import pandas as pd
import sys
import os

sys.stderr = open(snakemake.log[0], "w")


summary = snakemake.input.csv
bin_dir = snakemake.input.bin_dir
outdir = snakemake.output.outdir

mags_df = pd.read_csv(summary, index_col="MAG")
mags = mags_df.index.to_list()

os.system(f"mkdir -p {outdir}")
if len(mags) > 0:
    for mag in mags:
        cmd = f"mv {bin_dir}/{mag}.fa.gz {outdir}/"
        os.system(cmd)
