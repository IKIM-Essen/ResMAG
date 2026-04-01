import os
import sys
from concurrent.futures import ThreadPoolExecutor

sys.stderr = open(snakemake.log[0], "w")

prefix = snakemake.params.prefix
binner = snakemake.params.binner
output_tsv = snakemake.output[0]


def append_string_to_file(filepath, append_string, mod_ext):
    cmd = f"sed 's/$/\t{append_string}/' {filepath} > {filepath}{mod_ext}"
    os.system(cmd)


if __name__ == "__main__":
    # Directory containing the files
    directory = snakemake.input.outdir

    # List all files in the directory that start with the prefix specified in binning process
    files = [
        filename for filename in os.listdir(directory) if filename.startswith(prefix)
    ]

    num_threads = snakemake.threads
    mod_ext = ".modified"

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        # Submit the tasks to the ThreadPoolExecutor
        for filename in files:
            filepath = os.path.join(directory, filename)
            # String to append to each line
            append_string = "{0}_bin_{1}".format(
                binner, filename.split(f"{prefix}.")[-1]
            )

            executor.submit(append_string_to_file, filepath, append_string, mod_ext)

    modfiles = [
        os.path.join(directory, filename)
        for filename in os.listdir(directory)
        if filename.endswith(mod_ext)
    ]
    if len(modfiles) == 0:
        print("no metabat bins")
        # create an empty file
        cmd = f"touch {output_tsv}"
        os.system(cmd)
    else:
        print("starting to concatenate files")
        for modfile in modfiles:
            # Concatenate all modified files into a single file
            cmd = f"cat {modfile} >> {output_tsv}"
            os.system(cmd)

        print("all files concatenated")
        # remove all modified files after concatenating
        cmd = f"rm {directory}/*{mod_ext}"
        os.system(cmd)
