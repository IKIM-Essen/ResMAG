import sys
from pathlib import Path
import shutil
import os

sys.stderr = open(snakemake.log[0], "w")

# prefix rewrite caused by snakemake-fs
OLD_PREFIX = f"/local/work/{os.environ['USER']}/snakemake-scratch/fs/"
NEW_PREFIX = ""


def map_path(p):
    p = str(p)
    if p.startswith(OLD_PREFIX):
        return NEW_PREFIX + p[len(OLD_PREFIX) :]
    return p


for f in snakemake.input.to_delete:
    p = Path(map_path(f))

    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)  # remove folder + contents
    else:
        p.unlink(missing_ok=True)
