import pandas as pd
from ete4 import NCBITaxa

sys.stderr = open(snakemake.log[0], "w")

in_file = snakemake.input.tsv
out_file = snakemake.output.csv

ncbi = NCBITaxa()

rank_prefix = {
    "domain": "d__",
    "phylum": "p__",
    "class": "c__",
    "order": "o__",
    "family": "f__",
    "genus": "g__",
    "species": "s__",
}
VIRUS_TAXID = 10239

"""
Go from kaiju output format based on NCBI IDs to a csv with header:
taxon_id_ncbi,taxon_name_ncbi,absolute_abundance,relative_abundance

!GTDB and NCBI indeces are incompatible!

One caveat: get_name_translator() can return multiple taxids for the same name because some names are ambiguous in NCBI
"""


# transfer NCBI taxID to GTDB like string
def taxid_to_str_ncbi(taxid):
    try:

        lineage = ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        ranks = ncbi.get_rank(lineage)

        # to always have all ranks present
        result = {r: "" for r in rank_prefix}

        for tid in lineage:
            rank = ranks.get(tid)
            if rank in result:
                result[rank] = names[tid]

        tax_str = ";".join(
            f"{prefix}{result[rank]}" for rank, prefix in rank_prefix.items()
        )

        if not tax_str:
            return None
        return tax_str
    except:
        return None


# add GTDB-like taxon string to abundance
def kaiju_to_standard_abundance(input_tsv, use_percent=True):
    df = pd.read_csv(input_tsv, sep="\t")

    out = []
    unclassified_reads = 0
    unclassified_rel = 0

    total_reads = df["reads"].sum()

    for _, row in df.iterrows():

        taxid = row["taxon_id"]

        reads = int(row["reads"])

        rel = row["percent"] if use_percent else (reads / total_reads) * 1e6

        # handle Virus row seperatly
        if taxid == VIRUS_TAXID:
            if reads > 0:
                out.append(
                    {
                        "taxon_id_ncbi": VIRUS_TAXID,
                        "taxon_name_ncbi": "d__Viruses",
                        "absolute_abundance": reads,
                        "relative_abundance": rel,
                    }
                )
            continue

        tax_string = taxid_to_str_ncbi(taxid)

        if tax_string is None:
            unclassified_reads += reads
            unclassified_rel += rel
            continue

        out.append(
            {
                "taxon_id_ncbi": str(int(taxid)),  # save as integer in df
                "taxon_name_ncbi": tax_string,
                "absolute_abundance": reads,
                "relative_abundance": rel,
            }
        )

    df_out = pd.DataFrame(out)

    df_out = pd.concat(
        [
            df_out,
            pd.DataFrame(
                [
                    {
                        "taxon_id_ncbi": "NaN",
                        "taxon_name_ncbi": "unclassified",
                        "absolute_abundance": unclassified_reads,
                        "relative_abundance": unclassified_rel,
                    }
                ]
            ),
        ],
        ignore_index=True,
    )

    return df_out


result = kaiju_to_standard_abundance(in_file)
result.to_csv(out_file, index=False)
