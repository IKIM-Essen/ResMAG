import pandas as pd
from ete4 import GTDBTaxa

sys.stderr = open(snakemake.log[0], "w")

only_db_prep = snakemake.params.db_prep

if only_db_prep:
    print("Preparing of ete4 GTDB database..")
    gtdb = GTDBTaxa()
    gtdb.update_taxonomy_database()
    print("done")

else:
    in_file = snakemake.input.csv
    out_file = snakemake.output.csv

    gtdb = GTDBTaxa()

    rank_prefix = {
        "domain": "d__",
        "phylum": "p__",
        "class": "c__",
        "order": "o__",
        "family": "f__",
        "genus": "g__",
        "species": "s__",
    }

    """
    Go from bin/MAG summary (based on GTDB-Tk output) format to a csv with header:
    taxon_id_ncbi,taxon_name_ncbi,absolute_abundance,relative_abundance

    !GTDB and NCBI indeces are incompatible!
    """

    # transfer GTDB taxonomy string to GTDB index
    def str_to_taxid_gtdb(tax_string):
        try:
            for part in reversed(tax_string.split(";")):

                part = part.strip()

                if "__" not in part:
                    continue

                name = part  # IMPORTANT: keep full rank prefix (s__, g__, etc.)

                if not name:
                    continue

                result = gtdb._get_name_translator([name])

                if name in result:
                    taxid = str(result[name][0])
                    return taxid

            return "NaN"

        except Exception:
            return "NaN"

    def bins_to_standard_abundance_table(
        input_csv,
        taxonomy_col="classification",
    ):

        df = pd.read_csv(input_csv)

        # Collapse identical taxonomy strings
        result = (
            df.groupby(taxonomy_col)
            .size()
            .reset_index(name="absolute_abundance")
            .rename(columns={taxonomy_col: "taxon_name_gtdb"})
        )

        # Relative abundance (% of bins)
        total_bins = result["absolute_abundance"].sum()

        result["relative_abundance"] = result["absolute_abundance"] / total_bins * 100

        # Convert taxonomy strings to GTDB taxids
        result["taxon_id_gtdb"] = result["taxon_name_gtdb"].apply(
            lambda x: str_to_taxid_gtdb(x)
        )

        # Reorder columns
        result = result[
            [
                "taxon_id_gtdb",
                "taxon_name_gtdb",
                "absolute_abundance",
                "relative_abundance",
            ]
        ]

        # Sort by abundance
        result = result.sort_values("absolute_abundance", ascending=False)

        return result

    abundance_df = bins_to_standard_abundance_table(in_file)
    abundance_df.to_csv(out_file, index=False)
