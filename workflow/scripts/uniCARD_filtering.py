import json
import pandas as pd
import argparse

# identifier for CARD derived sequences
WT_PREFIX = "CARD_WT:"
CARD_PREFIX = "CARD:"

# columns to reduce DIAMOND output
DIAMOND_COLS = [
    "query",
    "target",
    "%identity",
    "length",
    "#mismatches",
    "#gaps",
    "query_start",
    "query_end",
    "target_start",
    "target_end",
    "evalue",
    "bitscore",
]


def write_filtered_ARGs(card_classes, unicard_in, outfile, filter_window):

    try:
        # read in CARD drugs and classes json file
        with open(card_classes) as f:
            data = json.load(f)

        card_df = pd.DataFrame.from_dict(data, orient="index")
        card_df.index.name = "ARO_ID"
        card_df = card_df.reset_index()

        if unicard_in.endswith(".gz"):
            df_unicard = pd.read_table(
                unicard_in, names=DIAMOND_COLS, compression="gzip"
            )
        else:
            df_unicard = pd.read_table(unicard_in, names=DIAMOND_COLS)

        ## top 10 ##
        # Step 1: Keep only 10 best hits per query
        df_top10 = df_unicard[df_unicard.groupby("query").cumcount() < 10]

        # Step 2: Keep all rows for a query **only if at least one Target starts with CARD_PREFIX,
        # so wild-type CARD entries are not kept (starting with WT_PREFIX)
        df_top10 = df_top10[
            df_top10["query"].isin(
                df_top10[df_top10["target"].str.startswith(CARD_PREFIX)]["query"]
            )
        ]

        ## filter besthit CARD entries

        # filter for the first row per query
        df_firsts = df_top10[df_top10.groupby("query").cumcount() < 1]
        # 1. check if best hit is CARD wild-type -> remove from df
        df_firsts_reduced = df_firsts[
            ~df_firsts["query"].isin(
                df_firsts[df_firsts["target"].str.startswith(WT_PREFIX)]["query"]
            )
        ]
        # 2. from remaining save those, which are CARD hits
        df_bests = df_firsts_reduced[
            df_firsts_reduced["query"].isin(
                df_firsts_reduced[
                    df_firsts_reduced["target"].str.startswith(CARD_PREFIX)
                ]["query"]
            )
        ]

        # remove first hit CARD queries from first row per query df
        df_firsts_remaining = df_firsts_reduced[
            ~df_firsts_reduced["query"].isin(df_bests["query"])
        ]

        df_filtered = df_top10[df_top10["query"].isin(df_firsts_remaining["query"])]

        filtered_card_rows = []
        # go through all remaining best hits and compare to CARD hits
        for idx, row in df_firsts_remaining.iterrows():
            # get all rows with same query to compare
            compare_df = df_filtered[df_filtered["query"] == row["query"]]
            compare_df = compare_df[
                compare_df["target"].str.startswith(CARD_PREFIX)
                & ~compare_df["target"].str.contains("_snp_")
            ]
            # check for each row in this comparison df if the CARD entry has same values than first match
            for idx_card, row_card in compare_df.iterrows():
                if (
                    row["%identity"] == row_card["%identity"]
                    and row["bitscore"] == row_card["bitscore"]
                ):
                    filtered_card_rows.append(row_card)
                    break
                # additional filtering in default 1% window below best hit
                elif (row["%identity"] * (1 - filter_window)) <= row_card[
                    "%identity"
                ] and (row["bitscore"] * (1 - filter_window)) <= row_card["bitscore"]:
                    filtered_card_rows.append(row_card)
                    break

        filtered_card_df = pd.DataFrame(filtered_card_rows)

        df_bests = pd.concat([df_bests, filtered_card_df])

        # Adding contig column
        df_bests.insert(
            1, "contig", df_top10["query"].map(lambda x: x.rsplit("_", 1)[0])
        )
        # extract ARO_ID, for variants it's eg 'CARD:3003394_snp_R123S'
        df_bests.insert(
            3,
            "ARO_ID",
            df_top10["target"].map(
                lambda x: x.split(CARD_PREFIX, 1)[-1].split("_", 1)[0]
            ),
        )

        # Adding CARD drug info:
        merged_df = pd.merge(df_bests, card_df, on="ARO_ID", how="left")

        # Write besthits to output file
        merged_df.to_csv(outfile, index=False)

    except Exception as e:
        print(f"Error filtering UniCARD results: {e}")
        raise


def main():
    parser = argparse.ArgumentParser(description="Filter UniCARD results")

    parser.add_argument(
        "--infile",
        required=True,
        help="DIAMOND output tsv file with UniCARD results to filter, may be gz",
    )
    parser.add_argument(
        "--card_hierarchy", required=True, help="Path to CARDs hierarchy json file"
    )
    parser.add_argument(
        "--outfile",
        required=True,
        help="Output csv file to store final filtered ARGs",
    )
    parser.add_argument(
        "--filter_window",
        type=float,
        default=0.01,
        help="Percent window for filtering (default 0.01 (1%%))",
    )

    args = parser.parse_args()

    infile = args.infile
    outfile = args.outfile
    card_hierarchy = args.card_hierarchy
    filter_window = args.filter_window

    print("Filtering UniCARD diamond results...")
    write_filtered_ARGs(card_hierarchy, infile, outfile, filter_window)

    print(f"Filtering done! Results saved to: {outfile}")


if __name__ == "__main__":
    main()
