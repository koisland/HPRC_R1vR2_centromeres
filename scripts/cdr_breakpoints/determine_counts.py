import sys
import polars as pl

LIFTOVER_COLS = ["qry", "qry_st", "qry_end", "ref", "ref_st", "ref_end"]


def main():
    liftover = sys.argv[1]
    fai = sys.argv[2]
    bed_complete_qry = sys.argv[3]
    bed_complete_ref = sys.argv[4]
    outfile = sys.argv[5]

    # Read in liftover.
    df_liftover = pl.read_csv(
        liftover,
        separator="\t",
        has_header=False,
        columns=range(len(LIFTOVER_COLS)),
        new_columns=LIFTOVER_COLS,
    )
    # Read in query contig (R1) length
    df_fai = pl.read_csv(
        fai,
        separator="\t",
        has_header=False,
        new_columns=["qry", "qry_length"],
        columns=[0, 1],
    )
    # As well as complete query contigs.
    df_complete_qry = pl.read_csv(
        bed_complete_qry,
        separator="\t",
        has_header=False,
        new_columns=[
            "chrom",
            "st",
            "end",
            "qry",
            "score",
            "stand",
            "qry_st",
            "qry_end",
            "item_rgb",
        ],
    ).select("qry", "qry_st", "qry_end")

    df_complete_ref = pl.read_csv(
        bed_complete_ref,
        separator="\t",
        has_header=False,
        new_columns=[
            "chrom",
            "st",
            "end",
            "ref",
            "score",
            "stand",
            "ref_st",
            "ref_end",
            "item_rgb",
        ],
    ).select("ref", "ref_st", "ref_end")

    df_final_liftover = (
        df_liftover.join(
            df_fai,
            on="qry",
            how="left",
        )
        # qry_right indicates that R1 is complete
        .join(df_complete_qry, on="qry", how="full")
        .select(
            "ref",
            "qry",
            "qry_length",
            "qry_right",
            "qry_st",
            "qry_end",
        )
        .join(df_complete_ref, on="ref", how="full")
        .with_columns(
            mtch_ref=pl.col("ref_right").str.extract_groups(
                r"^.*?_.*?(?<ref_chrom_name>chr[0-9XY]+)_.*?#(?<ref_hap>[12])#"
            ),
            mtch_qry=pl.col("qry_right").str.extract_groups(
                r"^.*?_.*?(?<qry_chrom_name>chr[0-9XY]+)_.*?#(?<qry_hap>[12])#"
            ),
        )
        .unnest("mtch_ref", "mtch_qry")
    )

    with open(outfile, "wt") as fh:
        for grp, df_grp in df_final_liftover.group_by(["ref"]):
            grp = grp[0]
            n_rows = df_grp.shape[0]
            # Multiple unique contigs
            if grp and n_rows > 1:
                # Multiple map but complete bed
                if df_grp["qry_right"].null_count() != n_rows:
                    status = "complete"
                else:
                    status = "cdr_break"
            # No mapping against R2 contig.
            # Some issues with p-arm fragments acrocentric chromosomes.
            # Try to recover by hap and sample
            elif grp is None:
                recovered_ovl = set()
                for row in df_grp.iter_rows(named=True):
                    ref = row["ref_right"]
                    qry = row["qry_right"]
                    status = "no_overlap"

                    if ref in recovered_ovl or qry in recovered_ovl:
                        continue

                    if ref:
                        df_qry = df_grp.filter(
                            pl.col("qry_chrom_name").eq(row["ref_chrom_name"])
                            & pl.col("qry_hap").eq(row["ref_hap"])
                        )
                        if not df_qry.is_empty():
                            qry = df_qry[0]["qry_right"][0]
                            status = "complete"
                            recovered_ovl.add(qry)
                    else:
                        df_ref = df_grp.filter(
                            pl.col("ref_chrom_name").eq(row["qry_chrom_name"])
                            & pl.col("ref_hap").eq(row["qry_hap"])
                        )
                        if not df_ref.is_empty():
                            ref = df_ref[0]["ref_right"][0]
                            status = "complete"
                        recovered_ovl.add(ref)
                    print(ref, qry, status, sep="\t", file=fh)
                continue
            # One contig
            else:
                row = df_grp.row(0, named=True)
                # complete bed
                if row["qry_right"] is not None:
                    status = "complete"
                else:
                    # break
                    # Double check and see if break at either edge.
                    if row["qry_st"] < 150_000 or row["qry_end"] > max(
                        0, row["qry_length"] - 150_000
                    ):
                        status = "cdr_break"
                    else:
                        status = "other_break"

            print(
                grp,
                ",".join(qry for qry in df_grp["qry"] if qry),
                status,
                sep="\t",
                file=fh,
            )


if __name__ == "__main__":
    raise SystemExit(main())
