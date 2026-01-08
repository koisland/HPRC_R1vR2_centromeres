import sys
import argparse
import polars as pl

EXP_COLS = (
    "chrom", "st", "end",
    "ochrom", "ost", "oend", "olen", "ohor", "operc",
    "nchrom",
)
OUT_COLS = (
    "ochrom", "ost", "oend", "olen", "ohor", "operc"
)
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    args = ap.parse_args()
    infile = args.infile
    df = pl.read_csv(
        infile,
        separator="\t",
        has_header=False,
        columns=range(len(EXP_COLS)),
        new_columns=EXP_COLS
    )
    for _, df_grp in df.group_by(["nchrom"]):
        nchroms = df_grp["nchrom"].unique()
        # Has an error in one or more arrays. Remove completely.
        if any(chrom != "." for chrom in nchroms):
            continue
        df_grp.select(OUT_COLS).write_csv(sys.stdout, separator="\t", include_header=False)

if __name__ == "__main__":
    raise SystemExit(main())
