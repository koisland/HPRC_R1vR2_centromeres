import argparse

import polars as pl
import seaborn as sns

from typing import Callable
from matplotlib.text import Text


DEF_IN_COLS = ["ctg", "start", "end", "length"]
DEF_OUT_SCHEMA = {"sample": pl.String, "len": pl.UInt32, "perc": pl.Float64}
DEF_CHRS = list(reversed([f"chr{i}" for i in [*range(1, 23), "X", "Y"]]))
DEF_N_CHR = (len(DEF_CHRS) * 2) - 2


def get_df_counts(
    infile: str,
    rgx_name_groups: str,
    n_chroms: int,
    fn_filter: Callable[[pl.DataFrame], pl.DataFrame] | None = None,
) -> pl.DataFrame:
    df = (
        pl.read_csv(
            infile,
            separator="\t",
            has_header=False,
            new_columns=DEF_IN_COLS,
        )
        .group_by("ctg")
        .agg(pl.sum("length").alias("length"))
        .with_columns(mtch_ctg=pl.col("ctg").str.extract_groups(rgx_name_groups))
        .unnest("mtch_ctg")
    )
    if fn_filter:
        df = fn_filter(df)

    df_cnts = (
        df.group_by("sample")
        .len()
        .with_columns(perc=((pl.col("len") / n_chroms) * 100).round(1))
        .sort("sample")
    )
    return df_cnts


def filter_R1(df: pl.DataFrame) -> pl.DataFrame:
    return df.filter(~(pl.col("chrom_name").str.contains("chr3"))).filter(
        pl.col("sample") != "HG03492"
    )


def filter_R2(df: pl.DataFrame) -> pl.DataFrame:
    df = df.filter(
        ~(
            pl.col("contig_name").is_in(
                [
                    "HG00621#2#JAHBCC020000017.1:2-689335",
                    "HG00438#2#JAHBCA020000020.1:2-703696",
                ]
            )
        )
    )
    # breakpoint()
    return df


def perc_cens(n: int) -> str:
    return f"{(n / DEF_N_CHR) * 100:,.0f}"


def main():
    ap = argparse.ArgumentParser(
        description="Count complete centromeres from HOR array length."
    )
    ap.add_argument(
        "-a",
        "--input_a",
        required=True,
        type=str,
        help=f"Input HOR array length by contig with sample, chromosome, and contig/haplotype in name (ex. HG00171_chr1_h2tg000057l#1-9773347:114672-7639070). Expects TSV with no header and the fields: {DEF_IN_COLS}",
    )
    ap.add_argument(
        "-b",
        "--input_b",
        required=True,
        type=str,
        help=f"Input HOR array length by contig with sample, chromosome, and contig/haplotype in name (ex. HG00171_chr1_h2tg000057l#1-9773347:114672-7639070). Expects TSV with no header and the fields: {DEF_IN_COLS}",
    )
    ap.add_argument(
        "-o",
        "--output",
        default="complete_centromeres.png",
        type=str,
        help=f"Output centromere counts plot.",
    )
    ap.add_argument(
        "-c",
        "--chroms",
        default=DEF_CHRS,
        type=str,
        nargs="+",
        help="Chromosome names.",
    )
    ap.add_argument(
        "-n",
        "--n_chroms",
        default=DEF_N_CHR,
        help="Number of chromosomes in diploid organism.",
    )

    args = ap.parse_args()
    chroms = args.chroms
    n_chroms = args.n_chroms
    no_chrom = "all" in args.chroms

    # Include rc- in pattern
    if not no_chrom:
        chroms.extend([f"rc-{chrom}" for chrom in chroms])
        rgx_chrom = "|".join([*chroms, "-"])
        rgx_name_groups = (
            r"^(?<sample>.*?)_(?<chrom_name>(" + rgx_chrom + r")*)_(?<contig_name>.*?)$"
        )
    else:
        rgx_name_groups = r"^(?<sample>.*?)_(?<contig_name>.*?)$"

    # Remove all chr3 from R1 due to only one HOR array when 2 expected.
    # No HG03492
    df_a = get_df_counts(
        args.input_a,
        rgx_name_groups,
        n_chroms,
        fn_filter=filter_R1,
    )
    # Scaffold fragments.
    # HG00621_rc-chr2_HG00621#2#JAHBCC020000017.1:2-689335
    # HG00438_chrX_HG00438#2#JAHBCA020000020.1:2-703696
    df_b = get_df_counts(args.input_b, rgx_name_groups, n_chroms, fn_filter=filter_R2)
    df_all = pl.concat(
        [
            df_a.with_columns(release=pl.lit("Release 1")),
            df_b.with_columns(release=pl.lit("Release 2")),
        ]
    )
    print(df_all.group_by(["release"]).agg(pl.col("len").sum()))

    g = sns.catplot(
        data=df_all,
        kind="bar",
        x="sample",
        y="len",
        hue="release",
        height=8.0,
        aspect=2.0,
    )
    g.set_axis_labels("", "# of centromeres completely assembled")
    g.tick_params(axis="x", rotation=45)

    yticks = [*range(0, 50, 10), DEF_N_CHR]
    yticklabels = [str(v) for v in yticks]
    g.ax.set_yticks(yticks, yticklabels)
    g.ax.set_ylim(0, 46)

    mean_a = df_a["len"].mean()
    mean_b = df_b["len"].mean()
    g.ax.axhline(mean_a, color="blue", linestyle="dotted")
    g.ax.axhline(mean_b, color="orange", linestyle="dotted")

    ax_2 = g.ax.secondary_yaxis(location="right")
    ax_2.set_yticks(
        [*yticks, mean_a, mean_b],
        [*[perc_cens(ytick) for ytick in yticks], perc_cens(mean_a), perc_cens(mean_b)],
    )
    yticklabels = ax_2.get_yticklabels()
    for i, color in [(-2, "blue"), (-1, "orange")]:
        yticklabels[i].set_color(color)

    ax_2.set_ylabel(r"% of centromeres completely assembled")

    sns.move_legend(
        g,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.075),
        ncol=2,
        title=None,
        frameon=False,
    )
    g.savefig(args.output, bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
