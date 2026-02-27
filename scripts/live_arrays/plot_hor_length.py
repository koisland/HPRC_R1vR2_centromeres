import os
import argparse
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from collections import OrderedDict
from typing import Any, Callable


DEF_COLS = ("chrom", "chrom_st", "chrom_end", "length")
DEF_CHR_ORDER = [f"chr{i}" for i in (*range(1, 23), "X", "Y")]
DEF_CHR_COLORS = dict(
    zip(
        DEF_CHR_ORDER,
        (
            "#403E80",
            "#2C477D",
            "#346B9B",
            "#3C80AA",
            "#4587A2",
            "#589D96",
            "#73ACA4",
            "#87BAAF",
            "#94BE9F",
            "#9FC38C",
            "#A1C27C",
            "#A4C165",
            "#C2C969",
            "#B5A957",
            "#DFD06C",
            "#F2D46C",
            "#E3C765",
            "#E5BA61",
            "#D89E56",
            "#C8824A",
            "#BB6B3E",
            "#AB5C40",
            "#9B4C41",
            "#8C3C42",
        ),
    )
)
LEGEND_KWARGS = dict(
    loc="upper right",
    alignment="left",
    handlelength=1.0,
    handleheight=1.0,
    frameon=False,
)


def filter_R1(df: pl.DataFrame) -> pl.DataFrame:
    return (
        df
        # Only if we're showing total length, otherwise we need to show small albeit incomplete arrays.
        .filter(~(pl.col("chrom_name").str.contains("chr3"))).filter(
            pl.col("sample") != "HG03492"
        )
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
    return df


def add_mean_lines(
    ax: Axes,
    df: pl.DataFrame,
    by_col: str,
    val_col: str,
    fn_label: Callable[[int], str],
    colors: dict[str, str],
) -> list[tuple[str, str, int]]:
    means = []
    for grp, df_grp in df.group_by([by_col]):
        grp = grp[0]
        mean = df_grp[val_col].mean()
        label = fn_label(mean)
        color = colors[grp]
        # Draw line
        ax.axhline(mean, linestyle="dotted", color=color)
        # Set new labels
        yticks = ax.get_yticks()
        yticklabels = ax.get_yticklabels()
        ax.set_yticks([*yticks, mean], [*yticklabels, label])
        yticklabels = ax.get_yticklabels()
        # Last one is color
        yticklabels[-1].set_color(color)
        means.append((label, color, mean))
    return means

def read_lengths(
    file: str,
    rgx_name_groups: str,
    fn_filter: Callable[[pl.DataFrame], pl.DataFrame] | None = None,
) -> pl.DataFrame:
    df_lengths = (
        pl.read_csv(
            file,
            has_header=False,
            separator="\t",
            columns=[0, 1, 2, 3],
            new_columns=DEF_COLS,
        )
        .with_columns(source=pl.lit("samples"))
        .with_columns(
            mtch_chrom=pl.col("chrom").str.extract_groups(rgx_name_groups),
        )
        .unnest("mtch_chrom")
        .with_columns(pl.col("chrom_name").str.extract("(chr[0-9XY]+)"))
    )
    if fn_filter:
        df_lengths = fn_filter(df_lengths)

    return df_lengths


def perc_chrom(x: int, total: int) -> str:
    return f"{(x / total) * 100:.0f}"


def draw_only_bar(
    df_all_lengths: pl.DataFrame,
    palette: dict[str, Any],
    palette_order: dict[str, Any],
    outfile: str,
):
    samples = df_all_lengths["sample"].unique()
    nchroms_max = len(samples) * 2
    
    fig, ax = plt.subplots(
        figsize=(30, 10),
        layout="constrained",
    )
    df_all_length_counts = df_all_lengths.group_by(
        ["chrom_name", "source", "release"]
    ).agg(count=pl.col("chrom").count())
    sns.barplot(
        x="chrom_name",
        y="count",
        hue="release",
        data=df_all_length_counts,
        order=palette_order.keys(),
        hue_order=palette.keys(),
        palette=palette,
        ax=ax,
        legend="full",
        edgecolor="black",
    )
    for cont in ax.containers:
        ax.bar_label(cont)

    ax.set_ylim(0, nchroms_max)
    ax.set_ylabel("Number of α-satellite HOR arrays")

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    # x-axis
    ax.set_xlabel("Chromosome")
    # Remove chr from x-ticks
    xtick_labels = [lbl.get_text().replace("chr", "") for lbl in ax.get_xticklabels()]
    ax.set_xticks(ax.get_xticks(), xtick_labels)
    means = add_mean_lines(
        ax=ax,
        df=df_all_length_counts,
        by_col="release",
        val_col="count",
        fn_label=lambda x: str(int(x)),
        colors=palette,
    )
    sns.move_legend(ax, title=None, **LEGEND_KWARGS)

    labels, colors, means = zip(*means)

    yticks = ax.get_yticks()
    ax_2 = ax.secondary_yaxis(location="right")
    ax_2.set_yticks(
        [*yticks, *means],
        [*[perc_chrom(ytick, nchroms_max) for ytick in yticks], *[perc_chrom(mean, nchroms_max) for mean in means]],
    )
    yticklabels = ax_2.get_yticklabels()
    for i, color in enumerate(reversed(colors), 1):
        yticklabels[-i].set_color(color)

    ax_2.set_ylabel(r"% of α-satellite HOR arrays")

    fig.savefig(outfile, dpi=600, bbox_inches="tight")


def main():
    ap = argparse.ArgumentParser(
        description="Plot cumulative centromere HOR array lengths."
    )
    ap.add_argument(
        "-a",
        "--infile_a",
        help="Input centromere HOR array lengths.",
        type=str,
    )
    ap.add_argument(
        "-b",
        "--infile_b",
        help="Second input centromere HOR array lengths.",
        type=str,
    )
    ap.add_argument(
        "-m",
        "--mode",
        type=str,
        choices=["total", "arr"],
        default="arr",
        help="Plotting mode. Either total live array length (total) or live array length (arr).",
    )
    ap.add_argument(
        "--color_a",
        type=str,
        help="Color for a input.",
        default="blue",
    )
    ap.add_argument(
        "--color_b",
        type=str,
        help="Color for b input.",
        default="red",
    )
    ap.add_argument(
        "-c",
        "--chroms",
        nargs="+",
        default=DEF_CHR_ORDER,
        help="Chromosome names.",
    )
    ap.add_argument("-o", "--output", help="Output plot file.", required=True, type=str)
    args = ap.parse_args()

    palette = {
        "Release 1": args.color_a,
        "Release 2": args.color_b,
    }

    # Reverse to prevent matching chr1 with both chr1 and chr11
    chroms = list(args.chroms)
    chroms.extend([f"rc-{chrom}" for chrom in chroms])
    rgx_chrom = "|".join([*chroms, "-"])
    rgx_name_groups = (
        r"^(?<sample>.*?)_(?<chrom_name>(" + rgx_chrom + r")*)_(?<contig_name>.*?)$"
    )
    df_lengths = pl.concat(
        [
            read_lengths(
                args.infile_a, rgx_name_groups, fn_filter=filter_R1
            ).with_columns(release=pl.lit("Release 1")),
            read_lengths(
                args.infile_b, rgx_name_groups, fn_filter=filter_R2
            ).with_columns(release=pl.lit("Release 2")),
        ]
    )

    # Merge asat HOR array lengths
    if args.mode == "total":
        df_all_lengths = df_lengths.group_by(["chrom", "source"]).agg(
            pl.col("chrom_name").first(),
            pl.col("chrom_st").min(),
            pl.col("chrom_end").max(),
            pl.col("length").sum(),
            pl.col("release").first(),
            pl.col("sample").first(),
        )
        ylabel = "Cumulative length of α-satellite HOR arrays (Mbp)"
    else:
        df_all_lengths = df_lengths
        ylabel = "Length of α-satellite HOR arrays (Mbp)"

    df_all_lengths = df_all_lengths.with_columns(
        color_key=pl.when(pl.col("source") != "samples")
        .then(pl.col("source"))
        .otherwise(pl.col("chrom_name"))
    )
    # Add remaining chroms to plot if multi-chroms.
    palette_order = {p: i for i, p in enumerate(DEF_CHR_ORDER)}
    df_all_lengths_pd = df_all_lengths.to_pandas()

    fig, axes = plt.subplots(
        nrows=2,
        ncols=1,
        sharex=True,
        sharey=False,
        figsize=(30, 10),
        height_ratios=(0.3, 0.7),
        layout="constrained",
    )
    fig: plt.Figure
    ax_bar: Axes = axes[0]
    ax_violin: Axes = axes[1]

    df_all_length_counts = df_all_lengths.group_by(
        ["chrom_name", "source", "release"]
    ).agg(count=pl.col("chrom").count())
    sns.barplot(
        x="chrom_name",
        y="count",
        hue="release",
        data=df_all_length_counts,
        order=palette_order.keys(),
        hue_order=palette.keys(),
        palette=palette,
        ax=ax_bar,
        legend="full",
        edgecolor="black",
    )
    for cont in ax_bar.containers:
        ax_bar.bar_label(cont)

    ax_bar.set_ylabel("Number of α-satellite HOR arrays")

    sns.violinplot(
        x="chrom_name",
        y="length",
        hue="release",
        data=df_all_lengths_pd,
        order=palette_order.keys(),
        hue_order=palette.keys(),
        palette=palette,
        density_norm="width",
        inner="quart",
        ax=ax_violin,
        legend=None,
    )
    sns.stripplot(
        x="chrom_name",
        y="length",
        hue="release",
        data=df_all_lengths_pd,
        linewidth=0.5,
        dodge=True,
        edgecolor="black",
        order=palette_order.keys(),
        hue_order=palette.keys(),
        palette=palette,
        size=4,
        ax=ax_violin,
        legend=None,
    )

    sns.move_legend(ax_bar, title=None, **LEGEND_KWARGS)

    # Hide spines
    for ax in axes:
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)
    # x-axis
    ax_violin.set_xlabel("Chromosome")
    # Remove chr from x-ticks
    xtick_labels = [
        lbl.get_text().replace("chr", "") for lbl in ax_violin.get_xticklabels()
    ]
    ax_violin.set_xticks(ax_violin.get_xticks(), xtick_labels)
    # Set units of y-axis
    ax_violin.yaxis.minorticks_on()
    # Remove sci notation
    ax_violin.yaxis.set_major_formatter(lambda v, _: str(round(v / 1_000_000, 3)))
    ax_violin.set_ylabel(ylabel)

    # Add line and ytick for mean length and count
    add_mean_lines(
        ax=ax_violin,
        df=df_all_lengths,
        by_col="release",
        val_col="length",
        fn_label=lambda x: str(round(x / 1_000_000, 1)),
        colors=palette,
    )
    add_mean_lines(
        ax=ax_bar,
        df=df_all_length_counts,
        by_col="release",
        val_col="count",
        fn_label=lambda x: str(int(x)),
        colors=palette,
    )
    output_fname, _ = os.path.splitext(args.output)
    draw_only_bar(
        df_all_lengths, palette, palette_order, f"{output_fname}_bar_only.png"
    )
    fig.savefig(args.output, dpi=600, bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
