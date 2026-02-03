import argparse
import numpy as np
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from collections import OrderedDict
from typing import Callable


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


def filter_R1(df: pl.DataFrame) -> pl.DataFrame:
    return (
        df
        # .filter(~(pl.col("chrom_name").str.contains("chr3")))
        .filter(pl.col("sample") != "HG03492")
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
        "--added_inputs",
        help="Additional centromere HOR array lengths. First column should be be a subset of --chroms",
        nargs="*",
        metavar="{lbl}={path}={color}",
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
        "-c",
        "--chroms",
        nargs="+",
        default=DEF_CHR_ORDER,
        help="Chromosome names.",
    )
    ap.add_argument(
        "--chrom_colors",
        default=None,
        help="Chromosome colors as TSV file with chrom to color mapping.",
    )
    ap.add_argument("-o", "--output", help="Output plot file.", required=True, type=str)
    args = ap.parse_args()

    # Reverse to prevent matching chr1 with both chr1 and chr11
    chroms = args.chroms
    chroms.extend([f"rc-{chrom}" for chrom in chroms])
    rgx_chrom = "|".join([*chroms, "-"])
    rgx_name_groups = (
        r"^(?<sample>.*?)_(?<chrom_name>(" + rgx_chrom + r")*)_(?<contig_name>.*?)$"
    )
    if args.chrom_colors:
        _chroms = set(args.chroms)
        with open(args.chrom_colors) as fh:
            chrom_colors = {}
            for line in fh:
                chrom, color = line.strip().split("\t")
                if chrom not in _chroms:
                    continue
                chrom_colors[chrom] = color
    else:
        chrom_colors = DEF_CHR_COLORS

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

    added_palettes = OrderedDict()
    dfs_added_lengths = []
    if args.added_inputs:
        for elems in (lbl_path.split("=") for lbl_path in args.added_inputs):
            # Allow color to be optional.
            try:
                lbl, path, color = elems
            except ValueError:
                lbl, path = elems
                # Generate random color.
                color = np.random.rand(3)

            df = pl.read_csv(
                path,
                has_header=False,
                separator="\t",
                columns=[0, 1, 2, 3],
                new_columns=DEF_COLS,
            ).with_columns(source=pl.lit(lbl))

            df = df.with_columns(
                chrom_name=pl.col("chrom").str.extract(f"^({rgx_chrom})$")
            )

            added_palettes[lbl] = color
            dfs_added_lengths.append(df)

    # Get order of chromosomes
    palettes = chrom_colors | added_palettes
    df_all_lengths: pl.DataFrame = pl.concat([df_lengths, *dfs_added_lengths])

    # Merge asat HOR array lengths
    if args.mode == "total":
        df_all_lengths = df_all_lengths.group_by(["chrom", "source"]).agg(
            pl.col("chrom_name").first(),
            pl.col("chrom_st").min(),
            pl.col("chrom_end").max(),
            pl.col("length").sum(),
            pl.col("release").first(),
        )
        ylabel = "Cumulative length of α-satellite HOR arrays (Mbp)"
    else:
        ylabel = "Length of α-satellite HOR arrays (Mbp)"

    df_all_lengths = df_all_lengths.with_columns(
        color_key=pl.when(pl.col("source") != "samples")
        .then(pl.col("source"))
        .otherwise(pl.col("chrom_name"))
    )
    # Add remaining chroms to plot if multi-chroms.
    uncovered_chroms = set(df_all_lengths["chrom_name"].unique()).difference(
        palettes.keys()
    )
    palettes = palettes | {chrom: "#FFFFFF" for chrom in uncovered_chroms}
    palette_order = {p: i for i, p in enumerate(palettes.keys())}
    df_all_lengths_pd = df_all_lengths.to_pandas()
    hue_order = ["Release 1", "Release 2"]
    sns.violinplot(
        x="chrom_name",
        y="length",
        hue="release",
        data=df_all_lengths_pd,
        order=palette_order.keys(),
        hue_order=hue_order,
        inner="quart",
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
        hue_order=hue_order,
        size=4,
    )

    ax = plt.gca()
    legend_kwargs = dict(
        loc="upper right",
        alignment="left",
        frameon=False,
    )
    try:
        # Sort legend elements
        handles_labels = ax.get_legend_handles_labels()
        # Only display dots.
        handles, labels = zip(
            *sorted(
                (
                    (handle, length)
                    for handle, length in zip(*handles_labels)
                    if isinstance(handle, Line2D)
                ),
                key=lambda x: palette_order.get(x[1], -1),
            )
        )
        # Place outside of figure.
        ax.legend(handles, labels, **legend_kwargs)
    except ValueError:
        ax.legend(
            handles=[
                Patch(color=color, label=chrom) for chrom, color in chrom_colors.items()
            ],
            **legend_kwargs,
        )

    # Hide spines
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    # x-axis
    ax.set_xlabel("Chromosome")
    # Remove chr from x-ticks
    xtick_labels = [lbl.get_text().replace("chr", "") for lbl in ax.get_xticklabels()]
    ax.set_xticks(ax.get_xticks(), xtick_labels)

    # Set units of y-axis
    ax.yaxis.minorticks_on()
    # Remove sci notation
    ax.yaxis.set_major_formatter("plain")
    new_xtick_labels = []
    _, ymax = ax.get_ylim()
    yticks, yticklabels = ax.get_yticks(), ax.get_yticklabels()
    for txt in yticklabels:
        _, y = txt.get_position()
        # Convert units and round.
        new_y_txt = str(round(y / 1_000_000, 3))
        txt.set_text(new_y_txt)
        new_xtick_labels.append(txt)

    # Add line and ytick for mean length
    mean_length_r1 = df_all_lengths.filter(pl.col("release") == "Release 1")[
        "length"
    ].mean()
    mean_length_r2 = df_all_lengths.filter(pl.col("release") == "Release 2")[
        "length"
    ].mean()
    mean_length_r1_label = str(round(mean_length_r1 / 1_000_000, 1))
    mean_length_r2_label = str(round(mean_length_r2 / 1_000_000, 1))
    ax.axhline(mean_length_r1, linestyle="dotted", color="blue")
    ax.axhline(mean_length_r2, linestyle="dotted", color="orange")
    ax.set_yticks(
        [*yticks, mean_length_r1, mean_length_r2],
        [*new_xtick_labels, mean_length_r1_label, mean_length_r2_label],
    )

    yticklabels = ax.get_yticklabels()
    for i, color in [(-2, "blue"), (-1, "orange")]:
        yticklabels[i].set_color(color)

    ax.set_ylabel(ylabel)
    ax.set_ylim(0, ymax)

    plt.gcf().set_size_inches(30, 8)
    plt.savefig(args.output, dpi=600, bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
