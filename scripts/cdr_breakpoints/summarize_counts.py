import os
import sys
import numpy as np
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt

from os.path import join
from typing import Literal as L
from matplotlib.axes import Axes
from matplotlib.patches import Patch
from scipy.stats import fisher_exact, false_discovery_control

# ORDER = [*[str(i) for i in range(1, 23)], "X", "Y"]
ORDER = ["All", *[str(i) for i in range(1, 23)], "X", "Y"]
HUE_ORDER = {
    "No break in centromeric region": "green",
    "Break in CDR": "red",
    "Break in non-CDR centromeric region": "orange",
}
HUE_ORDER_BREAKS = ("Break in CDR", "Break in non-CDR centromeric region")
ALPHA = 0.05
HEADER_BREAKS = ("chrom", "status", "cnt", "prop", "type", "pvalue")


# https://rowannicholls.github.io/python/graphs/ax_based/boxplots_significance.html
def draw_signif_brackets(
    ax: Axes, x1: int, x2: int, level: int, ylim: tuple[float, float], p: float
):
    # What level is this bar among the bars above the plot?
    # Plot the bar
    y_range = ylim[1] - ylim[0]
    bar_height = (y_range * 0.09 * level) + ylim[1]
    bar_tips = bar_height - (y_range * 0.02)
    ax.plot(
        [x1, x1, x2, x2],
        [bar_tips, bar_height, bar_height, bar_tips],
        lw=1,
        c="k",
        clip_on=False,
    )
    # Significance level
    if p < 0.001:
        sig_symbol = "***"
    elif p < 0.01:
        sig_symbol = "**"
    elif p < 0.05:
        sig_symbol = "*"
    text_height = bar_height + (y_range * 0.01)
    ax.text((x1 + x2) * 0.5, text_height, sig_symbol, ha="center", va="bottom", c="k")


def main():
    infile = sys.argv[1]
    output_dir = sys.argv[2]

    df = pl.read_csv(
        infile, separator="\t", has_header=False, new_columns=["r2", "r1", "status"]
    )
    df_filtered = (
        df
        # We need to filter out any that have None on either ref (missing in R2) or qry (Missing in R1)
        # no_overlap
        .filter(~pl.col("r2").eq(pl.lit("None")) & ~pl.col("r1").eq(pl.lit("None")))
        .with_columns(
            mtch=pl.col("r2").str.extract_groups(
                r"^(?<sm>.*?)_.*?chr(?<chrom>[0-9XY]+)_"
            )
        )
        .unnest("mtch")
        .cast({"chrom": pl.Enum(ORDER)})
    )
    df_chrom_zeroes = pl.DataFrame(
        [(chrom, status, 0) for chrom in ORDER for status in HUE_ORDER],
        orient="row",
        schema=["chrom", "status", "cnt"],
    ).cast({"chrom": pl.Enum(ORDER)})

    # Generate summary stats
    df_agg = (
        df_filtered.group_by(["chrom", "status"])
        .agg(cnt=pl.col("sm").count())
        .sort(by=["chrom", "status"])
    )
    df_agg_all = (
        df_filtered.group_by(["status"])
        .agg(chrom=pl.lit("All"), cnt=pl.col("sm").count())
        .cast({"chrom": pl.Enum(ORDER)})
        .select("chrom", "status", "cnt")
    )
    df_agg = (
        pl.concat([df_agg, df_agg_all])
        .join(df_chrom_zeroes, on=["chrom", "status"], how="full")
        .with_columns(
            chrom=pl.when(pl.col("chrom").is_null())
            .then(pl.col("chrom_right"))
            .otherwise(pl.col("chrom")),
            status=pl.when(pl.col("status").is_null())
            .then(pl.col("status_right"))
            .otherwise(pl.col("status")),
            cnt=pl.when(pl.col("cnt").is_null())
            .then(pl.col("cnt_right"))
            .otherwise(pl.col("cnt")),
        )
        .select("chrom", "status", "cnt")
        .with_columns(prop=(pl.col("cnt") / pl.col("cnt").sum().over("chrom")) * 100)
    )

    fig, ax = plt.subplots(layout="constrained", figsize=(20, 5))
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    sns.barplot(
        df_agg.cast({"chrom": pl.String}),
        x="chrom",
        y="prop",
        hue="status",
        hue_order=HUE_ORDER.keys(),
        palette=HUE_ORDER,
        order=ORDER,
        ax=ax,
    )
    sns.move_legend(
        ax,
        loc="upper right",
        ncol=1,
        handlelength=1.0,
        handleheight=1.0,
        title=None,
        frameon=False,
    )
    for status, cont in zip(HUE_ORDER, ax.containers, strict=True):
        def get_lbl(x):
            cnt = df_agg.filter(
                pl.col("status").eq(pl.lit(status))
                & pl.col("prop").eq(pl.lit(x))
            )["cnt"][0]
            return f"{x:.0f}%\n({cnt})"
        ax.bar_label(cont, fmt=get_lbl, fontsize="x-small")

    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Proportion of HPRC release 1 centromeres (%)")

    fig.savefig(
        join(output_dir, "cdr_breaks_summary.png"), dpi=600, bbox_inches="tight"
    )
    df_agg.write_csv(join(output_dir, "cdr_breaks_summary.tsv"), separator="\t")

    fig, axes = plt.subplots(
        ncols=6,
        nrows=4,
        figsize=(18, 12),
        sharex=True,
        sharey=True,
        layout="constrained",
    )
    axes: np.ndarray[tuple[L[4], L[6]], np.dtype[Axes]]
    # Then Fisher's exact test
    # * H_o: The proportion of breaks in the CDR is independent of chromosome
    # https://www.statology.org/fishers-exact-test-python/
    # https://xkcd.com/882/
    break_default_counts = {s: 0 for s in HUE_ORDER_BREAKS}
    break_colors = {s: HUE_ORDER[s] for s in HUE_ORDER_BREAKS}
    # {ax_idx: p_value}
    breaks_p_values = {}
    breaks_data = {}
    chroms = df_filtered["chrom"].unique().sort()
    for i, chrom in enumerate(chroms):
        row, col = divmod(i, 6)
        ax: Axes = axes[row, col]

        df_chrom = df_filtered.filter(
            pl.col("chrom").eq(chrom)
            & pl.col("status").ne(pl.lit("No break in centromeric region"))
        )
        df_chrom_other = df_filtered.filter(
            ~pl.col("chrom").eq(chrom)
            & pl.col("status").ne(pl.lit("No break in centromeric region"))
        )
        # Fill with 0 if not present.
        chrom_statuses = break_default_counts | dict(
            df_chrom["status"].value_counts(sort=True).iter_rows()
        )
        chrom_other_statuses = break_default_counts | dict(
            df_chrom_other["status"].value_counts(sort=True).iter_rows()
        )

        data = [
            list(chrom_statuses.values()),
            list(chrom_other_statuses.values()),
        ]
        res = fisher_exact(data, alternative="two-sided")
        breaks_p_values[(row, col)] = res.pvalue

        # Convert to proportion to plot.
        df_agg_status_prop = pl.DataFrame(
            [
                *((s, v, "Chromosome") for s, v in chrom_statuses.items()),
                *((s, v, "Other Chromosomes") for s, v in chrom_other_statuses.items()),
            ],
            orient="row",
            schema=["status", "cnt", "type"],
        ).with_columns(prop=(pl.col("cnt") / pl.col("cnt").sum().over("type")))
        breaks_data[(row, col)] = df_agg_status_prop
        ax.set_title(f"Chromosome {chrom}")
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

        sns.barplot(
            df_agg_status_prop,
            x="type",
            y="prop",
            hue="status",
            order=("Chromosome", "Other Chromosomes"),
            hue_order=HUE_ORDER_BREAKS,
            palette=break_colors,
            ax=ax,
            legend=None,
        )
        for cont in ax.containers:
            ax.bar_label(cont, fmt=lambda pct: f"{pct * 100:.1f}")

        ax.yaxis.set_major_formatter(lambda x, pos: round(x * 100, 1))
        ax.set_xlabel(None)
        ax.set_ylabel(None)

    # FDR adjust
    # Add signif bar.
    adj_pvalues = false_discovery_control(list(breaks_p_values.values()), method="bh")
    try:
        os.remove(join(output_dir, "fisher_exact_by_chrom.tsv"))
    except Exception:
        pass
    with open(join(output_dir, "fisher_exact_by_chrom.tsv"), "at") as fh:
        print(*HEADER_BREAKS, sep="\t", file=fh)
        for chrom, (idx, pvalue), (_, df_data), adj_pvalue in zip(
            chroms,
            breaks_p_values.items(),
            breaks_data.items(),
            adj_pvalues,
            strict=True,
        ):
            df_data = df_data.with_columns(
                chrom=pl.lit(chrom), pvalue=pl.lit(adj_pvalue)
            ).select(HEADER_BREAKS)

            df_data.write_csv(fh, separator="\t", include_header=False)
            if adj_pvalue > ALPHA:
                continue
            ax: Axes = axes[*idx]
            draw_signif_brackets(
                ax=ax, x1=0, x2=1, level=0, ylim=[0.0, 1.0], p=adj_pvalue
            )

    fig.legend(
        labels=HUE_ORDER_BREAKS,
        handles=[Patch(color=c) for s, c in break_colors.items()],
        ncol=2,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.01),
        handlelength=1.0,
        handleheight=1.0,
        title=None,
        frameon=False,
    )
    fig.supylabel("Proportion of breaks in HPRC release 1 centromeres (%)")
    fig.savefig(
        join(output_dir, "fisher_exact_by_chrom.png"), dpi=600, bbox_inches="tight"
    )


if __name__ == "__main__":
    raise SystemExit(main())
