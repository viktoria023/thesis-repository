import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from statsmodels.stats.proportion import proportions_ztest
from statsmodels.stats.proportion import proportion_confint

# make the default font size point 7
plt.rcParams.update({"font.size": 7})


def plot_compensation_barplot(
    major_percentage,
    minor_percentage,
    save_name="contingency_barplot",
    save_figure=False,
):

    # plot the two percentages as a nice bar plot
    fig = plt.figure(figsize=(2.75, 2.5))
    axis = fig.gca()

    axis.bar(
        ["Homogeneous\nsamples", "Heterogeneous\nsamples"],
        [major_percentage, minor_percentage],
        color=["#377eb8", "#e41a1c"],
        edgecolor=None,
    )
    axis.set_ylabel("% with compensatory mutation")
    axis.set_ylim(0, 0.5)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)

    # display percentages inside bars
    for i in range(2):
        axis.text(
            i,
            [major_percentage, minor_percentage][i] + 0.01,
            f"{[major_percentage, minor_percentage][i]:.1%}",
            ha="center",
            va="bottom",
            color="black",
        )

    # Set axes to display values as percentages without the percent sign
    axis.yaxis.set_major_formatter(FuncFormatter(lambda y, _: "{:.0f}".format(y * 100)))

    if save_figure:
        # save figure in high resolution
        fig.savefig(f"pdf/{save_name}.svg", format='svg', dpi=1200, bbox_inches='tight')

    plt.show()


def plot_stats_histogram(
    stats,
    significant=True,
    save_figure=False,
    save_name="sensitivity_specificity_histogram",
    colours=[
        "#ff7f0e",
        "#1f77b4",
    ],
    metric_labels=["senstivity", "specificity"],
    legend_label_1="FRS cutoff >=0.9",
    legen_label_2="No FRS cutoff",
):
    # Adjusting the plot to group bars by criteria and display percentage inside the bar
    fig, ax = plt.subplots(figsize=(2.5, 2.75))

    # Plotting bars next to each other
    x = [0, 0.7]  # positions for FRS cutoff >=0.9 (major allele definition)
    bar_width = 0.25

    bars1 = ax.bar(
        x,
        stats["FRS cutoff >=0.9 (major allele definition)"],
        width=bar_width,
        color=colours,
        edgecolor=None,
        linewidth=1.2,
        label=legend_label_1,
    )
    bars2 = ax.bar(
        [p + bar_width for p in x],
        stats["No FRS cutoff"],
        width=bar_width,
        color=["white", "white"],
        edgecolor=colours,
        linewidth=1.2,
        label=legen_label_2,
    )

    # bars1[0].set_color('#ff7f0e')

    # Adding percentage values inside the bars
    for bars in [bars1, bars2]:
        for bar in bars:
            yval = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                yval + 0.02,
                f"{yval:.1%}",
                ha="center",
                va="bottom",
                fontsize=7,
                color="k",
            )

    # Customize the plot
    # Set y-axis to display values as percentages without the percent sign
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: "{:.0f}".format(y * 100)))

    # ax.set_title('Sensitivity and Specificity for Rifampicin Resistance', fontsize=18)
    # ax.set_ylabel('%', fontsize=7)
    ax.set_ylim(0, 1.1)

    # Setting custom x-axis labels
    ax.set_xticks([r + bar_width / 2 for r in x])
    ax.set_xticklabels(metric_labels, fontsize=7)

    # Increase axis fontsizes
    ax.tick_params(axis="both", which="major", labelsize=7)

    # Adding the legend
    # ax.legend(loc='lower right', fontsize=12)

    # remove unnecessary spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # add significance indicator with horizontal bar between bars of sensitivity
    if significant:
        ax.plot([0, 0], [1.05, 1.07], color="black", lw=2)
        ax.plot([0.2, 0.2], [1.05, 1.07], color="black", lw=2)
        ax.plot([0, 0.2], [1.07, 1.07], color="black", lw=2)
        ax.text(0.1, 1.08, "*", fontsize=20, color="black", ha="center", va="center")

    plt.tight_layout()
    plt.show()

    # save figure if print_stats is True
    if save_figure:
        fig.savefig(f"pdf/{save_name}.svg", format='svg', dpi=1200, bbox_inches='tight')


def plot_performance_frs(
    df,
    save_figure=False,
    save_name="fig-sensitivity_specificity_FRS_cutoff",
    metrics=["sensitivity", "specificity"],
    colours=["#ff7f0e", "#1f77b4"],
):

    fig, ax1 = plt.subplots(figsize=(5, 3))

    # plot sensitivity and specificity using one y axis
    for i in [0, 1]:
        ax1.errorbar(
            df.min_FRS,
            df[metrics[i]],
            yerr=df[metrics[i] + "_error"],
            label=metrics[i],
            color=colours[i],
            elinewidth=1.0,
            linewidth=0.5,
            marker="o",
            markersize=3,
        )
    ax1.set_xlabel("Fraction of reads (FRS) supporting RAV")

    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.set_ylim(0.94, 0.99)
    ax1.set_yticks(np.arange(0.94, 0.99, 0.01))
    ax1.set_xlim(0, 1)

    # Set axes to display values as percentages without the percent sign
    plt.gca().yaxis.set_major_formatter(
        FuncFormatter(lambda y, _: "{:.0f}".format(y * 100))
    )
    plt.gca().xaxis.set_major_formatter(
        FuncFormatter(lambda x, _: "{:.0f}".format(x * 100))
    )

    if save_figure:
        fig.savefig(
            f"pdf/{save_name}.svg",
            format='svg', dpi=1200, bbox_inches='tight'
        )

    plt.show()


def plot_histogram(
    df,
    save_figure=False,
    save_name="sensitivity_specificity_histogram",
    colours=[
        "#ff7f0e",
        "#1f77b4",
    ],
    metrics=None,
):
    # Adjusting the plot to group bars by criteria and display percentage inside the bar
    fig, ax = plt.subplots(figsize=(2.5, 2.75))

    # Plotting bars next to each other
    x = [0, 0.7]  # positions for FRS cutoff >=0.9 (major allele definition)
    bar_width = 0.25
    b1 = ax.bar(0, df.iloc[0][metrics[0]], width=bar_width, color=colours[0])
    b2 = ax.bar(0.7, df.iloc[0][metrics[1]], width=bar_width, color=colours[1])
    b3 = ax.bar(
        0.25,
        df.iloc[1][metrics[0]],
        width=bar_width,
        color="white",
        edgecolor=colours[0],
    )
    b4 = ax.bar(
        0.95,
        df.iloc[1][metrics[1]],
        width=bar_width,
        color="white",
        edgecolor=colours[1],
    )

    # Adding percentage values inside the bars
    for bar in [b1, b2, b3, b4]:
        yval = bar[0].get_height()
        ax.text(
            bar[0].get_x() + bar[0].get_width() / 2,
            yval + 0.02,
            f"{yval:.1%}",
            ha="center",
            va="bottom",
            fontsize=7,
            color="k",
        )

    # Customize the plot
    # Set y-axis to display values as percentages without the percent sign
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: "{:.0f}".format(y * 100)))

    ax.set_ylim(0, 1.1)

    # Setting custom x-axis labels
    ax.set_xticks([r + bar_width / 2 for r in x])
    ax.set_xticklabels(metrics, fontsize=7)

    ax.tick_params(axis="both", which="major", labelsize=7)

    # remove unnecessary spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # perform statistical test for significant difference
    a = df.iloc[0]
    b = df.iloc[1]

    for i in range(len(metrics)):
        metric = metrics[i]
        if metric == "sensitivity":
            counts = [int(a.TP), int(b.TP)]
            nobs = [int(a.P), int(b.P)]
        elif metric == "specificity":
            counts = [int(a.TN), int(b.TN)]
            nobs = [int(a.N), int(b.N)]
        elif metric == "PPV":
            counts = [int(a.TP), int(b.TP)]
            nobs = [int(a.TP) + int(a.FP), int(b.TP) + int(b.FP)]
        elif metric == "NPV":
            counts = [int(a.TN), int(b.TN)]
            nobs = [int(a.TN) + int(a.FN), int(b.TN) + int(b.FN)]

        z, p = proportions_ztest(counts, nobs)
        if p < 0.05:
            print(metric, z, p)
            x_value = x[i]
            ax.plot(
                [x_value + 0.025, x_value + 0.025], [1.08, 1.1], color="black", lw=2
            )
            ax.plot(
                [x_value + 0.225, x_value + 0.225], [1.08, 1.1], color="black", lw=2
            )
            ax.plot([x_value + 0.025, x_value + 0.225], [1.1, 1.1], color="black", lw=2)
            ax.text(
                x_value + 0.125,
                1.125,
                "*",
                fontsize=20,
                color="black",
                ha="center",
                va="center",
            )

    plt.tight_layout()
    plt.show()

    # save figure if print_stats is True
    if save_figure:
        fig.savefig(f"pdf/{save_name}.svg", format='svg', dpi=1200, bbox_inches='tight')


def plot_frs_histogram(df, savefig=False, save_name="hist-FRS-RAV"):

    axes = df.MAX_RAV_FRS.hist(bins=np.arange(0.05, 1, 0.05), color="grey")

    axes.set_ylabel("Number of samples")
    axes.set_xlabel("Fraction of reads supporting RAV")
    axes.spines["top"].set_visible(False)
    axes.spines["right"].set_visible(False)
    axes.set_xlim(0, 1)
    axes.grid(False)
    fig = axes.get_figure()
    fig.set_size_inches(3, 2)
    if savefig:
        fig.savefig(f"pdf/{save_name}.svg", format='svg', dpi=1200, bbox_inches='tight')

def plot_het_frs_histogram(df, RAV_FRS=None, lineage_count=None, savefig=False, i=None):
    """
    Plot the FRS distribution for a given DataFrame with optional RAV FRS line.

    Parameters:
    df (pd.DataFrame): DataFrame containing 'FRS' column.
    RAV_FRS (float): RAV FRS value to mark with vertical line.
    lineage_count (int): Number of lineages for the RAV variant.
    savefig (bool): Whether to save the figure as a PDF.
    """
    # Create figure with better styling
    fig, ax = plt.subplots(figsize=(3.5, 1.8))

    # Create histogram with improved styling
    n, bins, patches = ax.hist(
        df["FRS"], bins=np.arange(0.05, 1, 0.025), color="red", alpha=0.5, linewidth=0.5
    )

    # Add vertical line(s) for RAV FRS if provided
    if RAV_FRS is not None:
        label_text = "RAV FRS"
        ax.axvline(
            x=RAV_FRS,
            color="black",
            linestyle="--",
            linewidth=1.5,
            alpha=0.9,
            label=label_text,
        )

        # Add legend with better styling
        # legend = ax.legend(loc='upper right', framealpha=0.9)
        # legend.get_frame().set_facecolor('#f8f9fa')

    # Improve axis labels and title
    ax.set_xlabel("FRS")
    # ax.set_ylabel('Number of Mutations', fontsize=12)

    # Remove top and right spines for cleaner look
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    # ax.spines['left'].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.set_yticks([])  # Remove y-ticks for cleaner look

    # Improve tick styling
    ax.tick_params(axis="both", which="major", labelsize=7, length=4, width=0.8)
    ax.tick_params(axis="both", which="minor", length=2, width=0.5)

    # Set x-axis limits and add grid
    ax.set_xlim(0, 1)
    # ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)

    # Improve layout
    plt.tight_layout()

    if savefig:
        if RAV_FRS is None:
            plt.savefig(
                f"pdf/het-samples/not-annotated/hist-frs-distribution-{i}.pdf",
                bbox_inches="tight",
                transparent=True,
                dpi=300,
            )
        else:
            plt.savefig(
                f"pdf/het-samples/annotated/hist-frs-distribution-{i}.pdf",
                bbox_inches="tight",
                transparent=True,
                dpi=300,
            )

    if i < 5 and RAV_FRS is not None:
        plt.show()
    else:
        plt.close()


def calculate_results_row(df, description=None, min_FRS=None):

    true_positives = int(
        df[(df["IS_RESISTANT"]) & (df["HAS_RESISTANT_MUTATION"])].shape[0]
    )
    positives = int(df[df["IS_RESISTANT"]].shape[0])
    false_positives = int(
        df[(~df["IS_RESISTANT"]) & (df["HAS_RESISTANT_MUTATION"])].shape[0]
    )

    true_negatives = int(
        df[(~df["IS_RESISTANT"]) & (~df["HAS_RESISTANT_MUTATION"])].shape[0]
    )
    negatives = int(df[~df["IS_RESISTANT"]].shape[0])
    false_negatives = int(
        df[(df["IS_RESISTANT"]) & (~df["HAS_RESISTANT_MUTATION"])].shape[0]
    )

    sensitivity = true_positives / positives

    # only need to calculate one direction as the numbers are sufficient large the limits are symmetric
    sensitivity_error = [
        abs(i - sensitivity)
        for i in proportion_confint(
            true_positives, positives, alpha=0.05, method="normal"
        )
    ][0]

    specificity = true_negatives / negatives

    specificity_error = [
        abs(i - specificity)
        for i in proportion_confint(
            true_negatives, negatives, alpha=0.05, method="normal"
        )
    ][0]

    ppv = true_positives / (true_positives + false_positives)

    ppv_error = [
        abs(i - ppv)
        for i in proportion_confint(
            true_positives,
            true_positives + false_positives,
            alpha=0.05,
            method="normal",
        )
    ][0]

    npv = true_negatives / (true_negatives + false_negatives)

    npv_error = [
        abs(i - npv)
        for i in proportion_confint(
            true_negatives,
            true_negatives + false_negatives,
            alpha=0.05,
            method="normal",
        )
    ][0]

    return pd.Series(
        [
            description,
            min_FRS,
            true_positives,
            false_positives,
            positives,
            true_negatives,
            false_negatives,
            negatives,
            sensitivity,
            specificity,
            ppv,
            npv,
            specificity_error,
            sensitivity_error,
            ppv_error,
            npv_error,
        ],
        index=[
            "description",
            "min_FRS",
            "TP",
            "FP",
            "P",
            "TN",
            "FN",
            "N",
            "sensitivity",
            "specificity",
            "PPV",
            "NPV",
            "sensitivity_error",
            "specificity_error",
            "PPV_error",
            "NPV_error",
        ],
    )
