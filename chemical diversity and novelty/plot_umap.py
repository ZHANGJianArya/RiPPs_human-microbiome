#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import umap
import pandas as pd


def plot_umap(df_umap, file, title='UMAP'):
    """
    Create and save a UMAP plot from a given DataFrame.

    Parameters:
    df_umap (DataFrame): DataFrame containing UMAP-transformed data.
    file (str): Path to save the UMAP plot.
    title (str, optional): Title of the plot. Defaults to 'UMAP'.

    Returns:
    None
    """
    plt.figure(figsize=(16, 16))

    labels = ['RiPP_from_HM', 'RiPP_from_MIBIG']
    colors = ['#000000', '#ff001d']

    for color, i, label in zip(colors, range(1, 10), labels):
        plt.scatter(df_umap.loc[df_umap["prediction"] == label, ["umap1"]],
                    df_umap.loc[df_umap["prediction"] == label, ["umap2"]],
                    color=color, lw=5, label=label, s=5)

    plt.title(title, fontsize=50)
    lgnd = plt.legend(fontsize=20)
    for i in range(2):
        # set legend scatter size
        lgnd.legendHandles[i]._sizes = [60]
    plt.axis()
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.savefig(file)

    svg = str(file) + '.svg'
    png = str(file) + '.png'
    plt.savefig(svg, format="svg", transparent=False)
    plt.savefig(png, format="png", transparent=False)


def draw_umap(df, file, metric='jaccard', title='UMAP'):
    """
    Compute UMAP and create a UMAP plot from a given DataFrame.

    Parameters:
    df (DataFrame): DataFrame containing original data with ECFP6 fingerprints.
    file (str): Path to save the UMAP plot.
    metric (str, optional): Metric used for UMAP. Defaults to 'jaccard'.
    title (str, optional): Title of the plot. Defaults to 'UMAP'.

    Returns:
    None
    """
    fit = umap.UMAP(metric=metric)

    ecfps = list(df["ECFP6"])
    fingerprints = [[*i] for i in ecfps]
    fingerprints = np.asarray(fingerprints, dtype=int)
    u = fit.fit_transform(fingerprints)
    umap_df = pd.DataFrame(u, columns=["umap1", "umap2"])
    df_umap = pd.concat([df, umap_df], axis=1)
    plt.figure(figsize=(24, 24))

    labels = ['RiPP_from_HM', 'RiPP_from_MIBIG']
    colors = ['#000000', '#ff001d']

    for color, i, label in zip(colors, range(1, 10), labels):
        plt.scatter(df_umap.loc[df_umap["prediction"] == label, ["umap1"]],
                    df_umap.loc[df_umap["prediction"] == label, ["umap2"]],
                    color=color, lw=2, label=label, s=1)

    plt.title(title, fontsize=50)
    lgnd = plt.legend(fontsize=20)
    for i in range(4):
        # set legend scatter size
        lgnd.legendHandles[i]._sizes = [30]
    plt.axis()
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.savefig(file)


if __name__ == '__main__':
    import argparse

    argparser = argparse.ArgumentParser(description="plot UMAP")
    argparser.add_argument("-i", "--input", required=True,
                           help="the input table file of peptide sequences with ECFP6")
    argparser.add_argument("-o", "--output", help="the output umap plot", required=True)
    argparser.add_argument("-m", "--metric", choices=['euclidean', 'jaccard'], help="the metric used for umap",
                           default="jaccard")
    args = argparser.parse_args()

    print("reading table...")
    df = pd.read_csv(args.input, sep="\t")
    print("computing umap...")
    # draw_umap(df, args.output, metric=args.metric)
    plot_umap(df, args.output)
