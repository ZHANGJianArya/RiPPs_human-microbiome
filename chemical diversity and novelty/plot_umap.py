#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import umap


def plot_umap(df_umap, file, title='UMAP'):
    plt.figure(figsize=(16, 16))

    #labels = ['Autoinducing','Lanthipeptide','LAP','Lassopeptide','NtoC_cyclized_peptides','Other_Known_RiPP','PMCS_KP','PMCS_NL','RiPP_like','rSAM_modified_RiPP','unknown_RiPP']
    #colors = ['#ca4656', '#ff81fb', '#009afb', '#9445fa', '#02316b', '#018f56', '#7e1606', '#989999', '#ff9123','#303841','#7f7f7f']

    #labels = ['Known_RiPP','PMCS_KP','PMCS_NL','unknown_RiPP']
    #colors = ['#018f56','#9445fa','#ca4656','#7f7f7f']
    labels = ['RiPP_from_HM','RiPP_from_MIBIG']
    colors = ['#000000','#ff001d']


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
    #plt.savefig(svg, format="svg", transparent=True)
    #plt.savefig(png, format="png", transparent=True)
    plt.savefig(svg, format="svg", transparent=False)
    plt.savefig(png, format="png", transparent=False)


def draw_umap(df, file, metric='jaccard', title='UMAP'):
    fit = umap.UMAP(metric=metric)

    ecfps = list(df["ECFP6"])
    fingerprints = [[*i] for i in ecfps]
    fingerprints = np.asarray(fingerprints, dtype=int)
    u = fit.fit_transform(fingerprints)
    umap_df = pd.DataFrame(u, columns=["umap1", "umap2"])
    df_umap = pd.concat([df, umap_df], axis=1)
    plt.figure(figsize=(24, 24))

    #labels = ['Autoinducing','Lanthipeptide','LAP','Lassopeptide','NtoC_cyclized_peptides','Other_Known_RiPP','PMCS_KP','PMCS_NL','RiPP_like','rSAM_modified_RiPP','unknown_RiPP']
    #colors = ['#ca4656', '#ff81fb', '#009afb', '#9445fa', '#02316b', '#018f56', '#7e1606', '#989999', '#ff9123','#303841','#7f7f7f']
    labels = ['RiPP_from_HM','RiPP_from_MIBIG']
    colors = ['#000000','#ff001d']   

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
    import pandas as pd

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
