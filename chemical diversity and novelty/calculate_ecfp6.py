#!/usr/bin/env python
import argparse
import pandas as pd
import subprocess
import numpy as np
import umap
from pathlib import Path
from sklearn.decomposition import PCA
from plot_pca import plot_2d_pca
from plot_umap import plot_umap


def calc_ecfp6(peptides):
    p = subprocess.run(["java", "-jar", "./peptideECFP6.jar"], encoding='utf-8', input=peptides,
                       stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    return p.stdout.split("\n")[:-1]


def calc_pca(fingerprints: np.ndarray):
    pca = PCA(n_components=2)
    result = pca.fit_transform(fingerprints)
    print("explained variance ratio: ", pca.explained_variance_ratio_)
    return result


def calc_umap(fingerprints: np.ndarray, metric):
    fit = umap.UMAP(metric=metric)
    u = fit.fit_transform(fingerprints)
    return u


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(description="calculate ECFP6 fingerprint for peptides")
    argparser.add_argument("-i", "--input", required=True,
                           help="the input table file of peptide sequences, must has a column with name 'sequence'")
    argparser.add_argument("-o", "--output", help="the output result file", required=True)
    argparser.add_argument("-p", "--plot", help="plot pca and umap", action="store_true")
    argparser.add_argument("-m", "--metric", choices=['euclidean', 'jaccard'], help="the metric used for umap",
                           default="jaccard")
    args = argparser.parse_args()

    print("reading input file...")
    df = pd.read_csv(args.input, sep="\t")
    peptides = "\n".join(list(df["sequence"]))
    print("calculating ECFP6...")
    ecfps = calc_ecfp6(peptides)
    ecfp_df = pd.DataFrame(ecfps, columns=["ECFP6"])

    fingerprints = [[*i] for i in ecfps]
    fingerprints = np.asarray(fingerprints, dtype=int)

    if args.plot:
        print("calculating PCA...")
        pca_df = pd.DataFrame(calc_pca(fingerprints), columns=["pca1", "pca2"])
        print("calculating UMAP...")
        umap_df = pd.DataFrame(calc_umap(fingerprints, args.metric), columns=["umap1", "umap2"])
        out = pd.concat([df, ecfp_df, pca_df, umap_df], axis=1)
        out.to_csv(args.output, sep="\t", index=False)
        output_dir = Path(args.output).parent
        filename = Path(args.input).stem
        print("plotting PCA...")
        plot_2d_pca(out, str(output_dir.joinpath(filename + "_pca.png")))
        print("plotting UMAP...")
        plot_umap(out, str(output_dir.joinpath(filename + "_umap_" + args.metric + ".png")))
    else:
        out = pd.concat([df, ecfp_df], axis=1)
        out.to_csv(args.output, sep="\t", index=False)
