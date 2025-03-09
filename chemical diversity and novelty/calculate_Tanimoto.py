# Import modules
import argparse
# from rdkit import DataStructs
from sklearn.neighbors import DistanceMetric
import numpy as np
import pandas as pd
from pandas import DataFrame, Series
import collections
from multiprocessing import Pool
import multiprocessing
import itertools


def readIn(inTable):
    """
    Read input table and return a dictionary of RiPPs.

    Args:
        inTable (str): Path to the input table file.

    Returns:
        dict: A dictionary where keys are RiPP classes and values are dictionaries 
              with sequence IDs as keys and ECFP6 fingerprints as values.
    """
    dict_ripps = dict()
    with open(inTable, 'r') as fin:
        fin.readline()
        for line in fin:
            seqID, rippClass, seq, ecfp6, pca1, pca2, umap1, umap2 = line.rstrip("\n").split("\t")
            if rippClass not in dict_ripps:
                dict_ripps[rippClass] = dict()
            dict_ripps[rippClass][seqID] = ecfp6
    return dict_ripps


def runDist(inTable):
    """
    Calculate Tanimoto similarity within and between RiPP classes.

    Args:
        inTable (str): Path to the input table file.

    Returns:
        None
    """
    dict_ripps = readIn(inTable)
    # For each class, the sequence with highest median was taken as the representative of this RiPP class
    ripp_class_representative = dict()
    # For each class, average median Tanimoto coefficient was used to evaluate the similarity within RiPP class
    within_ripp_similarity = dict()
    # Get Tanimoto similarity within groups
    for rippstype, rippSeqs in dict_ripps.items():
        dict_oneRipp_dist = dict() 
        allEcfp6 = list(rippSeqs.values())
        in_dist_allEcfp6 = [list(x) for x in allEcfp6]
        # Pairwise Jaccard/Tanimoto distance
        dist = DistanceMetric.get_metric('jaccard')
        pair_dis = dist.pairwise(in_dist_allEcfp6)
        pair_dis_to_df = pd.DataFrame(pair_dis, columns=list(rippSeqs.keys()), index=list(rippSeqs.keys()))
        # Average median Tanimoto coefficient
        ave_median_tanimoto = pair_dis_to_df.median(axis=0).mean()
        within_ripp_similarity[rippstype] = ave_median_tanimoto
        # Sequence with highest median Tanimoto coefficient was taken as representative
        print(pair_dis_to_df.median(axis=0))
        representative_id = pair_dis_to_df.median(axis=0).idxmax(axis=0)
        ripp_class_representative[rippstype] = dict_ripps[rippstype][representative_id]
    # Get pairwise distance between different RiPPs
    ripps_representaive_ecfp6 = list(ripp_class_representative.values())
    ripps_representaive_ecfp6_in = [list(x) for x in ripps_representaive_ecfp6]
    ripps_pair_dis = dist.pairwise(ripps_representaive_ecfp6_in)
    ripps_pair_dis_to_df = pd.DataFrame(ripps_pair_dis, columns=list(ripp_class_representative.keys()), index=list(ripp_class_representative.keys()))
    ripps_pair_dis_to_df.to_csv("representaive_pairwise_dist.tsv", sep="\t")
    
    df_within = DataFrame.from_dict(within_ripp_similarity, orient="index")
    df_within.columns = ["within_Tanimoto"]
    df_within.to_csv("within_ripps_Tanimoto.tsv", sep="\t")


def main():
    """
    Main function to parse arguments and run distance calculation.
    """
    parse = argparse.ArgumentParser(description="Calculate Tanimoto similarity between molecules")

    parse.add_argument("--inTable", help="Input table", required=True)

    args = parse.parse_args()

    runDist(args.inTable)


if __name__ == "__main__":
    main()
