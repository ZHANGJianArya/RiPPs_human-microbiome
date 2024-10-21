
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


# main body
def readIn(inTable):
    dict_ripps = dict()
    with open(inTable, 'r') as fin:
        fin.readline()
        for line in fin:
            seqID, rippClass, seq, ecfp6, pca1, pca2, umap1, umap2 = line.rstrip("\n").split("\t")
            if rippClass not in dict_ripps:
                dict_ripps[rippClass] = dict()
            dict_ripps[rippClass][seqID] = ecfp6
    return dict_ripps

'''
def cldist(parameters):
    item = parameters[0]
    query_seq = parameters[1]
    subject_seq = parameters[2]
    list_dist = []
    for y in subject_seq:
            dist = DataStructs.FingerprintSimilarity(str.encode(query_seq),str.encode(y))
            list_dist.append(dist)
    dict_item = dict()
    dict_item[item] = list_dist
    return dict_item
'''


def runDist(inTable):
    dict_ripps = readIn(inTable)
    # For each class, the sequence with highest median was taken as the the representative of this RiPP class
    ripp_class_representative = dict()
    # For each class, average median Tanimoto coeffcient was used to evaluate the similarity within RiPP class
    within_ripp_similarity = dict()
    # get tanimoto similarity within groups
    for rippstype, rippSeqs in dict_ripps.items():
        dict_oneRipp_dist = dict() 
        allEcfp6 = list(rippSeqs.values())
        in_dist_allEcfp6 = [list(x) for x in allEcfp6]
        # pairwise jaccard/tanimoto distance
        dist = DistanceMetric.get_metric('jaccard')
        pair_dis = dist.pairwise(in_dist_allEcfp6)
        pair_dis_to_df = pd.DataFrame(pair_dis, columns = list(rippSeqs.keys()), index = list(rippSeqs.keys()))
        # average meadian Tanimoto coeffcient
        ave_median_tanimoto = pair_dis_to_df.median(axis=0).mean()
        within_ripp_similarity[rippstype] = ave_median_tanimoto
        # seq with highest meadian Tanimoto coeffcient was taken as representative
        print(pair_dis_to_df.median(axis=0))
        #representative_id = pair_dis_to_df.median(axis=0).max(axis=0)
        representative_id = pair_dis_to_df.median(axis=0).idxmax(axis=0)
        ripp_class_representative[rippstype] = dict_ripps[rippstype][representative_id]
    # get pairwise distance between different RiPPs
    ripps_representaive_ecfp6 = list(ripp_class_representative.values())
    ripps_representaive_ecfp6_in = [list(x) for x in ripps_representaive_ecfp6]
    ripps_pair_dis = dist.pairwise(ripps_representaive_ecfp6_in)
    ripps_pair_dis_to_df = pd.DataFrame(ripps_pair_dis, columns = list(ripp_class_representative.keys()), index = list(ripp_class_representative.keys()))
    ripps_pair_dis_to_df.to_csv("representaive_pairwise_dist.tsv", sep="\t")
    
    df_within = DataFrame.from_dict(within_ripp_similarity, orient="index")
    df_within.columns = ["within_Tanimoto"]
    df_within.to_csv("within_ripps_Tanimoto.tsv", sep="\t")


def main():
    parse = argparse.ArgumentParser(description="calculate Tanimoto similarity between moleculars")

    parse.add_argument("--inTable", help="input table", required=True)

    args = parse.parse_args()

    runDist(args.inTable)


if __name__ == "__main__": 
    main()

