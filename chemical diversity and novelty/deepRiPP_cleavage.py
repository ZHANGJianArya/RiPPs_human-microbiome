#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from numba import cuda
import os
#os.environ['CUDA_VISIBLE_DEVICES'] = '2'
from Bio import SeqIO
from nlpprecursor.classification.data import DatasetGenerator as CDG
from nlpprecursor.annotation.data import DatasetGenerator as ADG
from pathlib import Path
import json
import pandas as pd
import numpy as np
import sys

import nlpprecursor

# This allows for backwards compatibility of the pickled models.
sys.modules["protai"] = nlpprecursor
        

def get_info(fa_file):
    aa_list=[]
    fa_seq = SeqIO.parse(fa_file, "fasta")
    for record in fa_seq:
        gene_dict={}
        gene_dict['sequence']=record.seq[:-1]
        gene_dict['name']=[record.id.split()[0],record.seq]
        #gene_dict['name']=[record.id.split()[0]]
        aa_list.append(gene_dict)
    sequences=aa_list

    annot_model_dir = '/data/ZHANGJian/RiPP/RiPP_from_other_bodysites/prediction/DeepRiPP/NLPPrecursor-1.0/nlpprecursor/annotation/models/annotation'
    annot_model_path = '/data/ZHANGJian/RiPP/RiPP_from_other_bodysites/prediction/DeepRiPP/NLPPrecursor-1.0/nlpprecursor/annotation/models/annotation/model.p'
    annot_vocab_path = '/data/ZHANGJian/RiPP/RiPP_from_other_bodysites/prediction/DeepRiPP/NLPPrecursor-1.0/nlpprecursor/annotation/models/annotation/vocab.pkl'

    cleavage_predictions = ADG.predict(annot_model_path, annot_vocab_path, sequences)
    return cleavage_predictions

def check_ripp(class_p):
    ref = {}
    _sequence=[]
    start=[]
    stop=[]
    score=[]
    name=[]
    status=[]

    for class_dict in class_p: 
         _sequence.append(class_dict['cleavage_prediction']['sequence'])
         start.append(class_dict['cleavage_prediction']['start'])
         stop.append(class_dict['cleavage_prediction']['stop'])
         score.append(class_dict['cleavage_prediction']['score'])
         status.append(class_dict['cleavage_prediction']['status'])
         name.append(class_dict['cleavage_prediction']['name'])

    ref['gene_id']=name
    ref['sequence']=_sequence
    ref['start']=start
    ref['stop']=stop
    ref['score']=score  
    ref['status']=status
    ref=pd.DataFrame(ref)
    return ref


if __name__ == '__main__':
    cds_dir ="/data/ZHANGJian/RiPP/HM_prioritized/Cleavage_predictions/"
    res_dir='/data/ZHANGJian/RiPP/HM_prioritized/Cleavage_predictions/'
    for file in os.listdir(res_dir):
        if file.split('.')[-1] != "fas":continue
        if file +'.xlsx' in os.listdir(cds_dir):
            continue
        if 'empty' in os.listdir(cds_dir):
           if file in [line.rstrip('\n') for line in open(cds_dir+'empty','r')]:
              continue    
        print(file)
        a= get_info(res_dir + os.sep + file)
        if a == []:
           f=open(cds_dir+'empty','a')   
           f.write('%s\n'%file)
           f.close()
        else:
           out_file = check_ripp(a)
           writer=pd.ExcelWriter(cds_dir+"%s.xlsx"%file[:-4])
           out_file.to_excel(writer,'sheet1',index=False)
           writer.save()
        
