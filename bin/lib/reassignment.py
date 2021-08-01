#!/usr/bin/env python
# coding: utf-8
import argparse
import pandas as pd
import sys
import os
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import psutil
import functools
from multiprocessing import Pool
import pickle
import itertools

def get_combinations_counter(species_name_list):
    return Counter(itertools.combinations(species_name_list,2))

def Count(groupby_count,species_i):
    try:
        return groupby_count[species_i]
    except KeyError:
        return 0

def MCount(species_i,species_j,mcount_set_list):
    return Counter(mcount_set_list)[frozenset({species_i,species_j})]

def build_i_explains_j_dict(species_i,species_j,name_groupby_count,unique_alignment_name_groupby_count,Counter_mcount_set_list,i_explains_j_dict,ratio,error_rate):
    AllCount_i=Count(name_groupby_count,species_i)
    UCount_i=Count(unique_alignment_name_groupby_count,species_i)
    UCount_j=Count(unique_alignment_name_groupby_count,species_j)
    MCount_i_j=Counter_mcount_set_list[frozenset({species_i,species_j})]
    if AllCount_i-MCount_i_j>=ratio*AllCount_i and UCount_j<error_rate*UCount_i:
        if species_i in i_explains_j_dict:
            i_explains_j_dict[species_i].add(species_j)
        else:
            i_explains_j_dict[species_i]={species_j}

def reassign_alignment(align_list,i_explains_j_dict,AS_threshold):
    def is_in_explain_other(x):
        return True if x in i_explains_j_dict.keys() else False
    def second_index(x):
        read_id = x['read_id']
        name = x['name']
        first_alignment_score = x['alignment_score']
        search_pair = align_list[align_list['read_id'] == read_id]
        def second_exist(x):
            return True if (name in i_explains_j_dict.keys() and  x in i_explains_j_dict[name] and x != name) else False
        search_pair['second_exist'] = search_pair['name'].apply(second_exist)
        search_pair=search_pair[search_pair['second_exist']]
        if not search_pair.empty:
            search_pair['first_alignment_score'] = first_alignment_score
            def first_AS_passed(x):
                return True if x['alignment_score'] * AS_threshold <= x['first_alignment_score'] else False
            search_pair['first_AS_passed'] = search_pair[['alignment_score','first_alignment_score']].apply(first_AS_passed,axis=1)
            search_pair=search_pair[search_pair['first_AS_passed']]
            align_list.loc[search_pair.index,['name','sequence_id']]=[name,f'{name}_reassigned']
    align_list['is_in_explain_other'] = align_list['name'].apply(is_in_explain_other)
    uniq_key = align_list[align_list['is_in_explain_other']]
    #  sort by ranking order of name
    species_ranking_list=uniq_key['name'].value_counts().keys().tolist()
    uniq_key.loc[:,'name']=pd.Categorical(uniq_key['name'],categories=species_ranking_list,ordered=True)
    uniq_key=uniq_key.sort_values(ascending=True,by=['name'])
    uniq_key[['read_id','name','alignment_score','sequence_id']].apply( second_index ,axis=1)
    return align_list

def Reassign(align_list,db_folder,error_rate=0.05,ratio=0.05,threads=24,AS_threshold=0,level='species'):
    #TODO iterate groupby and explainlist
    taxon_df=pd.read_csv(f'{db_folder}/sequence_name', sep='\t',header=None,names=['sequence_id','name'])
    if level=='species':
        taxon_df['name']=taxon_df['name'].apply(lambda x: " ".join(x.split(" ",2)[0:2]) if ' sp. ' not in x else " ".join(x.split(" ",3)[0:3]))
    align_list=align_list.merge(right=taxon_df,on=['sequence_id'],how='inner')
    #  keep the highest AS for same read_id & name
    align_list=align_list.sort_values(by=['alignment_score']).drop_duplicates(subset=['read_id','name'],keep='last')

    read_id_groupby=align_list.groupby(['read_id'])
    #  higher abun explains first
    name_groupby_count=align_list.groupby(['name']).size().sort_values(ascending=False)
    species_list=name_groupby_count.index
    name_groupby_count=name_groupby_count.to_dict()
    unique_alignment_align_list=align_list.drop_duplicates(subset='read_id',keep=False)
    unique_alignment_name_groupby_count=unique_alignment_align_list.groupby(['name']).size().to_dict()

    mcount_set_list=[]
    mcount_read_id_list=read_id_groupby.size()[read_id_groupby.size()>1].index

    with Pool(processes=threads) as exec:
        species_name_list=[read_id_groupby.get_group(read_id)['name'].tolist() for read_id in mcount_read_id_list]
        results=exec.map(get_combinations_counter, species_name_list,chunksize=768)
        #results=exec.imap_unordered(build_mcount_set_list, species_name_list,chunksize=768)
        #mcount_set_list=[mcount_set for result in results for mcount_set in result]
    Counter_mcount_set_list=functools.reduce(lambda a,b : a+b,results)
    i_explains_j_dict={}
    with ThreadPoolExecutor(threads) as exec:
        for species_i in species_list:
            for species_j in species_list:
                if species_i==species_j or species_i in i_explains_j_dict.values():
                    continue
                exec.submit(build_i_explains_j_dict, species_i,species_j,name_groupby_count,unique_alignment_name_groupby_count,Counter_mcount_set_list,i_explains_j_dict,ratio,error_rate)

    if not i_explains_j_dict:
        return align_list

    with open('i_explains_j_dict.pickle', 'wb') as file:
        pickle.dump(i_explains_j_dict,file)

    align_list=reassign_alignment(align_list,i_explains_j_dict,AS_threshold)
    align_list.to_csv('alignlist_reassigned.csv')
    return align_list


if __name__ == '__main__': 
    NANO_DIR=os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
    parser = argparse.ArgumentParser(description='MegaPath-Nano: Reassignment')
    parser.add_argument('--align_list', required=True,help='Input align list')
    parser.add_argument('--level', default='species',choices=['species', 'strain'],help='Species- or strain- level resolution')
    parser.add_argument('--AS_threshold', default=0.0, type=float,help='Percentage of AS required for an alignment to be reassigned to another species from the original species')
    parser.add_argument('--ratio', default=0.05, type=float,help='Ratio of dissimilarity between species')
    parser.add_argument('--error_rate', default=0.05, type=float,help='Allowed error rate')
    parser.add_argument('--threads', default=psutil.cpu_count(logical=True), type=int, help='Num of threads')
    parser.add_argument('--db_folder', default=f'{NANO_DIR}/db', help='Db folder')
    FLAGS = parser.parse_args()
    align_list_col_type = {'read_id': str, 'read_length': int, 'read_from': int, 'read_to': int, 'strand': str, 'sequence_id': str, 'sequence_length': int, 'sequence_from' : int, 'sequence_to': int, 'match': int, 'mapq': int, 'edit_dist': int, 'alignment_score': int, 'assembly_id': str, 'tax_id': int, 'species_tax_id': int, 'genus_tax_id': int, 'alignment_score_tiebreaker': float,'alignment_block_length':int}
    align_list=pd.read_csv(FLAGS.align_list, index_col =0,dtype=align_list_col_type)
    Reassign(align_list,db_folder=FLAGS.db_folder,error_rate=FLAGS.error_rate,ratio=FLAGS.ratio,threads=FLAGS.threads,AS_threshold=FLAGS.AS_threshold,level=FLAGS.level)





