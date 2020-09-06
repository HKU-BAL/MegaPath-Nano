#!/usr/bin/env python
# coding: utf-8
import argparse
import pandas as pd
import sys
import os
from collections import Counter
from pandarallel import pandarallel
from concurrent.futures import ThreadPoolExecutor
import json

def build_mcount_set_list(species_name_list,mcount_set_list):
    list_n=len(species_name_list)
    i=0
    j=1
    while i <list_n-1:
        while j<list_n:
            if species_name_list[i]==species_name_list[j]:
                j+=1
                continue
            mcount_set_list.append(frozenset({species_name_list[i],species_name_list[j]}))
            j+=1
        i+=1
        j=i+1
def Count(groupby_count,species_i):
    try:
        return groupby_count[species_i]
    except KeyError:
        return 0

def MCount(species_i,species_j,mcount_set_list):
    return Counter(mcount_set_list)[frozenset({species_i,species_j})]

def build_i_explains_j_dict(species_i,species_j,name_groupby_count,unique_alignment_name_groupby_count,Counter_mcount_set_list,i_explains_j_dict,ratio,error_rate):
    #print(i_explains_j_dict.values())
    #if species_i in i_explains_j_dict.values():
    #    pass
    #else:
    AllCount_i=Count(name_groupby_count,species_i)
    UCount_i=Count(unique_alignment_name_groupby_count,species_i)
    UCount_j=Count(unique_alignment_name_groupby_count,species_j)
    MCount_i_j=Counter_mcount_set_list[frozenset({species_i,species_j})]
    if AllCount_i+MCount_i_j>=ratio*AllCount_i and UCount_j<error_rate*UCount_i:
        if species_i in i_explains_j_dict:
            i_explains_j_dict[species_i].add(species_j)
        else:
            i_explains_j_dict[species_i]={species_j}

def iterate_reassign(align_list,i_explains_j_dict,AS_threshold):
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

            def first_AS_passed(x,AS_threshold):
                return True if x['alignment_score'] * AS_threshold <= x['first_alignment_score'] else False

            search_pair['first_AS_passed'] = search_pair[['alignment_score','first_alignment_score']].apply(first_AS_passed,args=(AS_threshold),axis=1)
            search_pair=search_pair[search_pair["first_AS_passed"]]

            align_list.loc[search_pair.index,['name','sequence_id']]=[name,name+"_RA"]
            align_list.loc[search_pair.index].to_csv('alignlist_update.csv',mode='a',header=None)

    align_list['is_in_explain_other'] = align_list['name'].parallel_apply(is_in_explain_other)

    uk = align_list[align_list['is_in_explain_other']]
    #  sort by ranking order of name
    species_ranking_list=uk['name'].value_counts().keys().tolist()
    uk['name']=pd.Categorical(uk['name'],categories=species_ranking_list,ordered=True)
    uk=uk.sort_values(ascending=True,by=['name'])
    uk[['read_id','name','alignment_score','sequence_id']].parallel_apply( second_index ,axis=1)
    return align_list

FLAGS=None
def Reassign(align_list,iteration=1,error_rate=0.05,ratio=0.05,threads=96,AS_threshold=0):
    #TODO iterate groupby and explainlist
    #TODO generate
    TAXON_FILE="/mnt/bal13/wwlui/Nanopath_TB/db_new_refseq.sequence_name.first_two_words"
    taxon_df=pd.read_csv(TAXON_FILE,sep='\t',header=None,names=['sequence_id','name'])
    align_list=align_list.merge(right=taxon_df,on=['sequence_id'],how='inner')

    raw_align_list=align_list

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

    with ThreadPoolExecutor(threads) as exec:
        for read_id in mcount_read_id_list:
            species_name_list=read_id_groupby.get_group(read_id)['name'].tolist()
            exec.submit(build_mcount_set_list, species_name_list,mcount_set_list)
    Counter_mcount_set_list=Counter(mcount_set_list)
    i_explains_j_dict={}
    with ThreadPoolExecutor(threads) as exec:
        for species_i in species_list:
            for species_j in species_list:
                #print(i_explains_j_dict.values())
                if species_i==species_j or species_i in i_explains_j_dict.values():
                #if species_i==species_j :
                    continue
                exec.submit(build_i_explains_j_dict, species_i,species_j,name_groupby_count,unique_alignment_name_groupby_count,Counter_mcount_set_list,i_explains_j_dict,ratio,error_rate)
    
    if i_explains_j_dict==False:
        return raw_align_list

    with open('i_explains_j_dict', 'w') as file:
        file.write(str(i_explains_j_dict))

    pandarallel.initialize()

    for i in range(iteration):
        print('Iteration:',i)
        align_list=iterate_reassign(align_list,i_explains_j_dict,AS_threshold)
        #should be no need to update align_list
        if os.path.exists('alignlist_update.csv'):
            update=pd.read_csv('alignlist_update.csv',names=['read_id', 'read_length', 'read_from', 'read_to', 'strand', 'sequence_id', 'sequence_length', 'sequence_from', 'sequence_to','match','alignment_block_length', 'mapq', 'edit_dist', 'alignment_score','assembly_id', 'tax_id', 'species_tax_id', 'genus_tax_id', 'alignment_score_tiebreaker','name','is_in_explain_other'])
            raw_align_list.update(update[~update.index.duplicated(keep='first')])
            os.remove('alignlist_update.csv')
        else:
            break


    raw_align_list=raw_align_list.astype({'read_id': str, 'read_length': int, 'read_from': int, 'read_to': int, 'strand': str, 'sequence_id': str, 'sequence_length': int, 'sequence_from' : int, 'sequence_to': int,'match': int, 'mapq': int, 'edit_dist': int, 'alignment_score': int,'assembly_id': str, 'tax_id': int, 'species_tax_id': int, 'genus_tax_id': int, 'alignment_score_tiebreaker': float,'alignment_block_length': int,'name': str})
    raw_align_list.to_csv('alignlist_reassigned.csv')
    return align_list


if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description='MegaPath-Nano: Reassignment')
    parser.add_argument('--align_list', required=True,help='Input align list')
    parser.add_argument('--AS_threshold', default=0.0, help='Percentage of AS required for an alignment to be reassigned to another species from the original species')
    parser.add_argument('--ratio', default=0.05, help='Ratio of dissimilarity between species')
    parser.add_argument('--error_rate', default=0.05, help='Allowed error rate')
    parser.add_argument('--threads', default=48, help='Max num of threads')
    parser.add_argument('--iteration', default=1, help='Num of iterations')
    FLAGS = parser.parse_args()
    align_list=pd.read_csv(FLAGS.align_list, index_col =0)
    Reassign(align_list,iteration=FLAGS.iteration,error_rate=FLAGS.error_rate,ratio=FLAGS.ratio,threads=FLAGS.threads,AS_threshold=FLAGS.AS_threshold)





