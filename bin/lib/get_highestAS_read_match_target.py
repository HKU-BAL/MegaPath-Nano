import pandas 
import pysam
import os
import sys

cwd=os.path.dirname(os.path.realpath(__file__))
root_dir=os.path.dirname(os.path.dirname(cwd))
data_dir=f'{root_dir}/db/amplicon'
TARGET_SEQID_FILE=f'{data_dir}/TB_seqid.txt'
target_seqid_df=pandas.read_csv(TARGET_SEQID_FILE,sep='\t',header=None,names=['sequence_id'])

in_bam=sys.argv[1]
dict={'read_id':[],'sequence_id':[],'alignment_score':[]}
with pysam.AlignmentFile(in_bam, "r" ,) as fin:
    for entry in fin:
        if entry.flag & 4==0:
            dict['read_id'].append(entry.qname)
            dict['sequence_id'].append(entry.reference_name)
            dict['alignment_score'].append(entry.get_tag('AS'))
align_list=pandas.DataFrame.from_dict(dict)

align_list_groupby_read_id=align_list.groupby(['read_id'])
align_list['alm_max'] = align_list_groupby_read_id['alignment_score'].transform(max)
align_list=align_list[align_list['alignment_score']==align_list['alm_max']]
#import pdb;pdb.set_trace()
merged=align_list.merge(target_seqid_df,on=['sequence_id'])
merged['read_id'].drop_duplicates().to_csv(f'{in_bam}.highestAS_read_match_target.list',header=None,index=False)

