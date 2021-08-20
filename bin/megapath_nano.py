#!/usr/bin/env python3
import os
import sqlite3
import argparse
import tempfile
import subprocess
import shlex
import psutil
import shutil
import random
import atexit
import time
import inspect
import warnings
import numpy
import pandas
import pybedtools
from pybedtools.bedtool import BedTool      # all bed should have 'chrom' == sequence_id, 'name' == assembly_id

from lib.aligner import Align
from lib.reassignment import Reassign
from lib.assembly_metadata import AssemblyMetadata

FLAGS = None
UNPARSED = None


assembly_info_col_name = [ 'assembly_id', 'bioproject', 'biosample', 'wgs_master', 'refseq_category', 'taxid', 'species_taxid', 'organism_name', 'infraspecific_name',
                           'isolate', 'version_status', 'assembly_level', 'release_type', 'genome_rep', 'seq_rel_date', 'asm_name', 'submitter', 'gbrs_paired_asm', 'paired_asm_comp',
                           'ftp_path', 'excluded_from_refseq', 'relation_to_type_material']
assembly_info_col_type = { 'assembly_id': str, 'bioproject': str, 'biosample': str, 'wgs_master': str, 'refseq_category': str, 'taxid': int, 'species_taxid': int, 'organism_name': str, 'infraspecific_name': str,
                           'isolate': str, 'version_status': str, 'assembly_level': str, 'release_type': str, 'genome_rep': str, 'seq_rel_date': str, 'asm_name': str, 'submitter': str, 'gbrs_paired_asm': str, 'paired_asm_comp': str,
                           'ftp_path': str, 'excluded_from_refseq': str, 'relation_to_type_material': str}

read_info_col_name = ['read_id', 'read_length', 'average_quality', 'read_length_trimmed', 'average_quality_trimmed', 'passed_filter']
read_info_col_type = {'read_id': str, 'read_length': int, 'average_quality': float, 'read_length_trimmed': int, 'average_quality_trimmed': float, 'passed_filter': int}

aligned_read_list_col_name = ['read_id', 'read_length', 'aligned', 'human_read', 'decoy_read', 'microbe_read', 'unaligned']
aligned_read_list_col_type = {'read_id': str, 'read_length': int, 'aligned': int, 'human_read': int, 'decoy_read': int, 'microbe_read': int, 'unaligned': int}

consolidated_read_list_col_name = ['read_id', 'read_length', 'average_quality', 'read_length_trimmed', 'average_quality_trimmed', 'passed_filter', 'aligned', 'human_read', 'decoy_read', 'microbe_read']
consolidated_read_list_col_type = {'read_id': str, 'read_length': int, 'average_quality': float, 'read_length_trimmed': int, 'average_quality_trimmed': float, 'passed_filter': int, 'aligned': int, 'human_read': int, 'decoy_read': int, 'microbe_read': int}

read_stat_col_name = ['total_number_of_read', 'passed_filter', 'aligned', 'human_read', 'decoy_read', 'microbe_read', 'unaligned', 
                      'total_read_bp', 'total_passed_filter_read_bp', 'total_aligned_read_bp', 'total_human_read_bp', 'total_decoy_read_bp', 'total_microbe_read_bp', 'total_unaligned_read_bp']
read_stat_col_type = {'total_number_of_read': int, 'passed_filter': int, 'aligned': int, 'human_read': int, 'decoy_read': int, 'microbe_read': int, 
                      'total_read_bp': int, 'total_passed_filter_read_bp': int, 'total_aligned_read_bp': int, 'total_human_read_bp': int, 'total_decoy_read_bp': int, 'total_microbe_read_bp': int,  'total_unaligned_read_bp': int}

quality_score_histogram_col_name = ['quality_score_bin', 'read_count']
quality_score_histogram_col_type = {'quality_score_bin': float, 'read_count': int}

read_length_histogram_col_name = ['read_length_bin', 'read_count']
read_length_histogram_col_type = {'read_length_bin': int, 'read_count': int}

similar_species_marker_name = ['similar_genus_tax_id', 'is_similar', 'similarity_cutoff_1', 'aligned_percent_1', 'combine_logic', 'similarity_cutoff_2', 'aligned_percent_2']
similar_species_marker_type = {'similar_genus_tax_id': int, 'is_similar': int, 'similarity_cutoff_1': int, 'aligned_percent_1': float, 'combine_logic': str, 'similarity_cutoff_2': int, 'aligned_percent_2': float}

similar_species_marker_name_with_assembly_id = ['assembly_id'] + similar_species_marker_name
similar_species_marker_type_with_assembly_id = similar_species_marker_type.copy()
similar_species_marker_type_with_assembly_id.update({'assembly_id': str})

align_stat_unique_col_name = ['unique_total_aligned_bp', 'unique_average_depth', 'unique_covered_percent']
align_stat_unique_col_type = {'unique_total_aligned_bp': int, 'unique_average_depth': float, 'unique_covered_percent': float}
align_stat_pre_noise_col_name = ['pre_total_number_of_read', 'pre_total_read_bp', 'pre_average_read_length',
                                 'pre_total_aligned_bp', 'pre_average_depth', 'pre_covered_percent', 
                                 'pre_average_identity', 'pre_average_edit_dist', 'pre_average_alignment_score']
align_stat_pre_noise_col_type = {'pre_total_number_of_read': int, 'pre_total_read_bp': int, 'pre_average_read_length': float,
                                 'pre_total_aligned_bp': int, 'pre_average_depth': float, 'pre_covered_percent': float, 
                                 'pre_average_identity': float, 'pre_average_edit_dist': float, 'pre_average_alignment_score': float}
align_stat_raw_col_name = ['total_number_of_read', 'total_read_bp', 'average_read_length',
                           'total_aligned_bp', 'average_depth', 'covered_percent', 
                           'average_identity', 'average_edit_dist', 'average_alignment_score']
align_stat_raw_col_type = {'total_number_of_read': int, 'total_read_bp': int, 'average_read_length': float,
                           'total_aligned_bp': int, 'average_depth': float, 'covered_percent': float, 
                           'average_identity': float, 'average_edit_dist': float, 'average_alignment_score': float}
align_stat_raw_col_name_by_assembly_id = ['species_tax_id', 'species_tax_name', 'assembly_id', 'assembly_length'] + align_stat_raw_col_name + ['tax_id', 'tax_name']
align_stat_raw_col_type_by_assembly_id = {'species_tax_id': int, 'species_tax_name': str, 'assembly_id': str, 'assembly_length': int}
align_stat_raw_col_type_by_assembly_id.update(align_stat_raw_col_type)
align_stat_raw_col_type_by_assembly_id.update({'tax_id': int, 'tax_name': str})
align_stat_raw_col_name_by_sequence_id = ['sequence_id', 'sequence_name', 'sequence_length'] + align_stat_raw_col_name
align_stat_raw_col_type_by_sequence_id = {'sequence_id': str, 'sequence_name': str, 'sequence_length': int}
align_stat_raw_col_type_by_sequence_id.update(align_stat_raw_col_type)

#added
align_stat_raw_col_name_by_sequence_id_with_assembly_id = ['assembly_id', 'assembly_name', 'assembly_length'] + ['sequence_id', 'sequence_name', 'sequence_length', 'sequence_total_aligned_bp', 'sequence_average_depth', 'sequence_covered_percent', 'assembly_adjusted_total_aligned_bp', 'assembly_adjusted_average_depth', 'assembly_adjusted_covered_percent'] + ['species_tax_id', 'species_tax_name', 'tax_id', 'tax_name', 'genus_tax_id', 'genus_tax_name']
align_stat_raw_col_type_by_sequence_id_with_assembly_id = {'species_tax_id': int, 'species_tax_name': str, 'assembly_id': str, 'assembly_name':str, 'assembly_length': int, 'sequence_id':str, 'sequence_name':str, 'sequence_length':int, 'sequence_total_aligned_bp':int, 'sequence_average_depth':float, 'sequence_covered_percent':float, 'assembly_adjusted_total_aligned_bp': int, 'assembly_adjusted_average_depth': float, 'assembly_adjusted_covered_percent': float, 'tax_id': int, 'tax_name': str, 'genus_tax_id': int, 'genus_tax_name': str}
#end

align_stat_col_name = ['adjusted_total_aligned_bp', 'adjusted_average_depth', 'adjusted_covered_percent'] + align_stat_raw_col_name
align_stat_col_type = {'adjusted_total_aligned_bp': int, 'adjusted_average_depth': float, 'adjusted_covered_percent': float}
align_stat_col_type.update(align_stat_raw_col_type)
align_stat_col_name_by_assembly_id = ['species_tax_id', 'species_tax_name', 'assembly_id', 'assembly_length'] + align_stat_col_name + ['tax_id', 'tax_name', 'genus_tax_id', 'genus_tax_name']
align_stat_col_type_by_assembly_id = {'species_tax_id': int, 'species_tax_name': str, 'assembly_id': str, 'assembly_length': int}
align_stat_col_type_by_assembly_id.update(align_stat_col_type)
align_stat_col_type_by_assembly_id.update({'tax_id': int, 'tax_name': str, 'genus_tax_id': int, 'genus_tax_name': str})
align_stat_col_name_with_pre_noise_by_assembly_id = align_stat_col_name_by_assembly_id + align_stat_pre_noise_col_name
align_stat_col_type_with_pre_noise_by_assembly_id = align_stat_col_type_by_assembly_id.copy()
align_stat_col_type_with_pre_noise_by_assembly_id.update(align_stat_pre_noise_col_type)

align_stat_col_name_with_pre_noise_with_similar_species_marker_by_assembly_id = align_stat_col_name_with_pre_noise_by_assembly_id + similar_species_marker_name
align_stat_col_type_with_pre_noise_with_similar_species_marker_by_assembly_id = align_stat_col_type_with_pre_noise_by_assembly_id.copy()
align_stat_col_type_with_pre_noise_with_similar_species_marker_by_assembly_id.update(similar_species_marker_type)


align_stat_col_name_dot_report = ['species_tax_id', 'species_tax_name', 'adjusted_total_aligned_bp', 'adjusted_average_depth', 'adjusted_covered_percent']
align_stat_col_type_dot_report = {'species_tax_id': int, 'species_tax_name': str, 'adjusted_total_aligned_bp': int, 'adjusted_average_depth': float, 'adjusted_covered_percent': float}
align_list_col_name_bam_filter = ['read_id', 'sequence_id', 'sequence_from', 'sequence_to', 'assembly_id']
align_list_col_type_bam_filter = {'read_id': str, 'sequence_id': str, 'sequence_from': int, 'sequence_to': int, 'assembly_id': str}

align_list_col_name = ['read_id', 'read_length', 'read_from', 'read_to', 'strand', 'sequence_id', 'sequence_length', 'sequence_from', 'sequence_to',
                       'match', 'mapq', 'edit_dist', 'alignment_score',
                       'assembly_id', 'tax_id', 'species_tax_id', 'genus_tax_id', 'alignment_score_tiebreaker']
align_list_col_type = {'read_id': str, 'read_length': int, 'read_from': int, 'read_to': int, 'strand': str, 'sequence_id': str, 'sequence_length': int, 'sequence_from' : int, 'sequence_to': int,
                       'match': int, 'mapq': int, 'edit_dist': int, 'alignment_score': int,
                       'assembly_id': str, 'tax_id': int, 'species_tax_id': int, 'genus_tax_id': int, 'alignment_score_tiebreaker': float}

noise_projection_aligned_bp_name = ['assembly_id', 'as_noise_source', 'as_noise_target', 'projected_noise_aligned_bp']
noise_projection_aligned_bp_type = {'assembly_id': str, 'as_noise_source': int, 'as_noise_target': int, 'projected_noise_aligned_bp': int}

noise_detection_stat_col_name = ['noise_span_bp', 'noise_span_percent', 'spike_span_bp', 'spike_span_percent', 'human_span_bp', 'human_span_percent', 'microbe_span_bp', 'microbe_span_percent', 'closing_spike_span_bp', 'closing_spike_span_percent', ]
noise_detection_stat_col_type = {'noise_span_bp': int, 'noise_span_percent': float, 'human_span_bp': int, 'human_span_percent': float, 
                                 'spike_span_bp': int, 'spike_span_percent': float, 'microbe_span_bp': int, 'microbe_span_percent': float, 'closing_spike_span_bp': int, 'closing_spike_span_percent': float}
noise_detection_stat_col_name_by_assembly_id = ['assembly_id',] + noise_detection_stat_col_name
noise_detection_stat_col_type_by_assembly_id = {'assembly_id': str}
noise_detection_stat_col_type_by_assembly_id.update(noise_detection_stat_col_type)

noise_removal_stat_col_name = ['noise_total_number_of_read', 'noise_total_read_bp', 'noise_total_aligned_bp', 'spike_total_number_of_read', 'spike_total_read_bp', 'spike_total_aligned_bp', 
                               'human_total_number_of_read', 'human_total_read_bp', 'human_total_aligned_bp', 'microbe_total_number_of_read', 'microbe_total_read_bp', 'microbe_total_aligned_bp',
                               'closing_spike_total_number_of_read', 'closing_spike_total_read_bp', 'closing_spike_total_aligned_bp',
                               'short_total_number_of_read', 'short_total_read_bp', 'short_total_aligned_bp', 'all_total_number_of_read', 'all_total_read_bp', 'all_total_aligned_bp']
noise_removal_stat_col_type = {'noise_total_number_of_read': int, 'noise_total_read_bp': int, 'noise_total_aligned_bp': int, 'spike_total_number_of_read': int, 'spike_total_read_bp': int, 'spike_total_aligned_bp': int,
                               'human_total_number_of_read': int, 'human_total_read_bp': int, 'human_total_aligned_bp': int, 'microbe_total_number_of_read': int, 'microbe_total_read_bp': int, 'microbe_total_aligned_bp': int,
                               'short_total_number_of_read': int, 'short_total_read_bp': int, 'short_total_aligned_bp': int, 'all_total_number_of_read': int, 'all_total_read_bp': int, 'all_total_aligned_bp': int}
noise_removal_stat_col_name_by_assembly_id = ['assembly_id',] + noise_removal_stat_col_name
noise_removal_stat_col_type_by_assembly_id = {'assembly_id': str}
noise_removal_stat_col_type_by_assembly_id.update(noise_removal_stat_col_type)

noise_stat_col_name_with_tax_name = ['species_tax_id', 'species_tax_name', 'assembly_id', 'assembly_length'] + noise_detection_stat_col_name + noise_removal_stat_col_name
noise_stat_col_type_with_tax_name = {'species_tax_id': int, 'species_tax_name': str, 'assembly_id': str, 'assembly_length': int}
noise_stat_col_type_with_tax_name.update(noise_detection_stat_col_type)
noise_stat_col_type_with_tax_name.update(noise_removal_stat_col_type)

similar_noise_span_bp_col_name = ['target_assembly_id', 'source_assembly_id', 'microbe_span_bp']
similar_noise_span_bp_col_type = {'target_assembly_id': str, 'source_assembly_id': str, 'microbe_span_bp': int}
similar_noise_span_bp_col_name_with_tax_name = ['target_species_tax_id', 'target_species_tax_name', 'target_assembly_length', 'source_species_tax_id', 'source_species_tax_name', 'source_assembly_length'] + similar_noise_span_bp_col_name + ['total_similar_noise_bp','total_similar_noise_percent']
similar_noise_span_bp_col_type_with_tax_name = {'target_species_tax_id': int, 'target_species_tax_name': str, 'target_assembly_length': int, 'source_species_tax_id': int, 'source_species_tax_name': str, 'source_assembly_length': int}
similar_noise_span_bp_col_type_with_tax_name.update(similar_noise_span_bp_col_type)
similar_noise_span_bp_col_type_with_tax_name.update({'total_similar_noise_bp': int, 'total_similar_noise_percent': float})

noise_source_stat_col_name = ['target_assembly_id', 'source_assembly_id', 'noise_read_count', 'noise_aligned_bp', 'noise_read_bp']
noise_source_stat_col_type = {'target_assembly_id': str, 'source_assembly_id': str, 'noise_read_count': int, 'noise_aligned_bp': int, 'noise_read_bp': int}
noise_source_stat_col_name_with_tax_name = ['target_species_tax_id', 'target_species_tax_name', 'target_assembly_length', 'source_species_tax_id', 'source_species_tax_name', 'source_assembly_length'] + noise_source_stat_col_name
noise_source_stat_col_type_with_tax_name = {'target_species_tax_id': int, 'target_species_tax_name': str, 'target_assembly_length': int, 'source_species_tax_id': int, 'source_species_tax_name': str, 'source_assembly_length': int}
noise_source_stat_col_type_with_tax_name.update(noise_source_stat_col_type)


class Log:
    def __init__(self, log_file, program_name):
        self.log_file = log_file
        self.missing_data = False
        self.program_name = program_name

    def print(self, message):
        print(self.program_name + ':', inspect.stack()[1][3] + '() :', message, file=self.log_file, flush=True)

    def print_time(self):
        print(self.program_name + ':', time.strftime('%Y-%m-%d %H:%M:%S'), file=self.log_file, flush=True)

    def print_missing_data(self, data_list, data_list_with_missing_data, data_name, message):
        self.missing_data = True
        missing_data_list = list(set(data_list[data_name]) - set(data_list_with_missing_data[data_name]))
        if not missing_data_list:
            print(self.program_name + ':', inspect.stack()[1][3] + '() :', message + ':', 'missing', data_name, ':', 'not found', file=self.log_file, flush=True)
        for missing_data in missing_data_list:
            print(self.program_name + ':', inspect.stack()[1][3] + '() :', message + ':', 'missing', data_name, ':', missing_data, file=self.log_file, flush=True)


class Placeholder:
    pass

class InputContainer:
    pass

class OutputContainer:
    pass

class IO_Container:
    def __init__(self):
        self.I = InputContainer()
        self.O = OutputContainer()


def exit_with_cleanup(temp_root, RAM_root):

    if temp_root is not None:
        try:
            shutil.rmtree(temp_root)
        except (FileNotFoundError):
            pass

    if RAM_root is not None:
        try:
            shutil.rmtree(RAM_root)
        except (FileNotFoundError):
            pass

def remove_target_from_align_list(*, align_list,blast_bed):
    #TODO blast id score variable
    bed = BedTool.from_dataframe(align_list[['sequence_id', 'sequence_from', 'sequence_to','read_id','read_length']])
    intersect_bed=bed.intersect(blast_bed,wo=True)
    try:
        intersect_df=pandas.read_csv(intersect_bed.fn,sep='\t',header=None)
    except Exception as e:
        print(e)
        return align_list
    intersect_df.rename(columns={0:'align_sequence_id',1:'align_start',2:'align_end',3:'read_id',4:'read_length',5:'TB_MO_sequence_id',6:'TB_MO_align_start',7:'TB_MO_align_end',8:'num_overlap_bp'},inplace=True)
    TB_read_id=intersect_df[intersect_df['num_overlap_bp'].astype('int32')>=50]['read_id'].drop_duplicates()
    pybedtools.cleanup()
    return align_list[~align_list['read_id'].isin(TB_read_id)]

def similarity_option(*,
                      divergence):

        # 0.4 divergence could be properly supported with some modification
    if divergence == 0.8:
        # This is not ideal but minimap2 requires 2 * (O + E) + A < 128 and O must be positive
        return '-k15 -w15 -A1 -B124 -O1 -E62 -s50 -z50'
    elif divergence == 1:
        # This is not ideal but minimap2 requires 2 * (O + E) + A < 128
        return '-k15 -w15 -A1 -B99 -O4 -E59 -s50 -z50'
    elif divergence == 2:
        # This is not ideal but minimap2 requires 2 * (O + E) + A < 128
        return '-k15 -w15 -A1 -B49 -O39 -E24 -s50 -z50'
    elif divergence == 5:
        return '-k15 -w15 -A1 -B19 -O39 -E9 -s50 -z50'
    elif divergence == 10:
        return '-k15 -w15 -A1 -B9 -O16 -E4 -s50 -z50'
    elif divergence == 20:
        return '-k15 -w15 -A1 -B4 -O9 -E2 -s50 -z50'
    else:
        raise RuntimeError('Divergence not supported')


def align_list_to_best_align_list(*,
                                  assembly_metadata,
                                  log,
                                  align_list,
                                  noise_bed=None):

    max_alignment_score = align_list[['read_id', 'alignment_score']].groupby(['read_id'], as_index=False).max().rename(columns={'alignment_score': 'max_alignment_score'})

    assembly_best_align_list = align_list.sort_values(['read_id', 'assembly_id', 'alignment_score', 'alignment_score_tiebreaker']).drop_duplicates(subset=['read_id', 'assembly_id'], keep='last')
    align_list_with_max_alignment_score = assembly_best_align_list.merge(
                                                                         right=max_alignment_score.set_index(['read_id']),
                                                                         how='inner', 
                                                                         left_on='read_id', 
                                                                         right_index = True,
                                                                         suffixes=['', '_y'],
                                                                         validate='m:1',
                                                                        ).query('alignment_score == max_alignment_score').drop(['max_alignment_score'], axis=1)
    max_score_count = align_list_with_max_alignment_score[['read_id']].assign(max_score_count = lambda x: 1).groupby(['read_id'], as_index=False).sum()

    align_list_with_max_score_count = align_list_with_max_alignment_score.merge(
                                                                                right=max_score_count.set_index(['read_id']),
                                                                                how='inner', 
                                                                                left_on='read_id', 
                                                                                right_index = True,
                                                                                suffixes=['', '_y'],
                                                                                validate='m:1',
                                                                               )

    align_list_with_unique_best = align_list_with_max_score_count.query('max_score_count == 1').drop(['max_score_count'], axis=1)
    align_list_non_unique_best = align_list_with_max_score_count.query('max_score_count > 1').drop(['max_score_count'], axis=1)

    align_stat_with_unique_best = align_list_to_align_stat_by_assembly_id(
                                                                          assembly_metadata=assembly_metadata,
                                                                          log=log,
                                                                          align_list=align_list_with_unique_best,
                                                                          noise_bed=noise_bed
                                                                         )
    assembly_abundance = align_stat_with_unique_best[['assembly_id', 'adjusted_total_aligned_bp']].rename(columns={'adjusted_total_aligned_bp': 'assembly_abundance'})

    align_list_non_unique_best = align_list_non_unique_best.merge(
                                                                  right=assembly_abundance.set_index(['assembly_id']),
                                                                  how='left', 
                                                                  left_on='assembly_id', 
                                                                  right_index = True,
                                                                  suffixes=['', '_y'],
                                                                  validate='m:1',
                                                                 ).fillna(0)

    read_abundance = align_list_non_unique_best[['read_id', 'assembly_abundance']].groupby('read_id', as_index=False).sum().rename(columns={'assembly_abundance': 'read_abundance'})

    align_list_non_unique_best = align_list_non_unique_best.merge(
                                                                  right=read_abundance.set_index(['read_id']),
                                                                  how='left', 
                                                                  left_on='read_id', 
                                                                  right_index = True,
                                                                  suffixes=['', '_y'],
                                                                  validate='m:1',
                                                                 ).fillna(0)

    align_list_non_unique_best['alignment_score_tiebreaker'] = align_list_non_unique_best['alignment_score_tiebreaker'].apply(lambda x: random.random())
    align_list_non_unique_best = align_list_non_unique_best.assign(relative_abundance = lambda x: numpy.where(x['read_abundance'] <= 0, 1, x['assembly_abundance'] / x['read_abundance']))
    align_list_non_unique_best['alignment_score_tiebreaker'] = align_list_non_unique_best['alignment_score_tiebreaker'] * align_list_non_unique_best['relative_abundance']
    align_list_non_unique_best = align_list_non_unique_best.drop(['assembly_abundance', 'read_abundance', 'relative_abundance'], axis=1)

    combined_align_list = pandas.concat([align_list_with_unique_best, align_list_non_unique_best], axis=0)

    return combined_align_list.sort_values(['read_id', 'alignment_score', 'alignment_score_tiebreaker']).drop_duplicates(subset=['read_id'], keep='last').copy()


def align_list_to_bed(*,
                      align_list):

    temp_align_list = align_list.assign(assembly_id_sequence_id = lambda x: x['assembly_id'] + ',' + x['sequence_id'])
    temp_bed = BedTool.from_dataframe(temp_align_list[['assembly_id_sequence_id', 'sequence_from', 'sequence_to']])

    temp_merged_bed = temp_bed.sort().merge()
    if temp_merged_bed.count() > 0:
        temp_merged_bed_df = temp_merged_bed.to_dataframe()
        temp_merged_bed_df = pandas.concat([temp_merged_bed_df['chrom'].str.split(',', n=1, expand=True).rename(columns={0:'assembly_id', 1:'sequence_id'}), temp_merged_bed_df[['start', 'end']]], axis=1)
        bed = BedTool.from_dataframe(temp_merged_bed_df[['sequence_id', 'start', 'end', 'assembly_id']])
        os.remove(temp_merged_bed.fn)
    else:
        bed = BedTool('', from_string=True)

    os.remove(temp_bed.fn)

    return bed


def bed_to_covered_bp_by_assembly_id(*,
                                     bed):

    if bed.count() > 0:
        bed_df = bed.to_dataframe()
        if 'name' not in bed_df.columns:
            raise RuntimeError('Assembly_id not found in bed')
        else:
            covered_bp = bed_df.rename(columns={'name': 'assembly_id'}).assign(covered_bp = lambda x: x['end'] - x['start'])[['assembly_id', 'covered_bp']].groupby(['assembly_id'], as_index=False).sum().copy()
        return covered_bp
    else:
        covered_bp = pandas.DataFrame(columns=['assembly_id', 'covered_bp'])
        covered_bp['covered_bp'] = covered_bp['covered_bp'].astype('int')
        return covered_bp


def bed_to_covered_bp_by_sequence_id(*,
                                     bed):

    if bed.count() > 0:
        bed_df = bed.to_dataframe()
        covered_bp = bed_df.rename(columns={'chrom': 'sequence_id'}).assign(covered_bp = lambda x: x['end'] - x['start'])[['sequence_id', 'covered_bp']].groupby(['sequence_id'], as_index=False).sum().copy()
        return covered_bp
    else:
        covered_bp = pandas.DataFrame(columns=['sequence_id', 'covered_bp'])
        covered_bp['covered_bp'] = covered_bp['covered_bp'].astype('int')
        return covered_bp


def merge_bed_with_assembly_id(*,
                               bed_list):

    merged_bed = BedTool('', from_string=True)

    for bed in bed_list:
        if bed.count() > 0:
            bed_df = bed.to_dataframe().rename(columns={'chrom':'sequence_id', 'name':'assembly_id'})
            if 'assembly_id' not in bed_df.columns:
                bed_df = bed_df.assign(assembly_id = lambda x: '')
            bed_df = bed_df.assign(assembly_id_sequence_id = lambda x: x['assembly_id'] + ',' + x['sequence_id'])
            temp_bed = BedTool.from_dataframe(bed_df[['assembly_id_sequence_id', 'start', 'end']])
            merged_bed = merged_bed.cat(temp_bed)

    if merged_bed.count() == 0:
        return merged_bed
    else:
        merged_bed_df = merged_bed.sort().merge().to_dataframe()
        merged_bed_df = pandas.concat([merged_bed_df['chrom'].str.split(',', n=1, expand=True).rename(columns={0:'assembly_id', 1:'sequence_id'}), merged_bed_df[['start', 'end']]], axis=1)

        return BedTool.from_dataframe(merged_bed_df[['sequence_id', 'start', 'end', 'assembly_id']])


def intersect_bed_with_assembly_id(*,
                                   bed_list):

    if not bed_list:
        return BedTool('', from_string=True)

    intersect_bed = None

    for bed in bed_list:
        if bed.count() > 0:
            bed_df = bed.to_dataframe().rename(columns={'chrom':'sequence_id', 'name':'assembly_id'})
            if 'assembly_id' not in bed_df.columns:
                bed_df = bed_df.assign(assembly_id = lambda x: '')
            bed_df = bed_df.assign(assembly_id_sequence_id = lambda x: x['assembly_id'] + ',' + x['sequence_id'])
            temp_bed = BedTool.from_dataframe(bed_df[['assembly_id_sequence_id', 'start', 'end']]).sort()
            if intersect_bed == None:
                intersect_bed = temp_bed
            else:
                intersect_bed = intersect_bed.intersect(temp_bed).sort()

    if intersect_bed is None:
        return BedTool('', from_string=True)
    else:
        if intersect_bed.count() == 0:
            return intersect_bed
        else:
            intersect_bed_df = intersect_bed.to_dataframe()
            intersect_bed_df = pandas.concat([intersect_bed_df['chrom'].str.split(',', n=1, expand=True).rename(columns={0:'assembly_id', 1:'sequence_id'}), intersect_bed_df[['start', 'end']]], axis=1)

            return BedTool.from_dataframe(intersect_bed_df[['sequence_id', 'start', 'end', 'assembly_id']])


def align_list_to_depth_bed(*,
                            temp_dir_name,
                            align_list,
                            min_depth=None,
                            can_equal_to_min=True,
                            max_depth=None,
                            can_equal_to_max=True):
 
    bed_genome_filename = os.path.join(temp_dir_name, 'align_to_depth_bed.bed_genome')

    bed_genome = align_list.assign(assembly_id_sequence_id = lambda x: x['assembly_id'] + ',' + x['sequence_id']).sort_values(['assembly_id_sequence_id']).drop_duplicates(subset=['assembly_id_sequence_id'])
    bed_genome.to_csv(path_or_buf=bed_genome_filename, sep='\t', header=False, index=False, columns=['assembly_id_sequence_id', 'sequence_length'])

    temp_bed = BedTool.from_dataframe(align_list.assign(assembly_id_sequence_id = lambda x: x['assembly_id'] + ',' + x['sequence_id'])[['assembly_id_sequence_id', 'sequence_from', 'sequence_to']])
    temp_genome_coverage = temp_bed.sort().genome_coverage(bg=True, g=bed_genome_filename)

    genome_coverage = temp_genome_coverage.to_dataframe().rename(columns={'name': 'depth'})
    genome_coverage = pandas.concat([genome_coverage['chrom'].str.split(',', n=1, expand=True).rename(columns={0:'assembly_id', 1:'sequence_id'}), genome_coverage], axis=1)

    if max_depth is not None:
        genome_coverage = genome_coverage.merge(
                                                right=max_depth.set_index(['assembly_id']),
                                                how='inner', 
                                                left_on='assembly_id', 
                                                right_index = True,
                                                suffixes=['', '_y'],
                                                validate='m:1',
                                               ).fillna(99999999)

        query_string = 'depth <'
        if can_equal_to_max == True:
            query_string = query_string + '='
        query_string = query_string + ' max_depth'

    if min_depth is not None:
        genome_coverage = genome_coverage.merge(
                                                right=min_depth.set_index(['assembly_id']),
                                                how='inner', 
                                                left_on='assembly_id', 
                                                right_index = True,
                                                suffixes=['', '_y'],
                                                validate='m:1',
                                               ).fillna(-1)

        if max_depth is not None:
            query_string = query_string + ' and '
        else:
            query_string = ''
        query_string = query_string + 'depth >'
        if can_equal_to_min == True:
            query_string = query_string + '='
        query_string = query_string + ' min_depth'

    if max_depth is None and min_depth is None:
        depth_genome_coverage = genome_coverage
    else:
        depth_genome_coverage = genome_coverage.query(query_string)

    depth_span_bp = depth_genome_coverage.assign(span_bp = lambda x: x['end'] - x['start'])[['assembly_id', 'span_bp']].groupby(['assembly_id'], as_index=False).sum().copy()
    depth_bed = merge_bed_with_assembly_id(bed_list=[BedTool.from_dataframe(depth_genome_coverage[['sequence_id', 'start', 'end', 'assembly_id']])])

    os.remove(temp_bed.fn)
    os.remove(temp_genome_coverage.fn)
    os.remove(bed_genome_filename)

    return depth_bed, depth_span_bp
   

def summary_stat_1(*,
	              align_list,
	              groupby_col):
    align_list_plus_summary = align_list.assign(
                                                total_number_of_read = lambda x: 1,
                                                total_aligned_bp = lambda x: x['sequence_to'] - x['sequence_from'],
                                               ).rename(columns={'read_length': 'total_read_bp'})[[groupby_col, 'total_number_of_read', 'total_read_bp', 'total_aligned_bp', 'match', 'edit_dist', 'alignment_score', 'alignment_score_tiebreaker']]

    return align_list_plus_summary.groupby([groupby_col], as_index=False).sum().copy()

def summary_stat_2(*,
	              align_stat,
	              length_col):
    align_stat_plus_summary = align_stat.assign(
                            	                  average_read_length = lambda x: x['total_read_bp'] / x['total_number_of_read'],
                                                average_depth = lambda x: x['total_aligned_bp'] / x[length_col],
                                                covered_percent = lambda x: x['covered_bp'] / x[length_col],
                                                noise_span_percent = lambda x: x['noise_span_bp'] / x[length_col],
                                                adjusted_covered_percent = lambda x: x['covered_bp'] / (x[length_col] - x['noise_span_bp']),
                                                average_identity = lambda x: x['match'] / x['total_aligned_bp'],
                                                average_edit_dist = lambda x: x['edit_dist'] / x['total_aligned_bp'],
                                                average_alignment_score = lambda x: x['alignment_score'] / x['total_aligned_bp'],
                                               ).replace([numpy.inf, -numpy.inf], numpy.nan).fillna(0)
    align_stat_plus_summary = align_stat_plus_summary.assign(
                                                             log_uncovered_percent = lambda x: numpy.where(1 - x['adjusted_covered_percent'] > 0, numpy.log(1 - x['adjusted_covered_percent']), 0),
                                                            ).assign(
                                                                     adjusted_average_depth = lambda x: x['adjusted_covered_percent'] * x['total_aligned_bp'] / (x[length_col] - x['noise_span_bp']),
                                                                    ).drop(['log_uncovered_percent'], axis=1)
    align_stat_plus_summary = align_stat_plus_summary.assign(
                                                             adjusted_total_aligned_bp = lambda x: x['adjusted_average_depth'] * x[length_col],
                                                            ).replace([numpy.inf, -numpy.inf], numpy.nan).fillna(0)
    align_stat_plus_summary[['adjusted_total_aligned_bp']] = align_stat_plus_summary[['adjusted_total_aligned_bp']].round(0).astype('int')

    return align_stat_plus_summary


def align_list_to_align_stat_by_assembly_id(*,
                                            assembly_metadata,
                                            log,
                                            align_list,
                                            noise_bed=None):

    assembly_best_align_list = align_list.sort_values(['read_id', 'assembly_id', 'alignment_score', 'alignment_score_tiebreaker']).drop_duplicates(subset=['read_id', 'assembly_id'], keep='last')

    align_list_grouped = summary_stat_1(
    	                                  align_list=assembly_best_align_list,
    	                                  groupby_col='assembly_id',
    	                                 )

    align_stat = align_list_grouped

    align_stat = assembly_metadata.get_assembly_length(assembly_list=align_stat, how='left').fillna(0)
    align_stat['assembly_length'] = align_stat['assembly_length'].astype('int')

    align_stat = assembly_metadata.get_tax_id(assembly_list=align_stat, how='left').fillna(0)
    align_stat['tax_id'] = align_stat['tax_id'].astype('int')
    align_stat['species_tax_id'] = align_stat['species_tax_id'].astype('int')
    align_stat['genus_tax_id'] = align_stat['genus_tax_id'].astype('int')
    align_stat['genus_height'] = align_stat['genus_height'].astype('int')

    covered_bed = align_list_to_bed(align_list=assembly_best_align_list)
    if noise_bed != None:
        covered_bed = covered_bed.subtract(noise_bed)

    covered_bp = bed_to_covered_bp_by_assembly_id(bed=covered_bed)

    align_stat = align_stat.merge(
                                  right=covered_bp.set_index(['assembly_id']),
                                  how='left', 
                                  left_on='assembly_id', 
                                  right_index = True,
                                  suffixes=['', '_y'],
                                  validate='1:1',
                                 ).fillna(0)
    align_stat['covered_bp'] = align_stat['covered_bp'].astype(int)

    if noise_bed is None:
        align_stat = align_stat.assign(noise_span_bp = lambda x : 0)
    else:
        noise_span_bp = bed_to_covered_bp_by_assembly_id(bed=noise_bed).rename(columns={'covered_bp':'noise_span_bp'})
        align_stat = align_stat.merge(
                                  right=noise_span_bp.set_index(['assembly_id']),
                                  how='left', 
                                  left_on='assembly_id', 
                                  right_index = True,
                                  suffixes=['', '_y'],
                                  validate='m:1',
                                 ).fillna(0)
        align_stat['noise_span_bp'] = align_stat['noise_span_bp'].astype('int')

    if align_stat.shape[0] != align_list_grouped.shape[0]:
        log.print_missing_data(align_list_grouped, align_stat, 'assembly_id', 'assembly_length or assembly_id_tax_id')


    return summary_stat_2(
      	                  align_stat=align_stat,
      	                  length_col='assembly_length',
      	                 )


def align_list_to_align_stat_by_sequence_id(*,
                                            assembly_metadata,
                                            log,
                                            align_list,
                                            noise_bed=None):

    sequence_best_align_list = align_list.sort_values(['read_id', 'sequence_id', 'alignment_score', 'alignment_score_tiebreaker']).drop_duplicates(subset=['read_id', 'sequence_id'], keep='last')

    align_list_grouped = summary_stat_1(
    	                                  align_list=sequence_best_align_list.assign(assembly_id_sequence_id = lambda x: x['assembly_id'] + ',' + x['sequence_id']),
    	                                  groupby_col='assembly_id_sequence_id',
    	                                 )
    align_list_grouped = pandas.concat([align_list_grouped['assembly_id_sequence_id'].str.split(',', n=1, expand=True).rename(columns={0:'assembly_id', 1:'sequence_id'}), align_list_grouped.drop(['assembly_id_sequence_id'], axis=1)], axis=1)
    align_stat = align_list_grouped

    align_stat = assembly_metadata.get_sequence_length(sequence_list=align_stat, how='left').fillna(0)
    align_stat['sequence_length'] = align_stat['sequence_length'].astype('int')

    covered_bed = align_list_to_bed(align_list=sequence_best_align_list)
    if noise_bed != None:
        covered_bed = covered_bed.subtract(noise_bed)

    covered_bp = bed_to_covered_bp_by_sequence_id(bed=covered_bed)

    align_stat = align_stat.merge(
                                  right=covered_bp.set_index(['sequence_id']),
                                  how='left', 
                                  left_on='sequence_id', 
                                  right_index = True,
                                  suffixes=['', '_y'],
                                  validate='1:1',
                                 ).fillna(0)
    align_stat['covered_bp'] = align_stat['covered_bp'].astype(int)

    if noise_bed is None:
        align_stat = align_stat.assign(noise_span_bp = lambda x : 0)
    else:
        noise_span_bp = bed_to_covered_bp_by_sequence_id(bed=noise_bed).rename(columns={'covered_bp':'noise_span_bp'})
        align_stat = align_stat.merge(
                                  right=noise_span_bp.set_index(['sequence_id']),
                                  how='left', 
                                  left_on='sequence_id', 
                                  right_index = True,
                                  suffixes=['', '_y'],
                                  validate='m:1',
                                 ).fillna(0)
        align_stat['noise_span_bp'] = align_stat['noise_span_bp'].astype('int')

    if align_stat.shape[0] != align_list_grouped.shape[0]:
        log.print_missing_data(align_list_grouped, align_stat, 'sequence_id', 'sequence_length')

    return summary_stat_2(
    	                 align_stat=align_stat,
    	                 length_col='sequence_length'
    	                )


def good_align_list(*,
                    align_list,
                    good_align_threshold):

    assembly_best_align_list = align_list.sort_values(['read_id', 'assembly_id', 'alignment_score', 'alignment_score_tiebreaker']).drop_duplicates(subset=['read_id', 'assembly_id'], keep='last')

    best_align_score = assembly_best_align_list[['read_id', 'alignment_score']].sort_values(['read_id', 'alignment_score']).drop_duplicates(subset=['read_id'], keep='last').rename(columns={'alignment_score':'best_alignment_score'}).set_index('read_id')

    align_list_with_best_align_score = assembly_best_align_list.merge(
                                                                      right=best_align_score,
                                                                      how='inner', 
                                                                      left_on='read_id', 
                                                                      right_index = True,
                                                                      suffixes=['', '_y'],
                                                                      validate='m:1',
                                                                     )
    if align_list_with_best_align_score.shape[0] != assembly_best_align_list.shape[0]:
        raise RuntimeError('Cannot join align_list')

    good_align_threshold_percent = good_align_threshold / 100

    return align_list_with_best_align_score.query('alignment_score >= best_alignment_score * @good_align_threshold_percent').drop(['best_alignment_score'], axis=1)


def select_alignment_by_bed(*,
                            temp_dir_name,
                            align_list,
                            bed,
                            max_overlap=100,
                            can_equal_to_max=True,
                            min_overlap=0,
                            can_equal_to_min=True):

    if bed is None or bed.count() == 0:
        if (min_overlap == 0 and can_equal_to_min == True) and (max_overlap != 0 or can_equal_to_max == True):
            return align_list
        else:
            empty_align_list = pandas.DataFrame(columns=align_list_col_name)
            for key, value in align_list_col_type.items():
                empty_align_list[key] = empty_align_list[key].astype(value)
            return empty_align_list

    align_bed_filename = os.path.join(temp_dir_name, 'select_alignment_by_bed.align.bed')

    indexed_align_list = align_list.copy().reset_index()
    indexed_align_list['align_index'] = indexed_align_list.index
    indexed_align_list = indexed_align_list.assign(assembly_id_sequence_id = lambda x: x['assembly_id'] + ',' + x['sequence_id'])

    align_bed = BedTool.from_dataframe(indexed_align_list[['assembly_id_sequence_id', 'sequence_from', 'sequence_to', 'read_id', 'align_index', 'strand']])

    bed_df = bed.to_dataframe()
    temp_bed = BedTool.from_dataframe(bed_df.assign(assembly_id_sequence_id = lambda x: x['name'] + ',' + x['chrom'])[['assembly_id_sequence_id', 'start', 'end']])

    annotated_align_bed = align_bed.annotate(files=temp_bed.fn)

    annotated_align_list = indexed_align_list.merge(
                                                    right=annotated_align_bed.to_dataframe()[['chrom', 'start', 'end', 'name', 'score', 'thickStart']].set_index(['score']),
                                                    how='inner', 
                                                    left_on='align_index', 
                                                    right_index = True,
                                                    suffixes=['', '_y'],
                                                    validate='1:1',
                                                   ).drop(['align_index'], axis=1)
    if annotated_align_list.shape[0] != align_list.shape[0]:
        raise LookupError('Some alignment cannot be matched')

    max_overlap_percent = max_overlap / 100
    min_overlap_percent = min_overlap / 100
    query_string = 'thickStart <' + ('=' if can_equal_to_max == True else '') + ' ' + str(max_overlap_percent) + ' and thickStart >' + ('=' if can_equal_to_min == True else '') + ' ' + str(min_overlap_percent)
    filtered_align_list = annotated_align_list.query(query_string).drop(['assembly_id_sequence_id', 'chrom', 'start', 'end', 'name', 'thickStart'], axis=1).copy()

    os.remove(align_bed.fn)
    os.remove(temp_bed.fn)
    os.remove(annotated_align_bed.fn)

    return filtered_align_list


def read_genome_set(*,
                    global_options,
                    genome_set_name):

    genome_set_filename = os.path.join(global_options['config_folder'], genome_set_name)

    if os.path.exists(genome_set_filename) == False or os.path.isfile(genome_set_filename) == False:
        raise RuntimeError('Genome_set ' + genome_set_name + ' not found')

    return pandas.read_csv(genome_set_filename, sep='\t', names=('assembly_id',), usecols=(0,), index_col=False, header=None)


def get_human_noise_bed(*,
                        assembly_metadata,
                        assembly_list,
                        assembly_folder,
                        log):

    merged_bed = BedTool('', from_string=True)

    human_noise_bed_path = assembly_metadata.get_human_noise_bed_path(assembly_list=assembly_list)

    human_noise_bed_path['path'] = human_noise_bed_path['path'].map(lambda x: os.path.join(assembly_folder, x))

    for bed_path in human_noise_bed_path.itertuples():
        if os.path.exists(bed_path.path) == False or os.path.isfile(bed_path.path) == False:
	        log.print('Human noise bed ' + bed_path.path + ' not found')
        else:
            human_noise_bed = BedTool(bed_path.path)
            human_noise_bed_df = human_noise_bed.to_dataframe()
            human_noise_bed_df['chrom'] = human_noise_bed_df['chrom'].map(lambda x: bed_path.assembly_id + ',' + x)
            temp_bed = BedTool.from_dataframe(human_noise_bed_df)
            merged_bed = merged_bed.cat(temp_bed)
            os.remove(temp_bed.fn)

    merged_bed = merged_bed.sort().merge()

    if merged_bed.count() > 0:
        merged_bed_df = merged_bed.to_dataframe()
        merged_bed_df = pandas.concat([merged_bed_df['chrom'].str.split(',', n=1, expand=True).rename(columns={0:'assembly_id', 1:'sequence_id'}), merged_bed_df[['start', 'end']]], axis=1)
        os.remove(merged_bed.fn)

        noise_bed = BedTool.from_dataframe(merged_bed_df[['sequence_id', 'start', 'end', 'assembly_id']])

        covered_bp = bed_to_covered_bp_by_assembly_id(bed=noise_bed)

        noise_span_bp = assembly_list[['assembly_id']].merge(
                                                             right=covered_bp.rename(columns={'covered_bp':'noise_span_bp'}).set_index(['assembly_id']),
                                                             how='left', 
                                                             left_on='assembly_id',
                                                             right_index = True,
                                                             suffixes=['', '_y'],
                                                             validate='1:1',
                                                            ).fillna(0).copy()
        noise_span_bp['noise_span_bp'] = noise_span_bp['noise_span_bp'].astype('int')

    else:
        noise_bed = merged_bed
        noise_span_bp = assembly_list[['assembly_id']].assign(noise_span_bp = lambda x: 0).copy()

    return noise_bed, noise_span_bp


def read_db(*,
            db_conn,
            sql,
            input_key,
            output_record=None):

    sqlite_max_variable_number = 999  # cannot call sqlite3.Connection.get_limit() so it is hard-coded here

    num_iteration = (input_key.shape[0] + sqlite_max_variable_number - 1) // sqlite_max_variable_number
    
    for i in range(num_iteration):
        num_key = min(sqlite_max_variable_number, input_key.shape[0] - i * sqlite_max_variable_number)
        key_list = tuple(input_key['key'].iloc[i * sqlite_max_variable_number: i * sqlite_max_variable_number + num_key])
        sql_formatted = sql.format(key=','.join(['?'] * num_key))
        output_record = pandas.concat([output_record, pandas.read_sql(sql_formatted, db_conn, params=key_list)], axis=0)

    return output_record


def get_assembly_info(*,
                      assembly_metadata,
                      db_conn,
                      log,
                      assembly_list):

    sql = 'SELECT'
    for col_name in assembly_info_col_name:
        sql = sql + ' ' + col_name + ','
    sql = sql[:-1]

    sql = sql + " FROM assembly_summary WHERE assembly_id in ({key})"

    output_assembly_list = pandas.DataFrame(columns=assembly_info_col_name)
    for key, value in assembly_info_col_type.items():
        output_assembly_list[key] = output_assembly_list[key].astype(value)

    output_assembly_list = read_db(
                                   db_conn=db_conn,
                                   sql=sql,
                                   input_key=assembly_list[['assembly_id']].rename(columns={'assembly_id': 'key'}),
                                   output_record=output_assembly_list)

    if output_assembly_list.shape[0] != assembly_list.shape[0]:
        log.print_missing_data(assembly_list, output_assembly_list, 'assembly_id', 'assembly_summary')

    output_assembly_list = assembly_metadata.get_tax_id(assembly_list=output_assembly_list, how='left').fillna(0)
    output_assembly_list['tax_id'] = output_assembly_list['tax_id'].astype('int')
    output_assembly_list['species_tax_id'] = output_assembly_list['species_tax_id'].astype('int')
    output_assembly_list['genus_tax_id'] = output_assembly_list['genus_tax_id'].astype('int')

    tax_id_list = pandas.concat([
                                 output_assembly_list[['tax_id']],
                                 output_assembly_list[['species_tax_id']].rename(columns={'species_tax_id':'tax_id'}),
                                 output_assembly_list[['genus_tax_id']].rename(columns={'genus_tax_id':'tax_id'}),
                                ], axis=0).sort_values(['tax_id']).drop_duplicates().query('tax_id != 0')

    tax_id_list = get_tax_name(
                               db_conn=db_conn, 
                               log=log,
                               tax_id_list=tax_id_list,
                              ).set_index('tax_id')
    
    output_assembly_list = output_assembly_list.merge(
                                                      right=tax_id_list,
                                                      how='left', 
                                                      left_on='tax_id', 
                                                      right_index = True,
                                                      suffixes=['', '_y'],
                                                      validate='m:1',
                                                     )
    output_assembly_list = output_assembly_list.merge(
                                                      right=tax_id_list.rename(columns={'tax_name':'species_tax_name'}),
                                                      how='left', 
                                                      left_on='species_tax_id', 
                                                      right_index = True,
                                                      suffixes=['', '_y'],
                                                      validate='m:1',
                                                     )
    output_assembly_list = output_assembly_list.merge(
                                                      right=tax_id_list.rename(columns={'tax_name':'genus_tax_name'}),
                                                      how='left', 
                                                      left_on='genus_tax_id', 
                                                      right_index = True,
                                                      suffixes=['', '_y'],
                                                      validate='m:1',
                                                     )

    return output_assembly_list


def get_tax_name(*,
                 db_conn,
                 log,
                 tax_id_list):

    sql = "SELECT tax_id, tax_name FROM names WHERE tax_id in({key}) and is_primary = 1"
    
    output_tax_id_list = pandas.DataFrame(columns=['tax_id', 'tax_name'])
    output_tax_id_list['tax_id'] = output_tax_id_list['tax_id'].astype('int')

    output_tax_id_list = read_db(
                                 db_conn=db_conn,
                                 sql=sql,
                                 input_key=tax_id_list[['tax_id']].rename(columns={'tax_id': 'key'}),
                                 output_record=output_tax_id_list)
    
    output_tax_id_list['tax_id'] = output_tax_id_list['tax_id'].astype('int')

    if output_tax_id_list.shape[0] != tax_id_list.shape[0]:
        log.print_missing_data(tax_id_list, output_tax_id_list, 'tax_id', 'tax_name')

    return output_tax_id_list


def get_sequence_name(*,
                      db_conn,
                      log,
                      sequence_list):
    
    output_sequence_list = pandas.DataFrame(columns=['sequence_id', 'sequence_name'])

    sql = "SELECT sequence_id, sequence_name FROM sequence_name WHERE sequence_id in({key})"

    output_sequence_list = read_db(
                                   db_conn=db_conn,
                                   sql=sql,
                                   input_key=sequence_list[['sequence_id']].rename(columns={'sequence_id': 'key'}),
                                   output_record=pandas.DataFrame(columns=['sequence_id', 'sequence_name']))

    if output_sequence_list.shape[0] != sequence_list.shape[0]:
        log.print_missing_data(sequence_list, output_sequence_list, 'sequence_id', 'sequence_name')

    return output_sequence_list


def generate_align_stat_group_by_assembly_id(*,
                                             assembly_metadata,
                                             log,
                                             align_list,
                                             assembly_info,
                                             noise_bed=None):

    if align_list is None or align_list.empty == True:
        align_stat = pandas.DataFrame(columns=align_stat_col_name_by_assembly_id)
        for key, value in align_stat_col_type_by_assembly_id.items():
            align_stat[key] = align_stat[key].astype(value)

        return align_stat

    align_stat = align_list_to_align_stat_by_assembly_id(
                                                         assembly_metadata=assembly_metadata,
                                                         log=log,
                                                         align_list=align_list,
                                                         noise_bed=noise_bed,
                                                        )

    align_stat_with_assembly_info = align_stat.merge(
                                                     right=assembly_info[['assembly_id', 'tax_name', 'species_tax_name', 'genus_tax_name']].set_index(['assembly_id']),
                                                     how='left', 
                                                     left_on='assembly_id', 
                                                     right_index = True,
                                                     suffixes=['', '_y'],
                                                     validate='m:1',
                                                    )

    return align_stat_with_assembly_info[align_stat_col_name_by_assembly_id].copy()


def generate_align_stat_group_by_sequence_id(*,
                                             assembly_metadata,
                                             log,
                                             align_list,
                                             sequence_info,
                                             noise_bed=None):

    if align_list is None or align_list.empty == True:
        align_stat = pandas.DataFrame(columns=align_stat_raw_col_name_by_sequence_id)
        for key, value in align_stat_raw_col_type_by_sequence_id.items():
            align_stat[key] = align_stat[key].astype(value)

        return align_stat

    align_stat = align_list_to_align_stat_by_sequence_id(
                                                         assembly_metadata=assembly_metadata,
                                                         log=log,
                                                         align_list=align_list,
                                                         noise_bed=noise_bed,
                                                        )
    sequence_info.to_csv('sequence_info',sep='\t')
    align_stat_with_assembly_info = align_stat.merge(
                                                     right=sequence_info[['sequence_id', 'sequence_name']].set_index(['sequence_id']),
                                                     how='left', 
                                                     left_on='sequence_id', 
                                                     right_index = True,
                                                     suffixes=['', '_y'],
                                                     validate='m:1',
                                                    )

    return align_stat_with_assembly_info[align_stat_raw_col_name_by_sequence_id].copy()



def step_adaptor_trimming(megapath_nano, adaptor_trimming):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    adaptor_trimming.O.query_filename_list = adaptor_trimming.I.query_filename_list.copy()
    adaptor_trimming.O.query_filename_list['path'] = adaptor_trimming.O.query_filename_list['path'].map(lambda x: os.path.join(megapath_nano.temp_dir_name, os.path.split(x)[1] + '.adaptor_trimmed'))

    trimming_process_list = []

    for file_index, filename in enumerate(adaptor_trimming.I.query_filename_list['path']):
        trimming_command = ['porechop', megapath_nano.global_options['porechopThreadOption'], '--require_two_barcodes','-i', filename, '-o', adaptor_trimming.O.query_filename_list['path'][file_index]]
        trimming_process = subprocess.Popen(trimming_command, close_fds=True)
        trimming_process_list.append(trimming_process)

    for index, trimming_process in enumerate(trimming_process_list):
        trimming_process.wait()
        if trimming_process.returncode != 0:
	          raise RuntimeError('Error encountered in adaptor_trimming', adaptor_trimming.I.query_filename_list['path'][index])

    megapath_nano.log.print('end')


def step_read_trimming_and_filter(megapath_nano, read_trimming_and_filter):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()
    megapath_nano.log.print('Read trimming:')
    megapath_nano.log.print('- head:  {crop} '.format(crop=read_trimming_and_filter.I.head_crop))
    megapath_nano.log.print('- tail:  {crop} '.format(crop=read_trimming_and_filter.I.tail_crop))
    megapath_nano.log.print('Filtering criteria:')
    megapath_nano.log.print('- read length >= {min_read_length} '.format(min_read_length=read_trimming_and_filter.I.min_read_length))
    megapath_nano.log.print('- average quality >= {min_quality} '.format(min_quality=read_trimming_and_filter.I.min_read_quality))

    RAM_dir_name = tempfile.mkdtemp(prefix='read_trimming_and_filter.', dir=megapath_nano.RAM_dir_name)

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    # Intermediate files between steps are placed in megapath_nano.temp_dir_name
    read_trimming_and_filter.O.query_filename_list = read_trimming_and_filter.I.query_filename_list.copy()
    read_trimming_and_filter.O.query_filename_list['path'] = read_trimming_and_filter.O.query_filename_list['path'].map(lambda x: os.path.join(megapath_nano.temp_dir_name, os.path.split(x)[1] + '.trimmed_and_filtered'))

    read_trimming_and_filter.O.read_info = None
    read_info_filename_list = read_trimming_and_filter.I.query_filename_list.copy()
    read_info_filename_list['path'] = read_info_filename_list['path'].map(lambda x: os.path.join(RAM_dir_name, os.path.split(x)[1] + '.read_info'))

    filter_process_list = []
    output_query_file_list = []
    read_info_file_list = []
    input_query_file_list = []

    for file_index, filename in enumerate(read_trimming_and_filter.I.query_filename_list['path']):

        output_query_file = os.open(read_trimming_and_filter.O.query_filename_list['path'][file_index], flags=os.O_CREAT | os.O_EXCL | os.O_WRONLY, mode=0o644)
        read_info_file = os.open(read_info_filename_list['path'][file_index], flags=os.O_CREAT | os.O_EXCL | os.O_WRONLY, mode=0o644)

        input_query_file = os.open(filename, os.O_RDONLY)

        filter_command = [os.path.join(megapath_nano.global_options['tool_folder'], 'nanofastq'),]
        if read_trimming_and_filter.I.min_read_length is not None:
            filter_command.extend(['-l', str(read_trimming_and_filter.I.min_read_length)])
        if read_trimming_and_filter.I.min_read_quality is not None:
            filter_command.extend(['-q', str(read_trimming_and_filter.I.min_read_quality)])
        if read_trimming_and_filter.I.head_crop is not None:
            filter_command.extend(['-h', str(read_trimming_and_filter.I.head_crop)])
        if read_trimming_and_filter.I.head_crop is not None:
            filter_command.extend(['-t', str(read_trimming_and_filter.I.tail_crop)])
        if read_trimming_and_filter.I.reassign_read_id == True:
            filter_command.extend(['-r', os.path.split(filename)[1]])

        filter_process = subprocess.Popen(filter_command, close_fds=True, stdin=input_query_file, stdout=output_query_file, stderr=read_info_file)
  
        filter_process_list.append(filter_process)
        output_query_file_list.append(output_query_file)
        read_info_file_list.append(read_info_file)
        input_query_file_list.append(input_query_file)


    for index, filter_process in enumerate(filter_process_list):

        filter_process.wait()
        if filter_process.returncode != 0:
            os.sys.exit('Error encountered in read_trimming_and_filter:', input_query_file_list[index])

        os.close(output_query_file_list[index])
        os.close(read_info_file_list[index])
        os.close(input_query_file_list[index])

        read_trimming_and_filter.O.read_info = pandas.concat([read_trimming_and_filter.O.read_info, pandas.read_csv(read_info_filename_list['path'][index], sep='\t', names=(read_info_col_name), usecols=(read_info_col_name), index_col=False, header=None)], axis=0)

    read_trimming_and_filter.O.read_info=read_trimming_and_filter.O.read_info.drop_duplicates(subset="read_id")
    read_trimming_and_filter.O.num_read = read_trimming_and_filter.O.read_info.shape[0]
    read_trimming_and_filter.O.passed_read_id_list = read_trimming_and_filter.O.read_info.query('passed_filter == 1')[['read_id', 'read_length']].copy()
    read_trimming_and_filter.O.passed_read_id_list=read_trimming_and_filter.O.passed_read_id_list.drop_duplicates()
    read_trimming_and_filter.O.num_read_passed_filter = read_trimming_and_filter.O.passed_read_id_list.shape[0]

    if megapath_nano.global_options['debug'] == False:
        shutil.rmtree(RAM_dir_name)

    megapath_nano.log.print('Number of reads before filter: {num_read}'.format(num_read=read_trimming_and_filter.O.num_read))
    megapath_nano.log.print('Number of reads passed filter: {num_read}'.format(num_read=read_trimming_and_filter.O.num_read_passed_filter))
    megapath_nano.log.print('end')


def step_human_and_decoy_filter(megapath_nano, human_and_decoy_filter):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()
    megapath_nano.log.print('Filtering criteria for human:')
    megapath_nano.log.print('- alignment score >= {threshold} '.format(threshold=human_and_decoy_filter.I.human_min_alignment_score))
    megapath_nano.log.print('  or')
    megapath_nano.log.print('- alignment score / read length >= {threshold:.2%} '.format(threshold=human_and_decoy_filter.I.human_min_alignment_score_percent/100))
    megapath_nano.log.print('Filtering criteria for decoy:')
    megapath_nano.log.print('- alignment score >= {threshold} '.format(threshold=human_and_decoy_filter.I.decoy_min_alignment_score))
    megapath_nano.log.print('  or')
    megapath_nano.log.print('- alignment score / read length >= {threshold:.2%} '.format(threshold=human_and_decoy_filter.I.decoy_min_alignment_score_percent/100))

    temp_dir_name = tempfile.mkdtemp(prefix='human_and_decoy_filter.', dir=megapath_nano.temp_dir_name)
    RAM_dir_name = tempfile.mkdtemp(prefix='human_and_decoy_filter.', dir=megapath_nano.RAM_dir_name)

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    human_and_decoy_filter.O.query_filename_list = human_and_decoy_filter.I.query_filename_list.copy()

    target_assembly_list = pandas.concat([human_and_decoy_filter.I.human_assembly_list, human_and_decoy_filter.I.decoy_assembly_list], axis=0)

    align_list = Align(
                       assembly_metadata=megapath_nano.assembly_metadata,
                       global_options=megapath_nano.global_options,
                       temp_dir_name=RAM_dir_name,
                       log_file=megapath_nano.aligner_log,
                       query_filename_list=human_and_decoy_filter.I.query_filename_list,
                       target_assembly_list=target_assembly_list,
                       aligner_options=shlex.split(megapath_nano.global_options['alignerThreadOption'] + ' -x map-ont'),
                       mapping_only=megapath_nano.global_options['mapping_only'],
                    )

    # human filtering first, then decoy filtering 

    human_best_align_list = align_list.merge(
                                             right=human_and_decoy_filter.I.human_assembly_list.set_index('assembly_id'),
                                             how='inner', 
                                             left_on='assembly_id', 
                                             right_index = True,
                                             suffixes=['', '_y'],
                                             validate='m:1',
                                            ).sort_values(['read_id', 'alignment_score', 'alignment_score_tiebreaker']).drop_duplicates(subset=['read_id'], keep='last')

    human_and_decoy_filter.O.human_best_align_list = human_best_align_list.query('alignment_score >= @human_and_decoy_filter.I.human_min_alignment_score or alignment_score * 100 / read_length >= @human_and_decoy_filter.I.human_min_alignment_score_percent')

    num_human_read_with_alignment_score_above_threshold = human_and_decoy_filter.O.human_best_align_list.query('alignment_score >= @human_and_decoy_filter.I.human_min_alignment_score').shape[0]
    num_human_read_with_normalized_alignment_score_above_threshold = human_and_decoy_filter.O.human_best_align_list.query('alignment_score * 100 / read_length >= @human_and_decoy_filter.I.human_min_alignment_score_percent').shape[0]
    
    #if not covered by bed
    
    
    human_and_decoy_filter.O.human_read_id_list = human_and_decoy_filter.O.human_best_align_list[['read_id', 'read_length']].sort_values(['read_id', 'read_length']).drop_duplicates()

    remaining_align_list = align_list.merge(
                                            right=human_and_decoy_filter.O.human_read_id_list.set_index('read_id').rename(columns={'read_length': 'filtered'}),
                                            how='left', 
                                            left_on='read_id', 
                                            right_index = True,
                                            suffixes=['', '_y'],
                                            validate='m:1',
                                           ).fillna(0).query('filtered == 0').drop(['filtered'], axis=1)

    decoy_best_align_list = remaining_align_list.merge(
		                                               right=human_and_decoy_filter.I.decoy_assembly_list.set_index('assembly_id'),
		                                               how='inner', 
		                                               left_on='assembly_id', 
		                                               right_index = True,
		                                               suffixes=['', '_y'],
		                                               validate='m:1',
		                                              ).sort_values(['read_id', 'alignment_score', 'alignment_score_tiebreaker']).drop_duplicates(subset=['read_id'], keep='last')

    human_and_decoy_filter.O.decoy_best_align_list = decoy_best_align_list.query('alignment_score >= @human_and_decoy_filter.I.decoy_min_alignment_score or alignment_score * 100 / read_length >= @human_and_decoy_filter.I.decoy_min_alignment_score_percent')

    num_decoy_read_with_alignment_score_above_threshold = human_and_decoy_filter.O.decoy_best_align_list.query('alignment_score >= @human_and_decoy_filter.I.decoy_min_alignment_score').shape[0]
    num_decoy_read_with_normalized_alignment_score_above_threshold = human_and_decoy_filter.O.decoy_best_align_list.query('alignment_score * 100 / read_length >= @human_and_decoy_filter.I.decoy_min_alignment_score_percent').shape[0]
    
    human_and_decoy_filter.O.decoy_read_id_list = human_and_decoy_filter.O.decoy_best_align_list[['read_id', 'read_length']].sort_values(['read_id', 'read_length']).drop_duplicates()

    human_and_decoy_filter.O.microbe_best_align_list = remaining_align_list.merge(
								                                                right=human_and_decoy_filter.O.decoy_read_id_list.set_index('read_id').rename(columns={'read_length': 'filtered'}),
								                                                how='left', 
								                                                left_on='read_id', 
								                                                right_index = True,
								                                                suffixes=['', '_y'],
								                                                validate='m:1',
								                                               )
    human_and_decoy_filter.O.microbe_best_align_list = human_and_decoy_filter.O.microbe_best_align_list.fillna(0).query('filtered == 0').drop(['filtered'], axis=1).sort_values(['read_id', 'alignment_score', 'alignment_score_tiebreaker']).drop_duplicates(subset=['read_id'], keep='last')
    human_and_decoy_filter.I.read_id_list=human_and_decoy_filter.I.read_id_list.drop_duplicates()
    human_and_decoy_filter.O.microbe_read_id_list = human_and_decoy_filter.I.read_id_list.merge(
                                      					                                                right=pandas.concat([human_and_decoy_filter.O.human_read_id_list, human_and_decoy_filter.O.decoy_read_id_list], axis=0).sort_values(['read_id', 'read_length']).drop_duplicates().set_index('read_id').rename(columns={'read_length': 'filtered'}),
                                      					                                                how='left',
                                      					                                                left_on='read_id',
                                      					                                                right_index = True,
                                      					                                                suffixes=['', '_y'],
                                      					                                                validate='1:1',
                                      					                                               ).fillna(0).query('filtered == 0').drop(['filtered'], axis=1)[['read_id', 'read_length']]
    
    print("human_and_decoy_filter.O.human_read_id_list.shape[0]")
    print(human_and_decoy_filter.O.human_read_id_list.shape[0])
    
    print("human_and_decoy_filter.O.decoy_read_id_list.shape[0]")
    print(human_and_decoy_filter.O.decoy_read_id_list.shape[0])
    
    print("human_and_decoy_filter.O.microbe_read_id_list.shape[0]")
    print(human_and_decoy_filter.O.microbe_read_id_list.shape[0])
    
    print("human_and_decoy_filter.I.read_id_list.shape[0]")
    print(human_and_decoy_filter.I.read_id_list.shape[0])
    """
    if human_and_decoy_filter.O.human_read_id_list.shape[0] + human_and_decoy_filter.O.decoy_read_id_list.shape[0] + human_and_decoy_filter.O.microbe_read_id_list.shape[0] != human_and_decoy_filter.I.read_id_list.shape[0]:
        raise RuntimeError('Total number of reads not match')
    if human_and_decoy_filter.O.human_read_id_list.shape[0] != human_and_decoy_filter.O.human_best_align_list.shape[0]:
        raise RuntimeError('Human reads and align list not match')
    if human_and_decoy_filter.O.decoy_read_id_list.shape[0] != human_and_decoy_filter.O.decoy_best_align_list.shape[0]:
        raise RuntimeError('Decoy reads and align list not match')
    """

    # filter query file

    human_and_decoy_filter.O.query_filename_list['path'] = human_and_decoy_filter.O.query_filename_list['path'].map(lambda x: os.path.join(megapath_nano.temp_dir_name, os.path.split(x)[1] + '.human_and_decoy_filtered'))

    read_id_to_keep_filename = os.path.join(RAM_dir_name, 'read_id_to_keep')
    human_and_decoy_filter.O.microbe_read_id_list.to_csv(path_or_buf=read_id_to_keep_filename, sep='\t', header=False, index=False, columns=['read_id'])

    filter_process_list = []
    output_file_list = []

    for file_index, filename in enumerate(human_and_decoy_filter.I.query_filename_list['path']):
        output_file = os.open(human_and_decoy_filter.O.query_filename_list['path'][file_index], flags=os.O_CREAT | os.O_EXCL | os.O_WRONLY, mode=0o644)
        filter_command = ['seqtk',]
        filter_command.extend(('subseq', filename, read_id_to_keep_filename))
        filter_process = subprocess.Popen(filter_command, close_fds=True, stdout=output_file)
        filter_process_list.append(filter_process)
        output_file_list.append(output_file)

    for index, filter_process in enumerate(filter_process_list):
        filter_process.wait()
        if filter_process.returncode != 0:
            raise RuntimeError('Error encountered in human_and_decoy_filter', human_and_decoy_filter.I.query_filename_list['path'][index])
        os.close(output_file_list[index])


    if megapath_nano.global_options['debug'] == False:
        shutil.rmtree(temp_dir_name)
        shutil.rmtree(RAM_dir_name)

    megapath_nano.log.print('Number of input reads: {num_read}'.format(num_read=human_and_decoy_filter.I.read_id_list.shape[0]))
    megapath_nano.log.print('Human read filtering:')
    megapath_nano.log.print('- number of reads with alignment score >= {threshold} : {num_read}'.format(threshold=human_and_decoy_filter.I.human_min_alignment_score, num_read=num_human_read_with_alignment_score_above_threshold))
    megapath_nano.log.print('- number of reads with alignment score / read length >= {threshold:.0%} : {num_read}'.format(threshold=human_and_decoy_filter.I.human_min_alignment_score_percent/100, num_read=num_human_read_with_normalized_alignment_score_above_threshold))
    megapath_nano.log.print('- number of reads satisfying either criteria: {num_read}'.format(num_read=human_and_decoy_filter.O.human_read_id_list.shape[0]))
    megapath_nano.log.print('Decoy read filtering:')
    megapath_nano.log.print('- number of reads with alignment score >= {threshold} : {num_read}'.format(threshold=human_and_decoy_filter.I.decoy_min_alignment_score, num_read=num_decoy_read_with_alignment_score_above_threshold))
    megapath_nano.log.print('- number of reads with alignment score / read length >= {threshold:.0%} : {num_read}'.format(threshold=human_and_decoy_filter.I.decoy_min_alignment_score_percent/100, num_read=num_decoy_read_with_normalized_alignment_score_above_threshold))
    megapath_nano.log.print('- number of reads satisfying either criteria: {num_read}'.format(num_read=human_and_decoy_filter.O.decoy_read_id_list.shape[0]))
    megapath_nano.log.print('Number of reads remaining: {num_read}'.format(num_read=human_and_decoy_filter.O.microbe_read_id_list.shape[0]))
    megapath_nano.log.print('end')


def step_placement_to_species(megapath_nano, placement_to_species):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()
    megapath_nano.log.print('Species identification threshold (aligned bp): {min_aligned_bp} '.format(min_aligned_bp=placement_to_species.I.species_ID_min_aligned_bp))

    RAM_dir_name = tempfile.mkdtemp(prefix='placement_to_species.', dir=megapath_nano.RAM_dir_name)

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)
    placement_to_species.O.align_list = Align(
                                              assembly_metadata=megapath_nano.assembly_metadata,
                                              global_options=megapath_nano.global_options,
                                              temp_dir_name=RAM_dir_name,
                                              log_file=megapath_nano.aligner_log,
                                              query_filename_list=placement_to_species.I.query_filename_list,
                                              target_assembly_list=placement_to_species.I.target_assembly_list,
                                              aligner_options=shlex.split(megapath_nano.global_options['alignerThreadOption'] + ' -N 50 -p 1 -x map-ont --split-prefix tmp'),
                                              paf_path_and_prefix=placement_to_species.I.paf_path_and_prefix,
                                              mapping_only=megapath_nano.global_options['mapping_only'],
                                              taxon_and_AMR_module_option=FLAGS.taxon_and_AMR_module_option,
                                              align_concat_fa=True,
                                              AMR_output_folder=f'{megapath_nano.output_folder}/{file_prefix_with_path}_amr/')
    if megapath_nano.global_options['debug'] == True:
        placement_to_species.O.align_list.to_csv(path_or_buf=file_prefix_with_path + '.species_align_list', sep='\t', header=True, index=False)

    if FLAGS.taxon_and_AMR_module_option=='AMR_module_only':
        os.sys.exit('Finished alignment.')
    if FLAGS.reassignment == True:
        placement_to_species.O.align_list=Reassign(align_list=placement_to_species.O.align_list,
                                                   db_folder=megapath_nano.global_options["db_folder"],
                                                   threads=min(psutil.cpu_count(logical=True), megapath_nano.global_options['max_aligner_thread']),
                                                   level=FLAGS.resolution)

    placement_to_species.O.best_align_list = placement_to_species.O.align_list.sort_values(['read_id', 'alignment_score', 'alignment_score_tiebreaker']).drop_duplicates(subset=['read_id'], keep='last').copy()

    placement_to_species.O.aligned_species_list = placement_to_species.O.best_align_list.assign(aligned_bp = lambda x: x['sequence_to'] - x['sequence_from']).groupby(['species_tax_id'], as_index=False).sum()[['species_tax_id', 'aligned_bp']]
    placement_to_species.O.selected_species_list = placement_to_species.O.aligned_species_list.query('aligned_bp >= @placement_to_species.I.species_ID_min_aligned_bp')

    placement_to_species.O.read_id_species_id = placement_to_species.O.best_align_list.merge(
                                                                                             right=placement_to_species.O.selected_species_list.set_index(['species_tax_id']),
                                                                                             how='inner', 
                                                                                             left_on='species_tax_id', 
                                                                                             right_index = True,
                                                                                             suffixes=['', '_y'],
                                                                                             validate='m:1',
                                                                                            )[['read_id', 'species_tax_id']]

    if placement_to_species.O.selected_species_list.shape[0] <= 0:
        os.sys.exit('No species can be identified, or the alignment output is empty in case the memory is not enough for alignment.')

    if megapath_nano.global_options['debug'] == False:
        shutil.rmtree(RAM_dir_name)

    megapath_nano.log.print('number of alignment: {num_alignment}'.format(num_alignment=placement_to_species.O.best_align_list.shape[0]))
    megapath_nano.log.print('number of species aligned: {num_species_aligned}'.format(num_species_aligned=placement_to_species.O.aligned_species_list.shape[0]))
    megapath_nano.log.print('number of species selected: {num_species_selected}'.format(num_species_selected=placement_to_species.O.selected_species_list.shape[0]))
    megapath_nano.log.print('end')


def step_placement_to_assembly(megapath_nano, placement_to_assembly):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()

    temp_dir_name = tempfile.mkdtemp(prefix='placement_to_assembly.', dir=megapath_nano.temp_dir_name)
    RAM_dir_name = tempfile.mkdtemp(prefix='placement_to_assembly.', dir=megapath_nano.RAM_dir_name)

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    target_assembly_list = placement_to_assembly.I.target_assembly_list.merge(
                                                                              right=placement_to_assembly.I.species_id_assembly_id.assign(aligned = lambda x: 1).set_index('assembly_id'),
                                                                              how='left', 
                                                                              left_on='assembly_id',
                                                                              right_index = True,
                                                                              suffixes=['', '_y'],
                                                                              validate='1:1',
                                                                             ).fillna(0).query('aligned == 0').drop(['aligned'], axis=1)

    species_target_assembly_list = megapath_nano.assembly_metadata.get_tax_id(assembly_list=target_assembly_list).merge(
                                                                                                                   right=placement_to_assembly.I.species_list[['species_tax_id']].set_index('species_tax_id'),
                                                                                                                   how='inner', 
                                                                                                                   left_on='species_tax_id',
                                                                                                                   right_index = True,
                                                                                                                   suffixes=['', '_y'],
                                                                                                                   validate='m:1',
                                                                                                                  )
    
    species_list = species_target_assembly_list[['species_tax_id']].drop_duplicates()

    read_id_species_id = placement_to_assembly.I.read_id_species_id.merge(
                                                                          right=species_list.set_index('species_tax_id'),
                                                                          how='inner', 
                                                                          left_on='species_tax_id',
                                                                          right_index = True,
                                                                          suffixes=['', '_y'],
                                                                          validate='m:1',
                                                                         )

    read_query = read_id_species_id.rename(columns={'species_tax_id': 'query_filename'})
    read_query['query_filename'] = read_query['query_filename'].map(lambda x: os.path.join(temp_dir_name, str(x)))
    read_query_filename = os.path.join(RAM_dir_name, 'read_query')
    read_query.to_csv(path_or_buf=read_query_filename, sep='\t', na_rep='NaN', header=False, index=False, columns=['read_id', 'query_filename'])

    split_command = [os.path.join(megapath_nano.global_options['tool_folder'], 'nanosplit'),]
    split_command.append(read_query_filename)
    for query_filename in placement_to_assembly.I.query_filename_list['path']:
        split_command.append(query_filename)
    split_process = subprocess.Popen(split_command, close_fds=True)
    split_process.wait()
    if split_process.returncode != 0:
        raise RuntimeError('Error encountered in placement_to_assembly')

    placement_to_assembly.O.align_list = pandas.DataFrame(columns=align_list_col_name)
    for key, value in align_list_col_type.items():
        placement_to_assembly.O.align_list[key] = placement_to_assembly.O.align_list[key].astype(value)

    placement_to_assembly.O.num_assembly_candidate = species_target_assembly_list.shape[0]

    for species in species_list.itertuples():

        species_query_filename = os.path.join(temp_dir_name, str(species.species_tax_id))

        species_align_list = Align(
                                   assembly_metadata=megapath_nano.assembly_metadata,
                                   global_options=megapath_nano.global_options,
                                   temp_dir_name=RAM_dir_name,
                                   log_file=megapath_nano.aligner_log,
                                   query_filename_list=pandas.DataFrame([species_query_filename,], columns=['path']),
                                   target_assembly_list=species_target_assembly_list.query('species_tax_id == @species.species_tax_id'),
                                   aligner_options=shlex.split(megapath_nano.global_options['alignerThreadOption'] + ' -N 1000 -p 0 -x map-ont'),
                                   mapping_only=megapath_nano.global_options['mapping_only'],
                                  )

        placement_to_assembly.O.align_list = pandas.concat([placement_to_assembly.O.align_list, species_align_list, ], axis=0,sort=True)

    if megapath_nano.global_options['debug'] == True:
        placement_to_assembly.O.align_list.to_csv(path_or_buf=file_prefix_with_path + '.assembly_align_list', sep='\t', header=True, index=False)

    if megapath_nano.global_options['debug'] == False:
        shutil.rmtree(temp_dir_name)
        shutil.rmtree(RAM_dir_name)

    megapath_nano.log.print('number of additional assembly candidates: {num_candidate}'.format(num_candidate=placement_to_assembly.O.num_assembly_candidate))
    megapath_nano.log.print('end')


def step_assembly_selection(megapath_nano, assembly_selection):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()

    megapath_nano.log.print('Assembly ID will be done on species ID genome set and assembly ID genome set for a species if average depth >= {avg_depth}'.format(avg_depth=assembly_selection.I.assembly_ID_min_average_depth))
    megapath_nano.log.print('Assembly ID will be done on species ID genome set otherwise')
  
    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    species_align_stat = align_list_to_align_stat_by_assembly_id(
                                                                 assembly_metadata=megapath_nano.assembly_metadata,
                                                                 log=megapath_nano.log,
                                                                 align_list=assembly_selection.I.species_align_list,
                                                                )
    species_align_stat = species_align_stat.sort_values(['species_tax_id', 'adjusted_average_depth', 'alignment_score_tiebreaker']).drop_duplicates(subset=['species_tax_id'], keep='last')
    species_align_stat = species_align_stat.merge(
                                                  right=assembly_selection.I.species_list[['species_tax_id']].set_index('species_tax_id'),
                                                  how='inner',
                                                  left_on='species_tax_id',
                                                  right_index = True,
                                                  suffixes=['', '_y'],
                                                  validate='1:1',
                                                 )

    species_reached_min_average_depth = species_align_stat.query('adjusted_average_depth >= @assembly_selection.I.assembly_ID_min_average_depth')[['species_tax_id']]
    species_not_reached_min_average_depth = assembly_selection.I.species_list.merge(
                                                                                    right=species_reached_min_average_depth.assign(reached_min_average_depth = lambda x: 1).set_index('species_tax_id'),
                                                                                    how='left',
                                                                                    left_on='species_tax_id',
                                                                                    right_index = True,
                                                                                    suffixes=['', '_y'],
                                                                                    validate='1:1',
                                                                                   ).fillna(0).query('reached_min_average_depth == 0')[['species_tax_id']]

    assembly_align_list = assembly_selection.I.assembly_align_list.merge(
                                                                         right=species_reached_min_average_depth.set_index('species_tax_id'),
                                                                         how='inner',
                                                                         left_on='species_tax_id',
                                                                         right_index = True,
                                                                         suffixes=['', '_y'],
                                                                         validate='m:1',
                                                                        )

    species_align_list = assembly_selection.I.species_align_list.merge(
                                                                       right=assembly_selection.I.read_id_species_id.rename(columns={'species_tax_id': 'read_species_tax_id'}).set_index('read_id'),
                                                                       how='inner',
                                                                       left_on='read_id',
                                                                       right_index = True,
                                                                       suffixes=['', '_y'],
                                                                       validate='m:1',
                                                                      ).query('species_tax_id == read_species_tax_id').drop(['read_species_tax_id'], axis=1)

    assembly_selection.O.align_list = pandas.concat([assembly_align_list, species_align_list], axis=0,sort=True)

    if megapath_nano.global_options['debug'] == True:
        species_align_stat.to_csv(path_or_buf=file_prefix_with_path + '.species_align_stat', sep='\t', header=True, index=False)
        assembly_selection.O.align_list.to_csv(path_or_buf=file_prefix_with_path + '.selection_align_list', sep='\t', header=True, index=False)

    assembly_selection.O.best_align_list = assembly_selection.O.align_list.sort_values(['read_id', 'alignment_score', 'alignment_score_tiebreaker']).drop_duplicates(subset=['read_id'], keep='last')
    assembly_selection.O.good_align_list = good_align_list(
    	                                                     align_list=assembly_selection.O.align_list,
    	                                                     good_align_threshold=assembly_selection.I.good_align_threshold,
    	                                                    )

    assembly_selection.O.align_stat = align_list_to_align_stat_by_assembly_id(
        				                                                              assembly_metadata=megapath_nano.assembly_metadata,
        				                                                              log=megapath_nano.log,
        				                                                              align_list=assembly_selection.O.good_align_list,
        				                                                             )

    assembly_selection.O.assembly_list = assembly_selection.O.align_stat.sort_values(['species_tax_id', 'adjusted_average_depth', 'alignment_score_tiebreaker']).drop_duplicates(subset=['species_tax_id'], keep='last').copy()


    megapath_nano.log.print('Number of species with assembly selected from species ID genome set and assembly ID genome set: {num_assembly}'.format(num_assembly=species_reached_min_average_depth.shape[0]))
    megapath_nano.log.print('Number of species with assembly selected from species ID genome set only: {num_species}'.format(num_species=species_not_reached_min_average_depth.shape[0]))
    megapath_nano.log.print('end')

def step_align_assembly_set(megapath_nano, align_assembly_set):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()
    
    RAM_dir_name = tempfile.mkdtemp(prefix='placement_to_assembly.', dir=megapath_nano.RAM_dir_name)

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    target_assembly_list = align_assembly_set.I.assembly_list.merge(
                                                                    right=align_assembly_set.I.species_id_assembly_id[['assembly_id']].assign(aligned = lambda x: 1).set_index('assembly_id'),
                                                                    how='left', 
                                                                    left_on='assembly_id',
                                                                    right_index = True,
                                                                    suffixes=['', '_y'],
                                                                    validate='1:1',
                                                                   ).fillna(0).query('aligned == 0').drop(['aligned'], axis=1)

    if megapath_nano.global_options['debug'] == True:
        target_assembly_list.to_csv(path_or_buf=file_prefix_with_path + '.target_assembly_list', sep='\t', header=True, index=False)

    if target_assembly_list.empty == False:
        align_list = Align(
                           assembly_metadata=megapath_nano.assembly_metadata, 
                           global_options=megapath_nano.global_options,
                           temp_dir_name=RAM_dir_name,
                           log_file=megapath_nano.aligner_log,
                           query_filename_list=align_assembly_set.I.query_filename_list,
                           target_assembly_list=target_assembly_list,
                           aligner_options=shlex.split(megapath_nano.global_options['alignerThreadOption'] + ' -N 1000 -p 0 -x map-ont'),
                           paf_path_and_prefix=align_assembly_set.I.paf_path_and_prefix,
                           mapping_only=megapath_nano.global_options['mapping_only'],
                          )
    else:
        align_list = pandas.DataFrame(columns=align_list_col_name)
        for key, value in align_list_col_type.items():
            align_list[key] = align_list[key].astype(value)

    if megapath_nano.global_options['debug'] == True:
        align_list.to_csv(path_or_buf=file_prefix_with_path + '.align_assembly_align_list', sep='\t', header=True, index=False)

    species_align_list = align_assembly_set.I.species_align_list.merge(
                                                                       right=align_assembly_set.I.assembly_list[['assembly_id']].set_index('assembly_id'),
                                                                       how='inner', 
                                                                       left_on='assembly_id',
                                                                       right_index = True,
                                                                       suffixes=['', '_y'],
                                                                       validate='m:1',
                                                                      )

    align_assembly_set.O.align_list = pandas.concat([align_list, species_align_list], axis=0,sort=True)

    if megapath_nano.global_options['debug'] == True:
        align_assembly_set.O.align_list.to_csv(path_or_buf=file_prefix_with_path + '.align_assembly_concat_align_list', sep='\t', header=True, index=False)

    #align_assembly_set.O.best_align_list = align_assembly_set.O.align_list.sort_values(['read_id', 'alignment_score', 'alignment_score_tiebreaker']).drop_duplicates(subset=['read_id'], keep='last')
    align_assembly_set.O.best_align_list = align_list_to_best_align_list(
                                                                         assembly_metadata=megapath_nano.assembly_metadata,
                                                                         log=megapath_nano.log,
                                                                         align_list=align_assembly_set.O.align_list,
                                                                        )

    if megapath_nano.global_options['debug'] == False:
        shutil.rmtree(RAM_dir_name)

    megapath_nano.log.print('end')


def step_raw_stat(megapath_nano, raw_stat):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()
    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)
    raw_align_list = raw_stat.I.best_align_list.merge(
                                                      right=raw_stat.I.huamn_and_decoy_best_align_list[['read_id', 'alignment_score']].rename(columns={'alignment_score': 'human_and_decoy_best_alignment_score'}).set_index('read_id'),
                                                      how='left', 
                                                      left_on='read_id', 
                                                      right_index = True,
                                                      suffixes=['', '_y'],
                                                      validate='m:1',
                                                     ).fillna(0).query('alignment_score >= human_and_decoy_best_alignment_score').drop(['human_and_decoy_best_alignment_score'], axis=1)

    if megapath_nano.global_options['debug'] == True:
        raw_align_list.to_csv(path_or_buf=file_prefix_with_path + '.raw_align_list', sep='\t', header=True, index=False)

    raw_stat.O.best_align_list = align_list_to_best_align_list(
                                                               assembly_metadata=megapath_nano.assembly_metadata,
                                                               log=megapath_nano.log,
                                                               align_list=raw_align_list,
                                                              )

    megapath_nano.log.print('end')


def step_variable_region(megapath_nano, variable_region):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()

    variable_region_divergence = 10

    RAM_dir_name = tempfile.mkdtemp(prefix='variable_region.', dir=megapath_nano.RAM_dir_name)

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    col_name = ('sequence_id', 'sequence_length', 'sequence_from', 'sequence_to', 'alignment_score')
    col_type = {'sequence_id': str, 'sequence_length' : int, 'sequence_from' : int, 'sequence_to': int, 'alignment_score' : int}

    temp_align_list_filename = os.path.join(RAM_dir_name, 'temp_align_list_filename')

    bed_list_col_name = ['target_sequence_id', 'target_sequence_from', 'target_sequence_to', 'target_assembly_id']
    bed_list_col_type = {'target_sequence_id': str, 'target_sequence_from': int, 'target_sequence_to': int, 'target_assembly_id': str}

    empty_bed_list = pandas.DataFrame(columns=bed_list_col_name)
    for key, value in bed_list_col_type.items():
        empty_bed_list[key] = empty_bed_list[key].astype(value)

    all_bed_list = empty_bed_list.copy()

    species_id_selected_assembly_id = megapath_nano.assembly_metadata.get_assembly_path(assembly_list=megapath_nano.assembly_metadata.get_tax_id(assembly_list=variable_region.I.selected_assembly_list))
    species_id_target_assembly_id = megapath_nano.assembly_metadata.get_assembly_path(assembly_list=megapath_nano.assembly_metadata.get_tax_id(assembly_list=variable_region.I.target_assembly_list))


    num_assembly_id = species_id_selected_assembly_id.shape[0]
    for assembly_id_index in range(num_assembly_id):

        assembly_id = species_id_selected_assembly_id.iloc[assembly_id_index]['assembly_id']
        species_tax_id = species_id_selected_assembly_id.iloc[assembly_id_index]['species_tax_id']
        assembly_path = os.path.join(megapath_nano.global_options['assembly_folder'], species_id_selected_assembly_id.iloc[assembly_id_index]['path'])

        query_species_id_assembly_id = species_id_target_assembly_id.query('species_tax_id == @species_tax_id and assembly_id != @assembly_id')

        if query_species_id_assembly_id.empty == True:
            continue

        species_bed_list = None

        temp_index_filename = os.path.join(RAM_dir_name, assembly_id + '.index')

        num_query_assembly_id = query_species_id_assembly_id.shape[0]
        for query_assembly_id_index in range(num_query_assembly_id):

            query_assembly_id = query_species_id_assembly_id.iloc[query_assembly_id_index]['assembly_id']
            query_assembly_path = os.path.join(megapath_nano.global_options['assembly_folder'], query_species_id_assembly_id.iloc[query_assembly_id_index]['path'])

            bed_filename = os.path.join(os.path.split(assembly_path)[0], query_assembly_id + '-' + assembly_id) + '.' + 'asm' + str(variable_region_divergence)

            bed_list = None

            if os.path.exists(bed_filename) == True:
                try:
                    bed_list = pandas.read_csv(bed_filename, 
                                               sep='\t', names=bed_list_col_name, usecols=bed_list_col_name, index_col=False, dtype=bed_list_col_type, header=None, skiprows=1)
                except pandas.io.common.EmptyDataError:
                    pass

            else:

                temp_align_list_file = os.open(temp_align_list_filename, flags=os.O_CREAT | os.O_EXCL | os.O_WRONLY, mode=0o644)

                aligner_command = [megapath_nano.global_options['aligner'],]
                aligner_command.extend(shlex.split(similarity_option(divergence=variable_region_divergence) + ' -N 1000 -p 0 -c'))
                aligner_command.extend(shlex.split(megapath_nano.global_options['alignerThreadOption']))
                if os.path.exists(temp_index_filename) == False:
                    aligner_command.extend(('-d', temp_index_filename,))
                    aligner_command.append(assembly_path)
                else:
                    aligner_command.append(temp_index_filename)
                aligner_command.append(query_assembly_path)
                aligner_process = subprocess.Popen(aligner_command, close_fds=True, stdout=subprocess.PIPE, stderr=megapath_nano.aligner_log)

                awk_command = ['awk', ]
                awk_command.append('{OFS="\t"};{gsub("AS:i:","")};{print $6,$7,$8,$9,$15}')
                awk_process = subprocess.Popen(awk_command, close_fds=True, stdin=aligner_process.stdout, stdout=temp_align_list_file)

                aligner_process.wait()
                awk_process.wait()
                if aligner_process.returncode != 0:
                    raise RuntimeError('Error encountered when aligning assembly', aligner_command)
                if awk_process.returncode != 0:
                    raise RuntimeError('Error encountered when processing aligner result', awk_command)
                os.close(temp_align_list_file)

                try:
                    align_list = pandas.read_csv(temp_align_list_filename, 
                                                 sep='\t', names=col_name, index_col=False, dtype=col_type, header=None)
                    min_alignment_score = megapath_nano.global_options['min_alignment_score']
                    if align_list.empty == False:
                        align_list = align_list.query('sequence_length > 0 and alignment_score >= @min_alignment_score')
                    if align_list.empty == False:
                        temp_bed = BedTool.from_dataframe(align_list[['sequence_id', 'sequence_from', 'sequence_to']])
                        merged_bed = temp_bed.sort().merge()
                        bed_list = merged_bed.to_dataframe()
                        os.remove(merged_bed.fn)
                        bed_list = bed_list.rename(columns={'chrom': 'target_sequence_id', 'start': 'target_sequence_from', 'end': 'target_sequence_to'}).assign(target_assembly_id = lambda x: assembly_id)
                        bed_list.to_csv(path_or_buf=bed_filename, sep='\t', header=True, index=False, columns=bed_list_col_name)
                    else:
                        empty_bed_list.to_csv(path_or_buf=bed_filename, sep='\t', header=True, index=False, columns=bed_list_col_name)

                except pandas.errors.ParserError:
                    empty_bed_list.to_csv(path_or_buf=bed_filename, sep='\t', header=True, index=False, columns=bed_list_col_name)

                os.remove(temp_align_list_filename)
            
            if bed_list is not None and bed_list.empty == False:
                if species_bed_list is None:
                    species_bed_list = bed_list
                else:
                    species_bed_list = pandas.concat([species_bed_list, bed_list], axis=0)
        
        if os.path.exists(temp_index_filename) == True:
            os.remove(temp_index_filename)

        if species_bed_list is not None and species_bed_list.empty == False:
            if all_bed_list is None:
                all_bed_list = species_bed_list
            else:
                all_bed_list = pandas.concat([all_bed_list, species_bed_list], axis=0)

    variable_region.O.noise_stat = variable_region.I.selected_assembly_list[['assembly_id', 'assembly_length']].copy()

    if all_bed_list.empty == False:
        num_assembly_id = species_id_target_assembly_id.assign(num_assembly_id = lambda x: 1).groupby(['species_tax_id'], as_index=False).sum()
        num_assembly_id = num_assembly_id.merge(
                                                right=species_id_selected_assembly_id.set_index('species_tax_id'),
                                                how='inner', 
                                                left_on='species_tax_id', 
                                                right_index = True,
                                                suffixes=['', '_y'],
                                                validate='1:1',
                                               )
        max_depth = num_assembly_id.assign(max_depth = lambda x: x['num_assembly_id'] * variable_region.I.variable_region_percent / 100)

        all_bed_list = all_bed_list.rename(columns={
                                                    'target_assembly_id': 'assembly_id',
                                                    'target_sequence_id': 'sequence_id',
                                                    'target_sequence_from': 'sequence_from',
                                                    'target_sequence_to': 'sequence_to',
                                                   })

        # using max_sequence_length is actually incorrect as variable regions at end of seqeunces can be omitted; should use get_sequence_length instead
        max_sequence_length = all_bed_list.groupby(['sequence_id']).max()[['sequence_to']].rename(columns={'sequence_to': 'sequence_length'})
        all_bed_list = all_bed_list.merge(
                                          right=max_sequence_length,
                                          how='inner', 
                                          left_on='sequence_id', 
                                          right_index = True,
                                          suffixes=['', '_y'],
                                          validate='m:1',
                                         )

        variable_region.O.noise_bed, noise_span_bp = align_list_to_depth_bed(
                                                                             temp_dir_name=RAM_dir_name,
                                                                             align_list=all_bed_list,
                                                                             max_depth=max_depth,
                                                                             can_equal_to_max=False,
                                                                            )

        variable_region.O.noise_stat = variable_region.O.noise_stat.merge(
                                                                          right=noise_span_bp.set_index('assembly_id').rename(columns={'span_bp': 'variable_span_bp'}),
                                                                          how='inner', 
                                                                          left_on='assembly_id', 
                                                                          right_index = True,
                                                                          suffixes=['', '_y'],
                                                                          validate='1:1',
                                                                         )
    else:
        variable_region.O.noise_bed = BedTool('', from_string=True)
        variable_region.O.noise_stat = variable_region.O.noise_stat.assign(variable_span_bp = lambda x: 0)

    variable_region.O.noise_stat = variable_region.O.noise_stat.assign(variable_span_percent = lambda x: x['variable_span_bp'] / x['assembly_length'])

    variable_region.O.noise_summary = variable_region.O.noise_stat['variable_span_percent'].describe()

    if megapath_nano.global_options['debug'] == False:
        shutil.rmtree(RAM_dir_name)

    megapath_nano.log.print('Variable region summary - average variable span: {average:.2%} - max variable span: {max:.2%}'.format(average=variable_region.O.noise_summary.loc['mean'], max=variable_region.O.noise_summary.loc['max']))
    megapath_nano.log.print('end')


def step_spike_filter(megapath_nano, spike_filter):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()

    RAM_dir_name = tempfile.mkdtemp(prefix='spike_filter.', dir=megapath_nano.RAM_dir_name)

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    spike_filter.O.noise_stat = align_list_to_align_stat_by_assembly_id(
        				                                                        assembly_metadata=megapath_nano.assembly_metadata,
        				                                                        log=megapath_nano.log,
        				                                                        align_list=spike_filter.I.align_list,
        				                                                       )

    spike_filter.O.noise_stat = spike_filter.O.noise_stat.assign(expected_max_depth = lambda x: x['adjusted_average_depth'] + spike_filter.I.expected_max_depth_stdev * numpy.sqrt(x['adjusted_average_depth']))
    spike_filter.O.noise_stat['expected_max_depth'] = spike_filter.O.noise_stat['expected_max_depth'].astype('int')
    spike_filter.O.noise_stat['expected_max_depth'] = numpy.where(spike_filter.O.noise_stat['expected_max_depth'] < 1, 1, spike_filter.O.noise_stat['expected_max_depth'])

    if megapath_nano.global_options['debug'] == True:
        spike_filter.O.noise_stat.to_csv(path_or_buf=file_prefix_with_path + '.spike_noise_stat', sep='\t', header=True, index=False)


    spike_filter.O.noise_bed, noise_span_bp = align_list_to_depth_bed(
                                                                      temp_dir_name=RAM_dir_name,
                                                                      align_list=spike_filter.I.align_list,
                                                                      min_depth=spike_filter.O.noise_stat.rename(columns={'expected_max_depth':'min_depth'}),
                                                                      can_equal_to_min=False,
                                                                     )

    spike_filter.O.noise_stat = spike_filter.O.noise_stat.merge(
                                                                right=noise_span_bp.set_index('assembly_id').rename(columns={'span_bp': 'spike_span_bp'}),
                                                                how='inner', 
                                                                left_on='assembly_id', 
                                                                right_index = True,
                                                                suffixes=['', '_y'],
                                                                validate='1:1',
                                                               )

    spike_filter.O.noise_stat = spike_filter.O.noise_stat.assign(spike_span_percent = lambda x: x['spike_span_bp'] / x['assembly_length'])[['assembly_id', 'spike_span_bp', 'spike_span_percent']]

    spike_filter.O.noise_summary = spike_filter.O.noise_stat['spike_span_percent'].describe()

    if megapath_nano.global_options['debug'] == False:
        shutil.rmtree(RAM_dir_name)

    megapath_nano.log.print('filter summary - average noise span: {average:.2%} - max noise span: {max:.2%}'.format(average=spike_filter.O.noise_summary.loc['mean'], max=spike_filter.O.noise_summary.loc['max']))
    megapath_nano.log.print('end')


def step_human_repetitive_region_filter(megapath_nano, human_repetitive_region_filter):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()

    RAM_dir_name = tempfile.mkdtemp(prefix='human_repetitive_region_filter.', dir=megapath_nano.RAM_dir_name)

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    bed_list_col_name = ['target_sequence_id', 'target_sequence_from', 'target_sequence_to', 'target_assembly_id']
    bed_list_col_type = {'target_sequence_id': str, 'target_sequence_from': int, 'target_sequence_to': int, 'target_assembly_id': str}


    file_extension = 'asm20'
    aligner_option = shlex.split(similarity_option(divergence=20) + ' ' + megapath_nano.global_options['alignerThreadOption'] + ' -N 10000 -p 0')

    assembly_path = megapath_nano.assembly_metadata.get_assembly_path(assembly_list=human_repetitive_region_filter.I.assembly_list)
    assembly_path['path'] = assembly_path['path'].map(lambda x: os.path.join(megapath_nano.global_options['assembly_folder'], x))

    num_assembly_id = assembly_path.shape[0]

    assembly_path = assembly_path.assign(with_bed = lambda x: 0).copy()

    all_bed_list = pandas.DataFrame(columns=bed_list_col_name)
    for key, value in bed_list_col_type.items():
        all_bed_list[key] = all_bed_list[key].astype(value)

    for index in range(num_assembly_id):

        bed_filename = os.path.join(os.path.split(assembly_path.iloc[index]['path'])[0], megapath_nano.global_options['human_repetitive_region_filter_assembly_id'] + '-' + assembly_path.iloc[index]['assembly_id']) + '.' + file_extension

        bed_list = None
        if os.path.exists(bed_filename) == True:
            assembly_path.with_bed.iloc[index] = 1     # the update is done correctly but yet but give warning messages
            try:
                bed_list = pandas.read_csv(bed_filename, 
                                           sep='\t', names=bed_list_col_name, usecols=bed_list_col_name, index_col=False, dtype=bed_list_col_type, header=None, skiprows=1)
            except pandas.io.common.EmptyDataError:
                continue
        
            if bed_list is not None and bed_list.empty == False:
                all_bed_list = pandas.concat([all_bed_list, bed_list], axis=0)

    assembly_path_without_bed = assembly_path.query('with_bed == 0')

    num_assembly_id_without_bed = assembly_path_without_bed.shape[0]

    if num_assembly_id_without_bed > 0:
        align_list = Align(
                           assembly_metadata=megapath_nano.assembly_metadata,
                           global_options=megapath_nano.global_options,
                           temp_dir_name=RAM_dir_name,
                           log_file=megapath_nano.aligner_log,
                           query_assembly_list=pandas.DataFrame(data={'assembly_id': [megapath_nano.global_options['human_repetitive_region_filter_assembly_id'],]}),
                           target_assembly_list=assembly_path_without_bed,
                           aligner_options=aligner_option,
                           mapping_only=megapath_nano.global_options['mapping_only'],
                          )
        bed = align_list_to_bed(align_list=align_list)
        if bed.count() > 0:
            align_bed_list = bed.to_dataframe().rename(columns={
                                                                'chrom': 'target_sequence_id',
                                                                'start': 'target_sequence_from',
                                                                'end': 'target_sequence_to',
                                                                'name': 'target_assembly_id'
                                                               })
        else:
            align_bed_list = pandas.DataFrame(columns=bed_list_col_name)
            for key, value in bed_list_col_type.items():
                all_bed_list[key] = all_bed_list[key].astype(value)

        if align_bed_list is not None and align_bed_list.empty == False:
            all_bed_list = pandas.concat([all_bed_list, align_bed_list], axis=0)

        for index in range(num_assembly_id_without_bed):
            assembly_id = assembly_path_without_bed.iloc[index]['assembly_id']
            assembly_bed_list = align_bed_list.query('target_assembly_id == @assembly_id')
            align_bed_filename = os.path.join(os.path.split(assembly_path_without_bed.iloc[index]['path'])[0], megapath_nano.global_options['human_repetitive_region_filter_assembly_id'] + '-' + assembly_id) + '.' + file_extension
            print(align_bed_filename)

            assembly_bed_list.to_csv(path_or_buf=align_bed_filename, sep='\t', header=True, index=False, columns=bed_list_col_name)

    if all_bed_list.empty == False:
        human_repetitive_region_filter.O.noise_bed = align_list_to_bed(align_list=all_bed_list.rename(columns={
                                                                                                     'target_sequence_id': 'sequence_id',
                                                                                                     'target_sequence_from': 'sequence_from',
                                                                                                     'target_sequence_to': 'sequence_to',
                                                                                                     'target_assembly_id': 'assembly_id',
                                                                                                    }))
        covered_bp = bed_to_covered_bp_by_assembly_id(bed=human_repetitive_region_filter.O.noise_bed)
    else:
        human_repetitive_region_filter.O.noise_bed = BedTool('', from_string=True)
        covered_bp = pandas.DataFrame(columns=['assembly_id', 'covered_bp'])
        covered_bp['covered_bp'] = covered_bp['covered_bp'].astype('int')

    noise_span_bp = human_repetitive_region_filter.I.assembly_list[['assembly_id']].merge(
                                                                                right=covered_bp.rename(columns={'covered_bp':'noise_span_bp'}).set_index(['assembly_id']),
                                                                                how='left', 
                                                                                left_on='assembly_id',
                                                                                right_index = True,
                                                                                suffixes=['', '_y'],
                                                                                validate='1:1',
                                                                               ).fillna(0).copy()
    noise_span_bp['noise_span_bp'] = noise_span_bp['noise_span_bp'].astype('int')

    human_repetitive_region_filter.O.noise_stat = megapath_nano.assembly_metadata.get_assembly_length(assembly_list=human_repetitive_region_filter.I.assembly_list, how='left').fillna(0)
    human_repetitive_region_filter.O.noise_stat['assembly_length'] = human_repetitive_region_filter.O.noise_stat['assembly_length'].astype('int')

    human_repetitive_region_filter.O.noise_stat = human_repetitive_region_filter.O.noise_stat.merge(
                                                                                right=noise_span_bp.set_index('assembly_id').rename(columns={'noise_span_bp': 'human_span_bp'}),
                                                                                how='inner', 
                                                                                left_on='assembly_id', 
                                                                                right_index = True,
                                                                                suffixes=['', '_y'],
                                                                                validate='1:1',
                                                                               )
    human_repetitive_region_filter.O.noise_stat = human_repetitive_region_filter.O.noise_stat.assign(human_span_percent = lambda x: x['human_span_bp'] / x['assembly_length'])[['assembly_id', 'human_span_bp', 'human_span_percent']]

    human_repetitive_region_filter.O.noise_summary = human_repetitive_region_filter.O.noise_stat['human_span_percent'].describe()

    megapath_nano.log.print('filter summary - average noise span: {average:.2%} - max noise span: {max:.2%}'.format(average=human_repetitive_region_filter.O.noise_summary.loc['mean'], max=human_repetitive_region_filter.O.noise_summary.loc['max']))
    megapath_nano.log.print('end')


def step_approx_abundance(megapath_nano, approx_abundance):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()

    RAM_dir_name = tempfile.mkdtemp(prefix='approx_abundance.', dir=megapath_nano.RAM_dir_name)

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    align_list = select_alignment_by_bed(
                                         temp_dir_name=RAM_dir_name,
                                         align_list=approx_abundance.I.align_list,
                                         bed=approx_abundance.I.noise_bed,
                                         max_overlap=approx_abundance.I.max_align_noise_overlap,
                                        )

    best_align_list = align_list_to_best_align_list(
                                                    assembly_metadata=megapath_nano.assembly_metadata,
                                                    log=megapath_nano.log,
                                                    align_list=align_list,
                                                   )

    approx_abundance.O.align_stat = align_list_to_align_stat_by_assembly_id(
      			                                                                assembly_metadata=megapath_nano.assembly_metadata,
      			                                                                log=megapath_nano.log,
      			                                                                align_list=best_align_list,
      			                                                                noise_bed=approx_abundance.I.noise_bed,
      			                                                               )

    if megapath_nano.global_options['debug'] == False:
        shutil.rmtree(RAM_dir_name)

    megapath_nano.log.print('end')


def step_microbe_repetitive_region_filter(megapath_nano, microbe_repetitive_region_filter):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()

    megapath_nano.log.print('Apply 99.2% similarity filter when a species is {num_fold} times more abundant than another species (same genus only)'.format(num_fold=microbe_repetitive_region_filter.I.abundance_threshold_99_2))
    megapath_nano.log.print('Apply 99.0% similarity filter when a species is {num_fold} times more abundant than another species (same genus only)'.format(num_fold=microbe_repetitive_region_filter.I.abundance_threshold_99))
    megapath_nano.log.print('Apply 98.0% similarity filter when a species is {num_fold} times more abundant than another species (same genus only)'.format(num_fold=microbe_repetitive_region_filter.I.abundance_threshold_98))
    megapath_nano.log.print('Apply 95.0% similarity filter when a species is {num_fold} times more abundant than another species'.format(num_fold=microbe_repetitive_region_filter.I.abundance_threshold_95))
    megapath_nano.log.print('Apply 90.0% similarity filter when a species is {num_fold} times more abundant than another species'.format(num_fold=microbe_repetitive_region_filter.I.abundance_threshold_90))
    megapath_nano.log.print('Apply 80.0% similarity filter when a species is {num_fold} times more abundant than another species'.format(num_fold=microbe_repetitive_region_filter.I.abundance_threshold_80))
    megapath_nano.log.print('Targeted max span of similar region = {span:.0%}'.format(span=microbe_repetitive_region_filter.I.targeted_max_span_percent / 100))
    megapath_nano.log.print('Allowed max span of similar region = {span:.0%}'.format(span=microbe_repetitive_region_filter.I.allowed_max_span_percent / 100))


    RAM_dir_name = tempfile.mkdtemp(prefix='microbe_repetitive_region_filter.', dir=megapath_nano.RAM_dir_name)

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    abundance_measure = microbe_repetitive_region_filter.I.align_stat[['assembly_id', 'adjusted_average_depth']].sort_values(['adjusted_average_depth']).copy().reset_index()
    abundance_measure = megapath_nano.assembly_metadata.get_tax_id(assembly_list=abundance_measure, how='left').fillna(0)
    abundance_measure = megapath_nano.assembly_metadata.get_assembly_length(assembly_list=abundance_measure, how='left').fillna(0)

    num_assembly_id = abundance_measure.shape[0]

    microbe_repetitive_region_filter.O.noise_bed = BedTool('', from_string=True)
    source_target_noise_span_bp = pandas.DataFrame(columns=['source_assembly_id', 'target_assembly_id', 'microbe_span_bp'])
    source_target_noise_span_bp['microbe_span_bp'] = source_target_noise_span_bp['microbe_span_bp'].astype('int')


    assembly_path = megapath_nano.assembly_metadata.get_assembly_path(assembly_list=abundance_measure)
    assembly_path['compressed_path'] = assembly_path['path'].map(lambda x: os.path.join(megapath_nano.global_options['assembly_folder'], x))
    assembly_path['uncompressed_path'] = assembly_path['path'].map(lambda x: os.path.join(RAM_dir_name, os.path.split(x)[1]+ '.uncompressed'))
    assembly_path['index_path'] = assembly_path['path'].map(lambda x: os.path.join(RAM_dir_name, os.path.split(x)[1] + '.index'))

    col_name = ('sequence_id', 'sequence_length', 'sequence_from', 'sequence_to', 'alignment_score')
    col_type = {'sequence_id': str, 'sequence_length' : int, 'sequence_from' : int, 'sequence_to': int, 'alignment_score' : int}

    temp_align_list_filename = os.path.join(RAM_dir_name, 'temp_align_list_filename')

    bed_list_col_name = ['target_sequence_id', 'target_sequence_from', 'target_sequence_to', 'target_assembly_id']
    bed_list_col_type = {'target_sequence_id': str, 'target_sequence_from': int, 'target_sequence_to': int, 'target_assembly_id': str}

    all_bed_list = pandas.DataFrame(columns=bed_list_col_name + ['source_assembly_id',])
    for key, value in bed_list_col_type.items():
        all_bed_list[key] = all_bed_list[key].astype(value)
    all_bed_list['source_assembly_id'] = all_bed_list['source_assembly_id'].astype(value)
    empty_bed_list = all_bed_list.copy()

    # skip species with zero abundance after human similar and spike filter
    for min_low_abundance_index in range(num_assembly_id):
        if abundance_measure.iloc[min_low_abundance_index]['adjusted_average_depth'] > 0:
            break

    # iloc returns a transposed series if input is an integer; returns a dataframe if input is a list
    for high_abundance_index in range(num_assembly_id - 1, 0, -1):
        if abundance_measure.iloc[high_abundance_index]['adjusted_average_depth'] < abundance_measure.iloc[0]['adjusted_average_depth'] * microbe_repetitive_region_filter.I.abundance_threshold_98:
            break
        if abundance_measure.iloc[high_abundance_index]['adjusted_average_depth'] < microbe_repetitive_region_filter.I.min_average_depth:
            break

        high_abundance_uncompressed_exists = False

        for low_abundance_index in range(min_low_abundance_index, high_abundance_index):

            if abundance_measure.iloc[high_abundance_index]['adjusted_average_depth'] >= abundance_measure.iloc[low_abundance_index]['adjusted_average_depth'] * microbe_repetitive_region_filter.I.abundance_threshold_80:
                similarity_index = 1
            elif abundance_measure.iloc[high_abundance_index]['adjusted_average_depth'] >= abundance_measure.iloc[low_abundance_index]['adjusted_average_depth'] * microbe_repetitive_region_filter.I.abundance_threshold_90:
                similarity_index = 2
            elif abundance_measure.iloc[high_abundance_index]['adjusted_average_depth'] >= abundance_measure.iloc[low_abundance_index]['adjusted_average_depth'] * microbe_repetitive_region_filter.I.abundance_threshold_95:
                similarity_index = 3
            elif abundance_measure.iloc[high_abundance_index]['adjusted_average_depth'] >= abundance_measure.iloc[low_abundance_index]['adjusted_average_depth'] * microbe_repetitive_region_filter.I.abundance_threshold_98:
                if abundance_measure.iloc[high_abundance_index]['genus_tax_id'] != 0 and abundance_measure.iloc[high_abundance_index]['genus_height'] <= microbe_repetitive_region_filter.I.genus_height and abundance_measure.iloc[high_abundance_index]['genus_tax_id'] == abundance_measure.iloc[low_abundance_index]['genus_tax_id']:
                    similarity_index = 4
                else:
                    continue
            elif abundance_measure.iloc[high_abundance_index]['adjusted_average_depth'] >= abundance_measure.iloc[low_abundance_index]['adjusted_average_depth'] * microbe_repetitive_region_filter.I.abundance_threshold_99:
                if abundance_measure.iloc[high_abundance_index]['genus_tax_id'] != 0 and abundance_measure.iloc[high_abundance_index]['genus_height'] <= microbe_repetitive_region_filter.I.genus_height and abundance_measure.iloc[high_abundance_index]['genus_tax_id'] == abundance_measure.iloc[low_abundance_index]['genus_tax_id']:
                    similarity_index = 5
                else:
                    continue
            elif abundance_measure.iloc[high_abundance_index]['adjusted_average_depth'] >= abundance_measure.iloc[low_abundance_index]['adjusted_average_depth'] * microbe_repetitive_region_filter.I.abundance_threshold_99_2:
                if abundance_measure.iloc[high_abundance_index]['genus_tax_id'] != 0 and abundance_measure.iloc[high_abundance_index]['genus_height'] <= microbe_repetitive_region_filter.I.genus_height and abundance_measure.iloc[high_abundance_index]['genus_tax_id'] == abundance_measure.iloc[low_abundance_index]['genus_tax_id']:
                    similarity_index = 6
                else:
                    continue
            else:
                break


            while True:

                if similarity_index == 1:
                    option = similarity_option(divergence=20)
                    file_extension = 'asm' + str(20)
                elif similarity_index == 2:
                    option = similarity_option(divergence=10)
                    file_extension = 'asm' + str(10)
                elif similarity_index == 3:
                    option = similarity_option(divergence=5)
                    file_extension = 'asm' + str(5)
                elif similarity_index == 4:
                    option = similarity_option(divergence=2)
                    file_extension = 'asm' + str(2)
                elif similarity_index == 5:
                    option = similarity_option(divergence=1)
                    file_extension = 'asm' + str(1)
                elif similarity_index == 6:
                    option = similarity_option(divergence=0.8)
                    file_extension = 'asm' + '0_8'
                else:
                    raise RuntimeError('Logic Error')

                aligner_option = shlex.split(option + ' -N 1000 -p 0 -c')

                bed_filename = os.path.join(os.path.split(assembly_path.iloc[low_abundance_index]['compressed_path'])[0], abundance_measure.iloc[high_abundance_index]['assembly_id'] + '-' + abundance_measure.iloc[low_abundance_index]['assembly_id']) + '.' + file_extension

                bed_list = None

                if os.path.exists(bed_filename) == True:
                    try:
                        bed_list = pandas.read_csv(bed_filename, 
                                                   sep='\t', names=bed_list_col_name, usecols=bed_list_col_name, index_col=False, dtype=bed_list_col_type, header=None, skiprows=1)
                    except pandas.io.common.EmptyDataError:
                        pass
                else:
                    if high_abundance_uncompressed_exists == False:
                        uncompress_command = ['zcat', '-f', assembly_path.iloc[high_abundance_index]['compressed_path']]
                        uncompressed_assembly_file = os.open(assembly_path.iloc[high_abundance_index]['uncompressed_path'], flags=os.O_CREAT | os.O_EXCL | os.O_WRONLY, mode=0o644)

                        uncompress_process = subprocess.Popen(uncompress_command, close_fds=True, stdout=uncompressed_assembly_file)
                        uncompress_process.wait()
                        if uncompress_process.returncode != 0:
                            raise RuntimeError('Error encountered when uncompressing assembly', uncompress_command)
                        os.close(uncompressed_assembly_file)
                        high_abundance_uncompressed_exists = True

                    # align high abundance species to low abundance species

                    temp_align_list_file = os.open(temp_align_list_filename, flags=os.O_CREAT | os.O_EXCL | os.O_WRONLY, mode=0o644)

                    aligner_command = [megapath_nano.global_options['aligner'],]
                    aligner_command.extend(aligner_option)
                    aligner_command.extend(shlex.split(megapath_nano.global_options['alignerThreadOption']))
                    if os.path.exists(assembly_path.iloc[low_abundance_index]['index_path']) == False:
                        aligner_command.extend(('-d', assembly_path.iloc[low_abundance_index]['index_path'],))
                        aligner_command.append(assembly_path.iloc[low_abundance_index]['compressed_path'])
                    else:
                        aligner_command.append(assembly_path.iloc[low_abundance_index]['index_path'])
                    aligner_command.append(assembly_path.iloc[high_abundance_index]['uncompressed_path'])
                    aligner_process = subprocess.Popen(aligner_command, close_fds=True, stdout=subprocess.PIPE, stderr=megapath_nano.aligner_log)

                    awk_command = ['awk', ]
                    awk_command.append('{OFS="\t"};{gsub("AS:i:","")};{print $6,$7,$8,$9,$15}')
                    awk_process = subprocess.Popen(awk_command, close_fds=True, stdin=aligner_process.stdout, stdout=temp_align_list_file)

                    aligner_process.wait()
                    awk_process.wait()
                    if aligner_process.returncode != 0:
                        raise RuntimeError('Error encountered when aligning assembly', aligner_command)
                    if awk_process.returncode != 0:
                        raise RuntimeError('Error encountered when processing aligner result', awk_command)
                    os.close(temp_align_list_file)

                    try:
                        align_list = pandas.read_csv(temp_align_list_filename, 
                                                     sep='\t', names=col_name, index_col=False, dtype=col_type, header=None)
                        min_alignment_score = megapath_nano.global_options['min_alignment_score']
                        if align_list.empty == False:
                            align_list = align_list.query('sequence_length > 0 and alignment_score >= @min_alignment_score')
                        if align_list.empty == False:
                            temp_bed = BedTool.from_dataframe(align_list[['sequence_id', 'sequence_from', 'sequence_to']])
                            merged_bed = temp_bed.sort().merge()
                            bed_list = merged_bed.to_dataframe()
                            os.remove(merged_bed.fn)
                            bed_list = bed_list.rename(columns={'chrom': 'target_sequence_id', 'start': 'target_sequence_from', 'end': 'target_sequence_to'}).assign(target_assembly_id = lambda x: abundance_measure.iloc[low_abundance_index]['assembly_id'])
                            bed_list.to_csv(path_or_buf=bed_filename, sep='\t', header=True, index=False, columns=bed_list_col_name)
                        else:
                            empty_bed_list.to_csv(path_or_buf=bed_filename, sep='\t', header=True, index=False, columns=bed_list_col_name)

                    except pandas.errors.ParserError:
                        empty_bed_list.to_csv(path_or_buf=bed_filename, sep='\t', header=True, index=False, columns=bed_list_col_name)

                    os.remove(temp_align_list_filename)

                if bed_list is not None and bed_list.empty == False:
                    covered_bp = bed_list.assign(covered_bp = lambda x: x['target_sequence_to'] - x['target_sequence_from']).sum()['covered_bp']
                    covered_percent = covered_bp / abundance_measure.iloc[low_abundance_index]['assembly_length'] * 100
                    if covered_percent > microbe_repetitive_region_filter.I.targeted_max_span_percent:
                        if similarity_index >= 6:
                            if covered_percent > microbe_repetitive_region_filter.I.allowed_max_span_percent:
                                bed_list = None
                            break
                        else:
                            similarity_index += 1
                            continue

                break
 
            if bed_list is not None and bed_list.empty == False:
                bed_list = bed_list.assign(source_assembly_id = lambda x: abundance_measure.iloc[high_abundance_index]['assembly_id'])
                all_bed_list = pandas.concat([all_bed_list, bed_list], axis=0)

            # if bed_list is not None and bed_list.empty == False:
            #     temp_bed = BedTool.from_dataframe(bed_list[['target_sequence_id', 'target_sequence_from', 'target_sequence_to']])
            #     merged_bed = temp_bed.sort().merge()
            #     bed_list = merged_bed.to_dataframe()
            #     os.remove(merged_bed.fn)

            #     bed_list = bed_list.rename(columns={'chrom': 'target_sequence_id', 'start': 'target_sequence_from', 'end': 'target_sequence_to'}).assign(
            #                                                                                                                                              target_assembly_id = lambda x: abundance_measure.iloc[low_abundance_index]['assembly_id'],
            #                                                                                                                                              source_assembly_id = lambda x: abundance_measure.iloc[high_abundance_index]['assembly_id'],
            #                                                                                                                                             )
            #     covered_percent = covered_bp / abundance_measure.iloc[low_abundance_index]['assembly_length'] * 100
            #     if covered_percent > microbe_repetitive_region_filter.I.max_span_percent_overall:
            #         megapath_nano.log.print('Max span percent exceeded (overall): ' + abundance_measure.iloc[low_abundance_index]['assembly_id'] + ' - ' + str(covered_percent))
            #     else:
            #         all_bed_list = pandas.concat([all_bed_list, bed_list], axis=0)

        if high_abundance_uncompressed_exists == True:
            os.remove(assembly_path.iloc[high_abundance_index]['uncompressed_path'])

    if all_bed_list.empty == True:
        microbe_repetitive_region_filter.O.source_target_noise_span_bp = pandas.DataFrame(columns=['target_assembly_id', 'source_assembly_id', 'microbe_span_bp'])
        microbe_repetitive_region_filter.O.source_target_noise_span_bp['microbe_span_bp'] = microbe_repetitive_region_filter.O.source_target_noise_span_bp['microbe_span_bp'].astype('int')
        microbe_repetitive_region_filter.O.noise_bed = BedTool('', from_string=True)
        covered_bp = pandas.DataFrame(columns=['assembly_id', 'covered_bp'])
        covered_bp['covered_bp'] = covered_bp['covered_bp'].astype('int')
    else:
        if megapath_nano.global_options['debug'] == True:
            all_bed_list.to_csv(path_or_buf=file_prefix_with_path + '.all_bed_list', sep='\t', header=True, index=False)

        source_target_noise_bed = align_list_to_bed(align_list=all_bed_list.assign(
                                                                                   assembly_id = lambda x: x['target_assembly_id'] + ';' + x['source_assembly_id'],
                                                                                   sequence_id = lambda x: x['target_sequence_id'],
                                                                                   sequence_from = lambda x: x['target_sequence_from'],
                                                                                   sequence_to = lambda x: x['target_sequence_to'],
                                                                                  ))
        source_target_noise_bed_df = source_target_noise_bed.to_dataframe()
        os.remove(source_target_noise_bed.fn)
        source_target_noise_bed_df = pandas.concat([source_target_noise_bed_df['name'].str.split(';', n=1, expand=True).rename(columns={0:'target_assembly_id', 1:'source_assembly_id'}), source_target_noise_bed_df[['start', 'end']]], axis=1)

        microbe_repetitive_region_filter.O.source_target_noise_span_bp = source_target_noise_bed_df.assign(microbe_span_bp = lambda x: x['end'] - x['start'])
        microbe_repetitive_region_filter.O.source_target_noise_span_bp = microbe_repetitive_region_filter.O.source_target_noise_span_bp.groupby(['target_assembly_id', 'source_assembly_id'], as_index=False).sum()
        microbe_repetitive_region_filter.O.source_target_noise_span_bp = microbe_repetitive_region_filter.O.source_target_noise_span_bp[['target_assembly_id', 'source_assembly_id', 'microbe_span_bp']]

        microbe_repetitive_region_filter.O.noise_bed = align_list_to_bed(align_list=all_bed_list.rename(columns={
                                                                                                       'target_assembly_id': 'assembly_id',
                                                                                                       'target_sequence_id': 'sequence_id',
                                                                                                       'target_sequence_from': 'sequence_from',
                                                                                                       'target_sequence_to': 'sequence_to',
                                                                                                      }))
 
        covered_bp = bed_to_covered_bp_by_assembly_id(bed=microbe_repetitive_region_filter.O.noise_bed)
    
    microbe_repetitive_region_filter.O.noise_stat = megapath_nano.assembly_metadata.get_assembly_length(assembly_list=microbe_repetitive_region_filter.I.assembly_list, how='left').fillna(0)
    microbe_repetitive_region_filter.O.noise_stat['assembly_length'] = microbe_repetitive_region_filter.O.noise_stat['assembly_length'].astype('int')

    microbe_repetitive_region_filter.O.noise_stat = microbe_repetitive_region_filter.O.noise_stat.merge(
                                                                                  right=covered_bp.set_index('assembly_id').rename(columns={'covered_bp': 'microbe_span_bp'}),
                                                                                  how='left', 
                                                                                  left_on='assembly_id', 
                                                                                  right_index = True,
                                                                                  suffixes=['', '_y'],
                                                                                  validate='1:1',
                                                                                 ).fillna(0)
    microbe_repetitive_region_filter.O.noise_stat['microbe_span_bp'] = microbe_repetitive_region_filter.O.noise_stat['microbe_span_bp'].astype('int')

    microbe_repetitive_region_filter.O.noise_stat = microbe_repetitive_region_filter.O.noise_stat.assign(microbe_span_percent = lambda x: x['microbe_span_bp'] / x['assembly_length'])[['assembly_id', 'microbe_span_bp', 'microbe_span_percent']]

    microbe_repetitive_region_filter.O.noise_summary = microbe_repetitive_region_filter.O.noise_stat['microbe_span_percent'].describe()

    if megapath_nano.global_options['debug'] == False:
        shutil.rmtree(RAM_dir_name)

    megapath_nano.log.print('filter summary - average noise span: {average:.2%} - max noise span: {max:.2%}'.format(average=microbe_repetitive_region_filter.O.noise_summary.loc['mean'], max=microbe_repetitive_region_filter.O.noise_summary.loc['max']))
    megapath_nano.log.print('end')


def step_noise_removal(megapath_nano, noise_removal):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()
    megapath_nano.log.print('alignments with > {max_noise_overlap:.0%} overlap with noise regions are removed'.format(max_noise_overlap=noise_removal.I.max_align_noise_overlap / 100))

    RAM_dir_name = tempfile.mkdtemp(prefix='noise_removal.', dir=megapath_nano.RAM_dir_name)

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    noise_removal.O.num_align_before = noise_removal.I.align_list.shape[0]
    
    noise_removal.O.align_list = select_alignment_by_bed(
                                                         temp_dir_name=RAM_dir_name,
                                                         align_list=noise_removal.I.align_list,
                                                         bed=noise_removal.I.noise_bed,
                                                         max_overlap=noise_removal.I.max_align_noise_overlap,
                                                       )
    noise_removal.O.align_list = noise_removal.O.align_list.merge(
                                                                  right=noise_removal.I.non_zero_approx_abundence_assembly_id.set_index('assembly_id'),
                                                                  how='inner', 
                                                                  left_on='assembly_id', 
                                                                  right_index = True,
                                                                  suffixes=['', '_y'],
                                                                  validate='m:1',
                                                                 )

    if megapath_nano.global_options['debug'] == True:
        noise_removal.I.align_list.to_csv(path_or_buf=file_prefix_with_path + '.before_noise_align_list', sep='\t', header=True, index=False)
        noise_removal.O.align_list.to_csv(path_or_buf=file_prefix_with_path + '.after_noise_align_list', sep='\t', header=True, index=False)

    noise_removal.O.num_align_after = noise_removal.O.align_list.shape[0]

    noise_removal.O.noise_bed = noise_removal.I.noise_bed
    
    #  TODO fix num
    #megapath_nano.log.print('{num_removed} out of {num_align} alignments removed'.format(
    #                                                                                num_removed=noise_removal.O.num_align_before-noise_removal.O.num_align_after,
    #                                                                                num_align=noise_removal.O.num_align_before))

    if megapath_nano.global_options['debug'] == False:
        shutil.rmtree(RAM_dir_name)

    megapath_nano.log.print('end')


def step_short_alignment_removal(megapath_nano, short_alignment_removal):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()
    megapath_nano.log.print('length of alignment >= {threshold}'.format(threshold=short_alignment_removal.I.min_align_length))

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    #best_align_list = short_alignment_removal.I.align_list.sort_values(['read_id', 'alignment_score', 'alignment_score_tiebreaker']).drop_duplicates(subset=['read_id'], keep='last')
    best_align_list = align_list_to_best_align_list(
                                                    assembly_metadata=megapath_nano.assembly_metadata,
                                                    log=megapath_nano.log,
                                                    align_list=short_alignment_removal.I.align_list,
                                                   )

    short_alignment_removal.O.num_read_before = best_align_list.shape[0]

    read_to_remove = best_align_list.query('(sequence_to - sequence_from) < @short_alignment_removal.I.min_align_length')[['read_id']]

    short_alignment_removal.O.align_list = short_alignment_removal.I.align_list.merge(
                                                                                      right=read_to_remove.set_index('read_id').assign(read_to_remove = lambda x: 1),
                                                                                      how='left', 
                                                                                      left_on='read_id', 
                                                                                      right_index = True,
                                                                                      suffixes=['', '_y'],
                                                                                      validate='m:1',
                                                                                     ).fillna(0).query('read_to_remove == 0').drop(['read_to_remove'], axis=1)

    short_alignment_removal.O.num_read_after = short_alignment_removal.O.num_read_before - read_to_remove.shape[0]
    
    megapath_nano.log.print('{num_removed} out of {num_read} reads removed'.format(
                                                                               num_removed=read_to_remove.shape[0],
                                                                               num_read=short_alignment_removal.O.num_read_before))

    megapath_nano.log.print('end')


def step_closing_spike_filter(megapath_nano, closing_spike_filter):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()

    RAM_dir_name = tempfile.mkdtemp(prefix='closing_spike_filter.', dir=megapath_nano.RAM_dir_name)

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    #best_align_list = closing_spike_filter.I.align_list.sort_values(['read_id', 'alignment_score', 'alignment_score_tiebreaker']).drop_duplicates(subset=['read_id'], keep='last')
    best_align_list = align_list_to_best_align_list(
                                                    assembly_metadata=megapath_nano.assembly_metadata,
                                                    log=megapath_nano.log,
                                                    align_list=closing_spike_filter.I.align_list,
                                                   )

    #best_align_list_with_short_alignment = closing_spike_filter.I.align_list_with_short_alignment.sort_values(['read_id', 'alignment_score', 'alignment_score_tiebreaker']).drop_duplicates(subset=['read_id'], keep='last')
    best_align_list_with_short_alignment = align_list_to_best_align_list(
                                                                         assembly_metadata=megapath_nano.assembly_metadata,
                                                                         log=megapath_nano.log,
                                                                         align_list=closing_spike_filter.I.align_list_with_short_alignment,
                                                                        )

    closing_spike_filter.O.noise_stat = align_list_to_align_stat_by_assembly_id(
                                                                                assembly_metadata=megapath_nano.assembly_metadata,
                                                                                log=megapath_nano.log,
                                                                                align_list=best_align_list_with_short_alignment,
                                                                               )

    closing_spike_filter.O.noise_stat = closing_spike_filter.O.noise_stat.assign(expected_max_depth = lambda x: x['adjusted_average_depth'] + closing_spike_filter.I.expected_max_depth_stdev * numpy.sqrt(x['adjusted_average_depth']))
    closing_spike_filter.O.noise_stat['expected_max_depth'] = closing_spike_filter.O.noise_stat['expected_max_depth'].astype('int')
    closing_spike_filter.O.noise_stat['expected_max_depth'] = numpy.where(closing_spike_filter.O.noise_stat['expected_max_depth'] < 1, 1, closing_spike_filter.O.noise_stat['expected_max_depth'])

    if megapath_nano.global_options['debug'] == True:
        closing_spike_filter.O.noise_stat.to_csv(path_or_buf=file_prefix_with_path + '.closing_spike_noise_stat', sep='\t', header=True, index=False)

    closing_spike_filter.O.closing_spike_noise_bed, noise_span_bp = align_list_to_depth_bed(
                                                                                            temp_dir_name=RAM_dir_name,
                                                                                            align_list=best_align_list_with_short_alignment,
                                                                                            min_depth=closing_spike_filter.O.noise_stat.rename(columns={'expected_max_depth':'min_depth'}),
                                                                                            can_equal_to_min=False,
                                                                                           )

    closing_spike_filter.O.noise_stat = closing_spike_filter.O.noise_stat.merge(
                                                                                right=noise_span_bp.set_index('assembly_id').rename(columns={'span_bp': 'closing_spike_span_bp'}),
                                                                                how='inner', 
                                                                                left_on='assembly_id', 
                                                                                right_index = True,
                                                                                suffixes=['', '_y'],
                                                                                validate='1:1',
                                                                               )

    closing_spike_filter.O.noise_stat = closing_spike_filter.O.noise_stat.assign(closing_spike_span_percent = lambda x: x['closing_spike_span_bp'] / x['assembly_length'])[['assembly_id', 'closing_spike_span_bp', 'closing_spike_span_percent']]

    closing_spike_filter.O.noise_summary = closing_spike_filter.O.noise_stat['closing_spike_span_percent'].describe()


    closing_spike_filter.O.noise_bed = merge_bed_with_assembly_id(bed_list=[closing_spike_filter.I.noise_bed, closing_spike_filter.O.closing_spike_noise_bed])

    closing_spike_filter.O.num_read_before = best_align_list.shape[0]

    read_to_remove = select_alignment_by_bed(
                                             temp_dir_name=RAM_dir_name,
                                             align_list=best_align_list,
                                             bed=closing_spike_filter.O.noise_bed,
                                             min_overlap=closing_spike_filter.I.max_align_noise_overlap,
                                             can_equal_to_min=False,
                                            )[['read_id']]

    closing_spike_filter.O.align_list = closing_spike_filter.I.align_list.merge(
                                                                                right=read_to_remove.set_index('read_id').assign(read_to_remove = lambda x: 1),
                                                                                how='left', 
                                                                                left_on='read_id', 
                                                                                right_index = True,
                                                                                suffixes=['', '_y'],
                                                                                validate='m:1',
                                                                               ).fillna(0).query('read_to_remove == 0').drop(['read_to_remove'], axis=1)

    closing_spike_filter.O.num_read_after = closing_spike_filter.O.align_list.sort_values(['read_id']).drop_duplicates(subset=['read_id']).shape[0]

    if megapath_nano.global_options['debug'] == True:
        closing_spike_filter.I.align_list.to_csv(path_or_buf=file_prefix_with_path + '.before_closing_spike_filter_align_list', sep='\t', header=True, index=False)
        closing_spike_filter.O.align_list.to_csv(path_or_buf=file_prefix_with_path + '.after_closing_spike_filter_align_list', sep='\t', header=True, index=False)

    
    if megapath_nano.global_options['debug'] == False:
        shutil.rmtree(RAM_dir_name)

    megapath_nano.log.print('filter summary - average noise span: {average:.2%} - max noise span: {max:.2%}'.format(average=closing_spike_filter.O.noise_summary.loc['mean'], max=closing_spike_filter.O.noise_summary.loc['max']))
    megapath_nano.log.print('{num_removed} out of {num_read} reads removed'.format(
                                                                                    num_removed=closing_spike_filter.O.num_read_before-closing_spike_filter.O.num_read_after,
                                                                                    num_read=closing_spike_filter.O.num_read_before))
    megapath_nano.log.print('end')


def step_combine_with_human_and_decoy(megapath_nano, combine_with_human_and_decoy):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    align_list = combine_with_human_and_decoy.I.align_list.merge(
                                                                 right=combine_with_human_and_decoy.I.huamn_and_decoy_best_align_list[['read_id', 'alignment_score']].rename(columns={'alignment_score': 'human_and_decoy_best_alignment_score'}).set_index('read_id'),
                                                                 how='left', 
                                                                 left_on='read_id', 
                                                                 right_index = True,
                                                                 suffixes=['', '_y'],
                                                                 validate='m:1',
                                                                ).fillna(0).query('alignment_score > human_and_decoy_best_alignment_score').drop(['human_and_decoy_best_alignment_score'], axis=1)

    combine_with_human_and_decoy.O.align_list = pandas.concat([align_list, combine_with_human_and_decoy.I.huamn_and_decoy_best_align_list], axis=0,sort=True)

    if megapath_nano.global_options['debug'] == True:
        combine_with_human_and_decoy.O.align_list.to_csv(path_or_buf=file_prefix_with_path + '.combine_align_list', sep='\t', header=True, index=False)

    megapath_nano.log.print('end')


def step_best_alignment(megapath_nano, best_alignment):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    best_alignment.O.best_align_list = align_list_to_best_align_list(
                                                                     assembly_metadata=megapath_nano.assembly_metadata,
                                                                     log=megapath_nano.log,
                                                                     align_list=best_alignment.I.align_list,
                                                                    )

    megapath_nano.log.print('end')


def step_separate_human_and_decoy(megapath_nano, separate_human_and_decoy):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    if separate_human_and_decoy.I.human_assembly_list is not None:
        best_align_human_isin = separate_human_and_decoy.I.best_align_list['assembly_id'].isin(separate_human_and_decoy.I.human_assembly_list['assembly_id'])
        best_align_human_align = separate_human_and_decoy.I.best_align_list[best_align_human_isin]
    else:
        best_align_human_isin = None
        best_align_human_align = pandas.DataFrame(columns=align_list_col_name)
        for key, value in align_list_col_type.items():
            best_align_human_align[key] = best_align_human_align[key].astype(value)



    if separate_human_and_decoy.I.decoy_assembly_list is not None:
        best_align_decoy_isin = separate_human_and_decoy.I.best_align_list['assembly_id'].isin(separate_human_and_decoy.I.decoy_assembly_list['assembly_id'])
        best_align_decoy_align = separate_human_and_decoy.I.best_align_list[best_align_decoy_isin]
    else:
        best_align_decoy_isin = None
        best_align_decoy_align = pandas.DataFrame(columns=align_list_col_name)
        for key, value in align_list_col_type.items():
            best_align_decoy_align[key] = best_align_decoy_align[key].astype(value)

    if best_align_human_isin is not None or best_align_decoy_isin is not None:
        if best_align_human_isin is not None and best_align_decoy_isin is not None:
            best_align_microbe_isin = ~(best_align_human_isin | best_align_decoy_isin)
        elif best_align_human_isin is not None:
            best_align_microbe_isin = ~best_align_human_isin
        else:
            best_align_microbe_isin = ~best_align_decoy_isin
        separate_human_and_decoy.O.best_align_list = separate_human_and_decoy.I.best_align_list[best_align_microbe_isin]
    else:
        separate_human_and_decoy.O.best_align_list = separate_human_and_decoy.I.best_align_list

    separate_human_and_decoy.O.human_best_align_list = pandas.concat([best_align_human_align, separate_human_and_decoy.I.human_best_align_list], axis=0,sort=True)
    separate_human_and_decoy.O.decoy_best_align_list = pandas.concat([best_align_decoy_align, separate_human_and_decoy.I.decoy_best_align_list], axis=0,sort=True)

    separate_human_and_decoy.O.read_list = separate_human_and_decoy.I.read_id_list.copy()

    human_read = separate_human_and_decoy.O.human_best_align_list.assign(human_read = lambda x: 1)[['read_id', 'human_read']].set_index('read_id')
    separate_human_and_decoy.O.read_list=separate_human_and_decoy.O.read_list.drop_duplicates()
    separate_human_and_decoy.O.read_list = separate_human_and_decoy.O.read_list.merge(
                                                                                      right=human_read,
                                                                                      how='left', 
                                                                                      left_on='read_id', 
                                                                                      right_index = True,
                                                                                      suffixes=['', '_y'],
                                                                                      validate='1:1',
                                                                                     )
    separate_human_and_decoy.O.decoy_best_align_list=separate_human_and_decoy.O.decoy_best_align_list.drop_duplicates()
    decoy_read = separate_human_and_decoy.O.decoy_best_align_list.assign(decoy_read = lambda x: 1)[['read_id', 'decoy_read']].set_index('read_id')
    separate_human_and_decoy.O.read_list = separate_human_and_decoy.O.read_list.merge(
                                                                                      right=decoy_read,
                                                                                      how='left', 
                                                                                      left_on='read_id', 
                                                                                      right_index = True,
                                                                                      suffixes=['', '_y'],
                                                                                      validate='1:1',
                                                                                     )

    separate_human_and_decoy.O.best_align_list=separate_human_and_decoy.O.best_align_list.drop_duplicates()
    microbe_read = separate_human_and_decoy.O.best_align_list.assign(microbe_read = lambda x: 1)[['read_id', 'microbe_read']].set_index('read_id')
    separate_human_and_decoy.O.read_list = separate_human_and_decoy.O.read_list.merge(
                                                                                      right=microbe_read,
                                                                                      how='left', 
                                                                                      left_on='read_id', 
                                                                                      right_index = True,
                                                                                      suffixes=['', '_y'],
                                                                                      validate='1:1',
                                                                                     )

    separate_human_and_decoy.O.read_list = separate_human_and_decoy.O.read_list.fillna(0)
    separate_human_and_decoy.O.read_list['human_read'] = separate_human_and_decoy.O.read_list['human_read'].astype('int')
    separate_human_and_decoy.O.read_list['decoy_read'] = separate_human_and_decoy.O.read_list['decoy_read'].astype('int')
    separate_human_and_decoy.O.read_list['microbe_read'] = separate_human_and_decoy.O.read_list['microbe_read'].astype('int')

    separate_human_and_decoy.O.read_list['aligned'] = separate_human_and_decoy.O.read_list['human_read'] | separate_human_and_decoy.O.read_list['decoy_read'] | separate_human_and_decoy.O.read_list['microbe_read']
    separate_human_and_decoy.O.read_list = separate_human_and_decoy.O.read_list.assign(unaligned = lambda x: 1 - x['aligned'])

    separate_human_and_decoy.O.assembly_list = separate_human_and_decoy.O.best_align_list[['assembly_id']].sort_values(['assembly_id']).drop_duplicates()

    megapath_nano.log.print('end')


def step_unique_alignment(megapath_nano, unique_alignment):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()
    megapath_nano.log.print('unique alignment threshold = {threshold:.0%}'.format(threshold=unique_alignment.I.unique_align_threshold / 100))

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    best_alignment_score_assembly_id = unique_alignment.I.best_align_list[['read_id', 'assembly_id']].rename(columns={'assembly_id':'best_assembly_id'}).set_index('read_id')

    non_best_align_list = unique_alignment.I.align_list.merge(
                                                              right=best_alignment_score_assembly_id,
                                                              how='inner', 
                                                              left_on='read_id', 
                                                              right_index = True,
                                                              suffixes=['', '_y'],
                                                              validate='m:1',
                                                             ).query('assembly_id != best_assembly_id').drop(['best_assembly_id'], axis=1)

    all_align_list = pandas.concat([non_best_align_list, unique_alignment.I.human_best_align_list, unique_alignment.I.decoy_best_align_list], axis=0,sort=True)
    second_best_alignment_score = all_align_list.sort_values(['read_id', 'alignment_score']).drop_duplicates(subset=['read_id'], keep='last')
    second_best_alignment_score = second_best_alignment_score[['read_id', 'alignment_score']].rename(columns={'alignment_score':'second_best_alignment_score'}).set_index('read_id')

    best_align_list_with_second_best = unique_alignment.I.best_align_list.merge(
                                                                                right=second_best_alignment_score,
                                                                                how='left', 
                                                                                left_on='read_id', 
                                                                                right_index = True,
                                                                                suffixes=['', '_y'],
                                                                                validate='m:1',
                                                                               ).fillna(0).copy()
    best_align_list_with_second_best['second_best_alignment_score'] = best_align_list_with_second_best['second_best_alignment_score'].astype('int')

    unique_align_threshold_percent = unique_alignment.I.unique_align_threshold / 100
    unique_alignment.O.best_align_list = best_align_list_with_second_best.query('alignment_score * @unique_align_threshold_percent > second_best_alignment_score').copy()

    unique_alignment.O.num_read_before = unique_alignment.I.best_align_list.shape[0]
    unique_alignment.O.num_read_after = unique_alignment.O.best_align_list.shape[0]

    megapath_nano.log.print('{num_filtered} reads out of {num_total} removed'.format(num_filtered=unique_alignment.O.num_read_before - unique_alignment.O.num_read_after, num_total=unique_alignment.O.num_read_before))
    megapath_nano.log.print('end')


def step_noise_projection(megapath_nano, noise_projection):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()
    megapath_nano.log.print('Number of genus to perform noise projection = {num_genus}'.format(num_genus=noise_projection.I.num_genus))
    megapath_nano.log.print('Minimum relative abundance to perform noise projection = {relative_abundance:.0%}'.format(relative_abundance=noise_projection.I.min_percent_abundance / 100))
    megapath_nano.log.print('Read length bin size for simulated reads = {bin_size}'.format(bin_size=noise_projection.I.read_length_bin_size))
    megapath_nano.log.print('Read length multiplier for maximum read length over average read length = {multiplier}'.format(multiplier=noise_projection.I.read_length_multiplier))
    megapath_nano.log.print('Error profile for simulated reads = ' + noise_projection.I.error_profile)
    megapath_nano.log.print('Number of simulated reads to generate = {num_read}'.format(num_read=noise_projection.I.num_read_to_simulate))
    

    RAM_dir_name = tempfile.mkdtemp(prefix='noise_projection.', dir=megapath_nano.RAM_dir_name)

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    temp_align_list_filename = os.path.join(RAM_dir_name, 'temp_align_list_filename')

    best_align_stat = align_list_to_align_stat_by_assembly_id(
                                                              assembly_metadata=megapath_nano.assembly_metadata,
                                                              log=megapath_nano.log,
                                                              align_list=noise_projection.I.best_align_list,
                                                              noise_bed=noise_projection.I.noise_bed,
                                                             ).query('adjusted_total_aligned_bp > 0')
    print(best_align_stat)

    genus_with_more_than_1_species = best_align_stat[['genus_tax_id']].assign(species_count = lambda x: 1).groupby(['genus_tax_id'], as_index=False).sum()
    genus_with_more_than_1_species = genus_with_more_than_1_species.query('species_count > 1').drop(['species_count'], axis=1).copy()
    print(genus_with_more_than_1_species)
    best_align_stat = best_align_stat.merge(
                                            right=genus_with_more_than_1_species.set_index('genus_tax_id'),
                                            how='inner', 
                                            left_on='genus_tax_id', 
                                            right_index = True,
                                            suffixes=['', '_y'],
                                            validate='m:1',
                                           ).sort_values(['adjusted_total_aligned_bp'], ascending=[False])
    print(best_align_stat[['assembly_id', 'genus_tax_id', 'genus_height']])
    genus_list = best_align_stat.query('genus_height <= @noise_projection.I.genus_height')[['genus_tax_id', 'assembly_id', 'adjusted_total_aligned_bp']]
    print(genus_list)
    genus_list = genus_list.sort_values(['adjusted_total_aligned_bp'], ascending=[False]).drop_duplicates(subset=['genus_tax_id'], keep='first')[['genus_tax_id']].copy()

    if noise_projection.I.num_genus > genus_list.shape[0]:
        num_genus = genus_list.shape[0]
    else:
        num_genus = noise_projection.I.num_genus

    print(genus_list)

    noise_projection.O.noise_aligned_bp = pandas.DataFrame(columns=noise_projection_aligned_bp_name)
    for key, value in noise_projection_aligned_bp_type.items():
        noise_projection.O.noise_aligned_bp[key] = noise_projection.O.noise_aligned_bp[key].astype(value)

    temp_pipe_filename = os.path.join(RAM_dir_name, 'temp_pipe')
    os.mkfifo(temp_pipe_filename)

    for i in range(num_genus):

        genus_tax_id = genus_list.iloc[i]['genus_tax_id']
        print(genus_tax_id)

        genus_align_stat = best_align_stat.query('genus_tax_id == @genus_tax_id')
        print(genus_align_stat)
        highest_abundance = genus_align_stat.iloc[0]['adjusted_total_aligned_bp']
        print(highest_abundance)
        min_abundance = int(highest_abundance * noise_projection.I.min_percent_abundance / 100)
        print(min_abundance)

        print(genus_align_stat[['assembly_id', 'adjusted_total_aligned_bp', 'average_read_length']])

        noise_source_assembly_list = genus_align_stat.query('adjusted_total_aligned_bp >= @min_abundance')[['assembly_id', 'adjusted_total_aligned_bp', 'average_read_length']]
        print(noise_source_assembly_list)
        noise_source_assembly_list = megapath_nano.assembly_metadata.get_assembly_path(assembly_list=noise_source_assembly_list)
        noise_source_assembly_list['path'] = noise_source_assembly_list['path'].map(lambda x: os.path.join(megapath_nano.global_options['assembly_folder'], x))

        for j in range(noise_source_assembly_list['assembly_id'].shape[0]):

            noise_source_assembly_id = noise_source_assembly_list.iloc[j]['assembly_id']

            max_read_length = int((noise_source_assembly_list.iloc[j]['average_read_length'] * noise_projection.I.read_length_multiplier + noise_projection.I.read_length_bin_size / 2) / noise_projection.I.read_length_bin_size) * noise_projection.I.read_length_bin_size
            if max_read_length < noise_projection.I.read_length_bin_size:
                max_read_length = noise_projection.I.read_length_bin_size
            if max_read_length < noise_projection.I.min_read_length:
                max_read_length = noise_projection.I.min_read_length

            sim_read_filename = os.path.join(RAM_dir_name, noise_source_assembly_id + '.' + noise_projection.I.error_profile + '.' + str(noise_projection.I.min_read_length) + '-' + str(max_read_length) + '.' + str(noise_projection.I.num_read_to_simulate))
            
            megapath_nano.log.print(str(genus_tax_id) + '-' + noise_source_assembly_id + ' : ' + sim_read_filename)

            sim_command = [os.path.join(megapath_nano.global_options['tool_folder'], megapath_nano.global_options['read_simulator']), 'linear', '-r', temp_pipe_filename, '-c', os.path.join(megapath_nano.global_options['tool_folder'], os.path.join(megapath_nano.global_options['read_simulation_profiles'], noise_projection.I.error_profile)), '-n', str(noise_projection.I.num_read_to_simulate), '--min_len', str(noise_projection.I.min_read_length), '--max_len', str(max_read_length), '-o', sim_read_filename]
            megapath_nano.log.print(sim_command)
            sim_process = subprocess.Popen(sim_command, close_fds=True)

            temp_pipe = os.open(temp_pipe_filename, flags=os.O_WRONLY)

            zcat_command = ['zcat', '-f', noise_source_assembly_list.iloc[j]['path']]
            megapath_nano.log.print(zcat_command)
            zcat_process = subprocess.Popen(zcat_command, close_fds=True, stdout=temp_pipe)

            zcat_process.wait()
            os.close(temp_pipe)
            sim_process.wait()

            try:
                os.remove(sim_read_filename + '_error_profile')
            except(FileNotFoundError):
                pass
            try:
                os.remove(sim_read_filename + '.log')
            except(FileNotFoundError):
                pass

            if zcat_process.returncode != 0 or sim_process.returncode != 0:
                continue
            
            print(genus_align_stat[['assembly_id']])

            sim_align_list = Align(
                                   assembly_metadata=megapath_nano.assembly_metadata,
                                   global_options=megapath_nano.global_options,
                                   temp_dir_name=RAM_dir_name,
                                   log_file=megapath_nano.aligner_log,
                                   query_filename_list=pandas.DataFrame(data={'path': [sim_read_filename + '_reads.fasta',]}),
                                   target_assembly_list=genus_align_stat,
                                   aligner_options=shlex.split(megapath_nano.global_options['alignerThreadOption'] + ' -x map-ont -N 1000 -p 0'),
                                   mapping_only=megapath_nano.global_options['mapping_only'],
                                  )

            os.remove(sim_read_filename + '_reads.fasta')

            sim_best_alignment_score = sim_align_list.sort_values(['read_id', 'alignment_score']).drop_duplicates(subset='read_id', keep='last').rename(columns={'alignment_score': 'best_alignment_score'})
            sim_weighted_align_list = sim_align_list.merge(
                                                           right=sim_best_alignment_score.set_index('read_id'),
                                                           how='left', 
                                                           left_on='read_id', 
                                                           right_index = True,
                                                           suffixes=['', '_y'],
                                                           validate='m:1',
                                                          )
            sim_weighted_align_list['weighted_read_count'] = sim_weighted_align_list['alignment_score'] / sim_weighted_align_list['best_alignment_score']
            sim_weighted_align_stat = sim_weighted_align_list.groupby(['assembly_id'], as_index=False).sum()

            sim_best_align_list = align_list_to_best_align_list(
                                                                assembly_metadata=megapath_nano.assembly_metadata,
                                                                log=megapath_nano.log,
                                                                align_list=sim_align_list,
                                                               )
            sim_best_align_stat = align_list_to_align_stat_by_assembly_id(
                                                                          assembly_metadata=megapath_nano.assembly_metadata,
                                                                          log=megapath_nano.log,
                                                                          align_list=sim_best_align_list,
                                                                         )
            print(sim_best_align_stat[['assembly_id', 'adjusted_total_aligned_bp']])
            print(sim_weighted_align_stat[['assembly_id', 'weighted_read_count', 'alignment_score']])



    megapath_nano.log.print('end')


def step_similar_species_marker(megapath_nano, similar_species_marker):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()
    megapath_nano.log.print('Number of genus to check for similarity = {num_genus}'.format(num_genus=similar_species_marker.I.num_genus))
    if similar_species_marker.I.alignment_similarity_1 == similar_species_marker.I.alignment_similarity_2:
        megapath_nano.log.print('A species will be considered to be highly similar if:')
        megapath_nano.log.print('- Alignment similarity cutoff = {similarity:.0%} and aligned region >= {aligned_region:.0%} of the whole genome'.format(similarity=similar_species_marker.I.alignment_similarity_1 / 100, aligned_region=min(similar_species_marker.I.aligned_region_threshold_1, similar_species_marker.I.aligned_region_threshold_2) / 100))
    else:
        if similar_species_marker.I.similarity_combine_logic == 'and':
            megapath_nano.log.print('A species will be considered to be highly similar if both of the conditions below are satisfied:')
        else:
            megapath_nano.log.print('A species will be considered to be highly similar if either or both of the conditions below are satisfied:')
        megapath_nano.log.print('- Alignment similarity cutoff = {similarity:.0%} and aligned region >= {aligned_region:.0%} of the whole genome'.format(similarity=similar_species_marker.I.alignment_similarity_1 / 100, aligned_region=similar_species_marker.I.aligned_region_threshold_1 / 100))
        megapath_nano.log.print('- Alignment similarity cutoff = {similarity:.0%} and aligned region >= {aligned_region:.0%} of the whole genome'.format(similarity=similar_species_marker.I.alignment_similarity_2 / 100, aligned_region=similar_species_marker.I.aligned_region_threshold_2 / 100))
    
    RAM_dir_name = tempfile.mkdtemp(prefix='similar_species_marker.', dir=megapath_nano.RAM_dir_name)

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    temp_align_list_filename = os.path.join(RAM_dir_name, 'temp_align_list_filename')

    best_align_stat = align_list_to_align_stat_by_assembly_id(
                                                              assembly_metadata=megapath_nano.assembly_metadata,
                                                              log=megapath_nano.log,
                                                              align_list=similar_species_marker.I.best_align_list,
                                                              noise_bed=similar_species_marker.I.noise_bed,
                                                             )

    genus_with_more_than_1_species = best_align_stat[['genus_tax_id']].assign(species_count = lambda x: 1).groupby(['genus_tax_id'], as_index=False).sum()
    genus_with_more_than_1_species = genus_with_more_than_1_species.query('species_count > 1').drop(['species_count'], axis=1).copy()
    best_align_stat = best_align_stat.merge(
                                            right=genus_with_more_than_1_species.set_index('genus_tax_id'),
                                            how='inner', 
                                            left_on='genus_tax_id', 
                                            right_index = True,
                                            suffixes=['', '_y'],
                                            validate='m:1',
                                           )
    genus_list = best_align_stat.query('genus_height <= @similar_species_marker.I.genus_height')[['genus_tax_id', 'assembly_id', 'adjusted_total_aligned_bp']]
    genus_list = genus_list.sort_values(['adjusted_total_aligned_bp'], ascending=[False]).drop_duplicates(subset=['genus_tax_id'], keep='first')[['genus_tax_id']].copy()

    if similar_species_marker.I.num_genus > genus_list.shape[0]:
        num_genus = genus_list.shape[0]
    else:
        num_genus = similar_species_marker.I.num_genus

    similar_species_marker.O.similar_species_marker = pandas.DataFrame(columns=similar_species_marker_name_with_assembly_id)
    for key, value in similar_species_marker_type_with_assembly_id.items():
        similar_species_marker.O.similar_species_marker[key] = similar_species_marker.O.similar_species_marker[key].astype(value)

    temp_pipe_filename = os.path.join(RAM_dir_name, 'temp_pipe')
    os.mkfifo(temp_pipe_filename)

    is_similar_tag = 1

    for i in range(num_genus):

        genus_tax_id = genus_list.iloc[i]['genus_tax_id']

        genus_align_stat = best_align_stat.query('genus_tax_id == @genus_tax_id')[['assembly_id', 'genus_tax_id', 'adjusted_total_aligned_bp']].sort_values(['adjusted_total_aligned_bp'], ascending=[False])

        if similar_species_marker.I.alignment_similarity_1 == 80:
            aligner_options_1 = megapath_nano.global_options['alignerThreadOption']  + ' -N 1000 -p 0 '+ similarity_option(divergence=20)
        elif similar_species_marker.I.alignment_similarity_1 == 90:
            aligner_options_1 = megapath_nano.global_options['alignerThreadOption']  + ' -N 1000 -p 0 '+ similarity_option(divergence=10)
        elif similar_species_marker.I.alignment_similarity_1 == 95:
            aligner_options_1 = megapath_nano.global_options['alignerThreadOption']  + ' -N 1000 -p 0 '+ similarity_option(divergence=5)
        elif similar_species_marker.I.alignment_similarity_1 == 98:
            aligner_options_1 = megapath_nano.global_options['alignerThreadOption']  + ' -N 1000 -p 0 '+ similarity_option(divergence=2)
        elif similar_species_marker.I.alignment_similarity_1 == 99:
            aligner_options_1 = megapath_nano.global_options['alignerThreadOption']  + ' -N 1000 -p 0 '+ similarity_option(divergence=1)
        else:
            raise RuntimeError('alignment_similarity_1 is invalid')

        if similar_species_marker.I.alignment_similarity_2 == 80:
            aligner_options_2 = megapath_nano.global_options['alignerThreadOption']  + ' -N 1000 -p 0 '+ similarity_option(divergence=20)
        elif similar_species_marker.I.alignment_similarity_2 == 90:
            aligner_options_2 = megapath_nano.global_options['alignerThreadOption']  + ' -N 1000 -p 0 '+ similarity_option(divergence=10)
        elif similar_species_marker.I.alignment_similarity_2 == 95:
            aligner_options_2 = megapath_nano.global_options['alignerThreadOption']  + ' -N 1000 -p 0 '+ similarity_option(divergence=5)
        elif similar_species_marker.I.alignment_similarity_2 == 98:
            aligner_options_2 = megapath_nano.global_options['alignerThreadOption']  + ' -N 1000 -p 0 '+ similarity_option(divergence=2)
        elif similar_species_marker.I.alignment_similarity_2 == 99:
            aligner_options_2 = megapath_nano.global_options['alignerThreadOption']  + ' -N 1000 -p 0 '+ similarity_option(divergence=1)
        else:
            raise RuntimeError('alignment_similarity_2 is invalid')

        align_list_1 = Align(
                             assembly_metadata=megapath_nano.assembly_metadata,
                             global_options=megapath_nano.global_options,
                             temp_dir_name=RAM_dir_name,
                             log_file=megapath_nano.aligner_log,
                             query_assembly_list=genus_align_stat.iloc[0:1],
                             target_assembly_list=genus_align_stat.iloc[1:],
                             aligner_options=shlex.split(aligner_options_1),
                             mapping_only=megapath_nano.global_options['mapping_only'],
                            )
        if align_list_1 is not None and align_list_1.empty == False:
            covered_bed_1 = align_list_to_bed(align_list=align_list_1)
            covered_bp_1 = bed_to_covered_bp_by_assembly_id(bed=covered_bed_1)
            covered_bp_1 = megapath_nano.assembly_metadata.get_assembly_length(assembly_list=covered_bp_1)
            covered_percent_1 = covered_bp_1.assign(covered_percent = lambda x: x['covered_bp'] / x['assembly_length'])[['assembly_id', 'covered_percent']]

        if aligner_options_2 != aligner_options_1:
            align_list_2 = Align(
                                 assembly_metadata=megapath_nano.assembly_metadata,
                                 global_options=megapath_nano.global_options,
                                 temp_dir_name=RAM_dir_name,
                                 log_file=megapath_nano.aligner_log,
                                 query_assembly_list=genus_align_stat.iloc[0:1],
                                 target_assembly_list=genus_align_stat.iloc[1:],
                                 aligner_options=shlex.split(aligner_options_2),
                                 mapping_only=megapath_nano.global_options['mapping_only'],
                                )
            covered_bed_2 = align_list_to_bed(align_list=align_list_2)
            covered_bp_2 = bed_to_covered_bp_by_assembly_id(bed=covered_bed_2)
            covered_bp_2 = megapath_nano.assembly_metadata.get_assembly_length(assembly_list=covered_bp_2)
            covered_percent_2 = covered_bp_2.assign(covered_percent = lambda x: x['covered_bp'] / x['assembly_length'])[['assembly_id', 'covered_percent']]
        else:
            covered_percent_2 = covered_percent_1

        assembly_id_passed_condition_1 = covered_percent_1.query('covered_percent >= @similar_species_marker.I.aligned_region_threshold_1 / 100')[['assembly_id']]
        assembly_id_passed_condition_2 = covered_percent_2.query('covered_percent >= @similar_species_marker.I.aligned_region_threshold_2 / 100')[['assembly_id']]

        if similar_species_marker.I.similarity_combine_logic == 'and':
            assembly_id_passed_condition = assembly_id_passed_condition_1.merge(
                                                                                right=assembly_id_passed_condition_2.set_index('assembly_id'),
                                                                                how='inner', 
                                                                                left_on='assembly_id', 
                                                                                right_index = True,
                                                                                suffixes=['', '_y'],
                                                                                validate='1:1',
                                                                               )
        else:
            assembly_id_passed_condition = pandas.concat([assembly_id_passed_condition_1, assembly_id_passed_condition_2]).sort_values(['assembly_id']).drop_duplicates()

        marker = genus_align_stat[['assembly_id', 'genus_tax_id']].rename(columns={'genus_tax_id': 'similar_genus_tax_id'})
        marker = marker.merge(
                              right=covered_percent_1.rename(columns={'covered_percent': 'aligned_percent_1'}).set_index('assembly_id'),
                              how='left', 
                              left_on='assembly_id', 
                              right_index = True,
                              suffixes=['', '_y'],
                              validate='1:1',
                             ).fillna(0)
        marker = marker.merge(
                              right=covered_percent_2.rename(columns={'covered_percent': 'aligned_percent_2'}).set_index('assembly_id'),
                              how='left', 
                              left_on='assembly_id', 
                              right_index = True,
                              suffixes=['', '_y'],
                              validate='1:1',
                             ).fillna(0)

        if assembly_id_passed_condition.empty == False:
            marker = marker.merge(
                                  right=assembly_id_passed_condition.assign(is_similar = lambda x: 1).set_index('assembly_id'),
                                  how='left', 
                                  left_on='assembly_id', 
                                  right_index = True,
                                  suffixes=['', '_y'],
                                  validate='1:1',
                                 ).fillna(0)
            assembly_id = genus_align_stat.iloc[0]['assembly_id']
            marker.loc[marker['assembly_id'] == assembly_id, 'is_similar'] = 1
            marker.loc[marker['assembly_id'] == assembly_id, 'aligned_percent_1'] = 100
            marker.loc[marker['assembly_id'] == assembly_id, 'aligned_percent_2'] = 100

            marker['is_similar'] = marker['is_similar'] * is_similar_tag
            is_similar_tag += 1
        
        else:
            marker = marker.assign(is_similar = lambda x: 0)

        marker['similarity_cutoff_1'] = similar_species_marker.I.alignment_similarity_1
        marker['similarity_cutoff_2'] = similar_species_marker.I.alignment_similarity_2
        marker['combine_logic'] = similar_species_marker.I.similarity_combine_logic

        for key, value in similar_species_marker_type.items():
            marker[key] = marker[key].astype(value)

        similar_species_marker.O.similar_species_marker = pandas.concat([similar_species_marker.O.similar_species_marker, marker],sort=True)


    megapath_nano.log.print('end')


def step_noise_detection_statistics(megapath_nano, noise_detection_statistics):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    noise_detection_statistics.O.noise_stat = noise_detection_statistics.I.assembly_list.merge(
                                                                                               right=noise_detection_statistics.I.spike_noise_stat.set_index('assembly_id'),
                                                                                               how='left', 
                                                                                               left_on='assembly_id', 
                                                                                               right_index = True,
                                                                                               suffixes=['', '_y'],
                                                                                               validate='1:1',
                                                                                              )
    noise_detection_statistics.O.noise_stat = noise_detection_statistics.O.noise_stat.merge(
                                                                                            right=noise_detection_statistics.I.human_noise_stat.set_index('assembly_id'),
                                                                                            how='left', 
                                                                                            left_on='assembly_id', 
                                                                                            right_index = True,
                                                                                            suffixes=['', '_y'],
                                                                                            validate='1:1',
                                                                                           )
    noise_detection_statistics.O.noise_stat = noise_detection_statistics.O.noise_stat.merge(
                                                                                            right=noise_detection_statistics.I.microbe_noise_stat.set_index('assembly_id'),
                                                                                            how='left', 
                                                                                            left_on='assembly_id', 
                                                                                            right_index = True,
                                                                                            suffixes=['', '_y'],
                                                                                            validate='1:1',
                                                                                           )
    noise_detection_statistics.O.noise_stat = noise_detection_statistics.O.noise_stat.merge(
                                                                                            right=noise_detection_statistics.I.closing_spike_noise_stat.set_index('assembly_id'),
                                                                                            how='left', 
                                                                                            left_on='assembly_id', 
                                                                                            right_index = True,
                                                                                            suffixes=['', '_y'],
                                                                                            validate='1:1',
                                                                                           )

    covered_bp = bed_to_covered_bp_by_assembly_id(bed=noise_detection_statistics.I.noise_bed)

    noise_detection_statistics.O.noise_stat = noise_detection_statistics.O.noise_stat.drop(['noise_span_bp', 'noise_span_percent'], axis=1).merge(
                                                                                                                                                  right=covered_bp.rename(columns={'covered_bp': 'noise_span_bp'}).set_index('assembly_id'),
                                                                                                                                                  how='left', 
                                                                                                                                                  left_on='assembly_id', 
                                                                                                                                                  right_index = True,
                                                                                                                                                  suffixes=['', '_y'],
                                                                                                                                                  validate='1:1',
                                                                                                                                                 )
    noise_detection_statistics.O.noise_stat = noise_detection_statistics.O.noise_stat.assign(noise_span_percent = lambda x: x['noise_span_bp'] / x['assembly_length'])

    noise_detection_statistics.O.noise_stat = noise_detection_statistics.O.noise_stat.fillna(0)

    for key, value in noise_detection_stat_col_type.items():
        if value == int:
            noise_detection_statistics.O.noise_stat[key] = noise_detection_statistics.O.noise_stat[key].astype(value)
    
    noise_detection_statistics.O.noise_stat = noise_detection_statistics.O.noise_stat[noise_detection_stat_col_name_by_assembly_id]

    noise_detection_statistics.O.noise_summary = noise_detection_statistics.O.noise_stat['noise_span_percent'].describe()

    megapath_nano.log.print('Overall noise summary - average noise span: {average:.2%} - max noise span: {max:.2%}'.format(average=noise_detection_statistics.O.noise_summary.loc['mean'], max=noise_detection_statistics.O.noise_summary.loc['max']))
    megapath_nano.log.print('end')


def step_noise_removal_statistics(megapath_nano, noise_removal_statistics):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()

    RAM_dir_name = tempfile.mkdtemp(prefix='noise_removal_statistics.', dir=megapath_nano.RAM_dir_name)

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    read_id = noise_removal_statistics.I.best_align_list[['read_id', 'read_length']]

    spike_noise_align_list = select_alignment_by_bed(
                                                       temp_dir_name=RAM_dir_name,
                                                       align_list=noise_removal_statistics.I.best_align_list,
                                                       bed=noise_removal_statistics.I.spike_noise_bed,
                                                       min_overlap=noise_removal_statistics.I.max_align_noise_overlap,
                                                       can_equal_to_min=False,
                                                      )
    human_noise_align_list = select_alignment_by_bed(
                                                       temp_dir_name=RAM_dir_name,
                                                       align_list=noise_removal_statistics.I.best_align_list,
                                                       bed=noise_removal_statistics.I.human_noise_bed,
                                                       min_overlap=noise_removal_statistics.I.max_align_noise_overlap,
                                                       can_equal_to_min=False,
                                                      )
    microbe_noise_align_list = select_alignment_by_bed(
                                                         temp_dir_name=RAM_dir_name,
                                                         align_list=noise_removal_statistics.I.best_align_list,
                                                         bed=noise_removal_statistics.I.microbe_noise_bed,
                                                         min_overlap=noise_removal_statistics.I.max_align_noise_overlap,
                                                         can_equal_to_min=False,
                                                        )
    closing_spike_noise_align_list = select_alignment_by_bed(
                                                             temp_dir_name=RAM_dir_name,
                                                             align_list=noise_removal_statistics.I.best_align_list,
                                                             bed=noise_removal_statistics.I.closing_spike_noise_bed,
                                                             min_overlap=noise_removal_statistics.I.max_align_noise_overlap,
                                                             can_equal_to_min=False,
                                                            )
    noise_align_list = select_alignment_by_bed(
                                               temp_dir_name=RAM_dir_name,
                                               align_list=noise_removal_statistics.I.best_align_list,
                                               bed=noise_removal_statistics.I.noise_bed,
                                               min_overlap=noise_removal_statistics.I.max_align_noise_overlap,
                                               can_equal_to_min=False,
                                              )
    short_align_list = noise_removal_statistics.I.best_align_list.query('sequence_to - sequence_from < @noise_removal_statistics.I.min_align_length')
    all_align_list = pandas.concat([noise_align_list, short_align_list], axis=0,sort=True).sort_values(['read_id']).drop_duplicates(subset=['read_id'])

    spike_noise_align_stat = spike_noise_align_list.assign(spike_total_number_of_read = lambda x: 1, spike_total_read_bp = lambda x: x['read_length'], spike_total_aligned_bp = lambda x: x['sequence_to'] - x['sequence_from'])[['assembly_id', 'spike_total_number_of_read', 'spike_total_read_bp', 'spike_total_aligned_bp']].groupby(['assembly_id'], as_index=False).sum().copy()
    human_noise_align_stat = human_noise_align_list.assign(human_total_number_of_read = lambda x: 1, human_total_read_bp = lambda x: x['read_length'], human_total_aligned_bp = lambda x: x['sequence_to'] - x['sequence_from'])[['assembly_id', 'human_total_number_of_read', 'human_total_read_bp', 'human_total_aligned_bp']].groupby(['assembly_id'], as_index=False).sum().copy()
    microbe_noise_align_stat = microbe_noise_align_list.assign(microbe_total_number_of_read = lambda x: 1, microbe_total_read_bp = lambda x: x['read_length'], microbe_total_aligned_bp = lambda x: x['sequence_to'] - x['sequence_from'])[['assembly_id', 'microbe_total_number_of_read', 'microbe_total_read_bp', 'microbe_total_aligned_bp']].groupby(['assembly_id'], as_index=False).sum().copy()
    closing_spike_noise_align_stat = closing_spike_noise_align_list.assign(closing_spike_total_number_of_read = lambda x: 1, closing_spike_total_read_bp = lambda x: x['read_length'], closing_spike_total_aligned_bp = lambda x: x['sequence_to'] - x['sequence_from'])[['assembly_id', 'closing_spike_total_number_of_read', 'closing_spike_total_read_bp', 'closing_spike_total_aligned_bp']].groupby(['assembly_id'], as_index=False).sum().copy()
    noise_align_stat = noise_align_list.assign(noise_total_number_of_read = lambda x: 1, noise_total_read_bp = lambda x: x['read_length'], noise_total_aligned_bp = lambda x: x['sequence_to'] - x['sequence_from'])[['assembly_id', 'noise_total_number_of_read', 'noise_total_read_bp', 'noise_total_aligned_bp']].groupby(['assembly_id'], as_index=False).sum().copy()
    short_align_stat = short_align_list.assign(short_total_number_of_read = lambda x: 1, short_total_read_bp = lambda x: x['read_length'], short_total_aligned_bp = lambda x: x['sequence_to'] - x['sequence_from'])[['assembly_id', 'short_total_number_of_read', 'short_total_read_bp', 'short_total_aligned_bp']].groupby(['assembly_id'], as_index=False).sum().copy()
    all_align_stat = all_align_list.assign(all_total_number_of_read = lambda x: 1, all_total_read_bp = lambda x: x['read_length'], all_total_aligned_bp = lambda x: x['sequence_to'] - x['sequence_from'])[['assembly_id', 'all_total_number_of_read', 'all_total_read_bp', 'all_total_aligned_bp']].groupby(['assembly_id'], as_index=False).sum().copy()

    noise_removal_statistics.O.noise_best_align_list = all_align_list

    noise_removal_statistics.O.noise_align_stat = noise_removal_statistics.I.assembly_list
    noise_removal_statistics.O.noise_align_stat = noise_removal_statistics.O.noise_align_stat.merge(
                                                                                                    right=spike_noise_align_stat.set_index('assembly_id'),
                                                                                                    how='left', 
                                                                                                    left_on='assembly_id', 
                                                                                                    right_index = True,
                                                                                                    suffixes=['', '_y'],
                                                                                                    validate='1:1',
                                                                                                   )
    noise_removal_statistics.O.noise_align_stat = noise_removal_statistics.O.noise_align_stat.merge(
                                                                                                    right=human_noise_align_stat.set_index('assembly_id'),
                                                                                                    how='left', 
                                                                                                    left_on='assembly_id', 
                                                                                                    right_index = True,
                                                                                                    suffixes=['', '_y'],
                                                                                                    validate='1:1',
                                                                                                   )
    noise_removal_statistics.O.noise_align_stat = noise_removal_statistics.O.noise_align_stat.merge(
                                                                                                    right=microbe_noise_align_stat.set_index('assembly_id'),
                                                                                                    how='left', 
                                                                                                    left_on='assembly_id', 
                                                                                                    right_index = True,
                                                                                                    suffixes=['', '_y'],
                                                                                                    validate='1:1',
                                                                                                   )
    noise_removal_statistics.O.noise_align_stat = noise_removal_statistics.O.noise_align_stat.merge(
                                                                                                    right=closing_spike_noise_align_stat.set_index('assembly_id'),
                                                                                                    how='left', 
                                                                                                    left_on='assembly_id', 
                                                                                                    right_index = True,
                                                                                                    suffixes=['', '_y'],
                                                                                                    validate='1:1',
                                                                                                   )
    noise_removal_statistics.O.noise_align_stat = noise_removal_statistics.O.noise_align_stat.merge(
                                                                                                    right=noise_align_stat.set_index('assembly_id'),
                                                                                                    how='left', 
                                                                                                    left_on='assembly_id', 
                                                                                                    right_index = True,
                                                                                                    suffixes=['', '_y'],
                                                                                                    validate='1:1',
                                                                                                   )
    noise_removal_statistics.O.noise_align_stat = noise_removal_statistics.O.noise_align_stat.merge(
                                                                                                    right=short_align_stat.set_index('assembly_id'),
                                                                                                    how='left', 
                                                                                                    left_on='assembly_id', 
                                                                                                    right_index = True,
                                                                                                    suffixes=['', '_y'],
                                                                                                    validate='1:1',
                                                                                                   )
    noise_removal_statistics.O.noise_align_stat = noise_removal_statistics.O.noise_align_stat.merge(
                                                                                                    right=all_align_stat.set_index('assembly_id'),
                                                                                                    how='left', 
                                                                                                    left_on='assembly_id', 
                                                                                                    right_index = True,
                                                                                                    suffixes=['', '_y'],
                                                                                                    validate='1:1',
                                                                                                   )
    noise_removal_statistics.O.noise_align_stat = noise_removal_statistics.O.noise_align_stat.fillna(0)
    for key, value in noise_removal_stat_col_type.items():
        if value == int:
            noise_removal_statistics.O.noise_align_stat[key] = noise_removal_statistics.O.noise_align_stat[key].astype(value)


    if megapath_nano.global_options['debug'] == False:
        shutil.rmtree(RAM_dir_name)

    megapath_nano.log.print('end')


def step_noise_source_statistics(megapath_nano, noise_source_statistics):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    microbe_noise_align_list = noise_source_statistics.I.noise_best_align_list.rename(columns={'assembly_id': 'target_assembly_id'}).merge(
                                                                                                                                           right=noise_source_statistics.I.best_align_list[['read_id', 'assembly_id']].rename(columns={'assembly_id': 'source_assembly_id'}).set_index('read_id'),
                                                                                                                                           how='inner', 
                                                                                                                                           left_on='read_id', 
                                                                                                                                           right_index = True,
                                                                                                                                           suffixes=['', '_y'],
                                                                                                                                           validate='1:1',
                                                                                                                                          )

    if megapath_nano.global_options['debug'] == True:
        noise_source_statistics.I.noise_best_align_list.to_csv(path_or_buf=file_prefix_with_path + '.noise_source_noise_best', sep='\t', header=True, index=False)
        noise_source_statistics.I.best_align_list.to_csv(path_or_buf=file_prefix_with_path + '.noise_source_best', sep='\t', header=True, index=False)
        microbe_noise_align_list.to_csv(path_or_buf=file_prefix_with_path + '.noise_source_list', sep='\t', header=True, index=False)

    microbe_noise_align_list = microbe_noise_align_list.assign(noise_read_count = lambda x: 1, noise_aligned_bp = lambda x: x['sequence_to'] - x['sequence_from'], noise_read_bp = lambda x: x['read_length'])
    microbe_noise_align_stat = microbe_noise_align_list[['target_assembly_id', 'source_assembly_id', 'noise_read_count', 'noise_aligned_bp', 'noise_read_bp']].groupby(['target_assembly_id', 'source_assembly_id'], as_index=False).sum().copy()

    non_microbe_noise_align_list = noise_source_statistics.I.noise_best_align_list.rename(columns={'assembly_id': 'target_assembly_id'}).merge(
                                                                                                                                               right=noise_source_statistics.I.read_list.drop(['read_length'], axis=1).query('microbe_read == 0').set_index('read_id'),
                                                                                                                                               how='inner', 
                                                                                                                                               left_on='read_id', 
                                                                                                                                               right_index = True,
                                                                                                                                               suffixes=['', '_y'],
                                                                                                                                               validate='1:1',
                                                                                                                                              )

    non_microbe_noise_align_list = non_microbe_noise_align_list.assign(noise_read_count = lambda x: 1, noise_aligned_bp = lambda x: x['sequence_to'] - x['sequence_from'], noise_read_bp = lambda x: x['read_length'])
    non_microbe_noise_align_list['source_assembly_id'] = numpy.where(non_microbe_noise_align_list['aligned'] == 0, 'unidentified', numpy.where(non_microbe_noise_align_list['human_read'] == 1, 'human', numpy.where(non_microbe_noise_align_list['decoy_read'] == 1, 'decoy', 'error')))
    non_microbe_noise_align_stat = non_microbe_noise_align_list[['target_assembly_id', 'source_assembly_id', 'noise_read_count', 'noise_aligned_bp', 'noise_read_bp']].groupby(['target_assembly_id', 'source_assembly_id'], as_index=False).sum().copy()

    noise_source_statistics.O.noise_source_stat = pandas.concat([microbe_noise_align_stat, non_microbe_noise_align_stat], axis=0).query('target_assembly_id != source_assembly_id').copy()

    if megapath_nano.global_options['debug'] == True:
        noise_source_statistics.O.noise_source_stat.to_csv(path_or_buf=file_prefix_with_path + '.noise_source_stat', sep='\t', header=True, index=False)

    megapath_nano.log.print('end')


def step_max_adjusted_abundance(megapath_nano, max_adjusted_abundance):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    combined_align_list = pandas.concat([max_adjusted_abundance.I.noise_best_align_list, max_adjusted_abundance.I.best_align_list])

    combined_align_stat = align_list_to_align_stat_by_assembly_id(
                                                                  assembly_metadata=megapath_nano.assembly_metadata,
                                                                  log=megapath_nano.log,
                                                                  align_list=combined_align_list,
                                                                 )

    max_adjusted_abundance.O.abundance_stat = combined_align_stat[['assembly_id', 'adjusted_total_aligned_bp', 'adjusted_average_depth', 'adjusted_covered_percent']].rename(columns={
                                                                                                                                                                                      'adjusted_total_aligned_bp': 'max_adjusted_total_aligned_bp',
                                                                                                                                                                                      'adjusted_average_depth': 'max_adjusted_average_depth',
                                                                                                                                                                                      'adjusted_covered_percent': 'max_adjusted_covered_percent',
                                                                                                                                                                                     })

    if megapath_nano.global_options['debug'] == True:
        max_adjusted_abundance.O.abundance_stat.to_csv(path_or_buf=file_prefix_with_path + '.max_adjusted_abundance', sep='\t', header=True, index=False)

    megapath_nano.log.print('end')


def step_read_statistics(megapath_nano, read_statistics):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    read_statistics.O.read_list = read_statistics.I.read_info.merge(
                                                                    right=read_statistics.I.read_list.drop(['read_length'], axis=1).set_index('read_id'),
                                                                    how='left', 
                                                                    left_on='read_id', 
                                                                    right_index = True,
                                                                    suffixes=['', '_y'],
                                                                    validate='1:1',
                                                                   ).fillna(0)
    for key, value in aligned_read_list_col_type.items():
        if value == int:
            read_statistics.O.read_list[key] = read_statistics.O.read_list[key].astype(value)


    read_statistics.O.read_list['aligned'] = read_statistics.O.read_list['aligned'].astype('int')
    read_statistics.O.read_list['human_read'] = read_statistics.O.read_list['human_read'].astype('int')
    read_statistics.O.read_list['decoy_read'] = read_statistics.O.read_list['decoy_read'].astype('int')
    read_statistics.O.read_list['microbe_read'] = read_statistics.O.read_list['microbe_read'].astype('int')

    read_statistics.O.read_stat = read_statistics.O.read_list.rename(columns={'read_length': 'total_read_bp'}).assign(
                                                                                                                      dummy = lambda x: 1,
                                                                                                                      total_number_of_read = lambda x: 1,
                                                                                                                      total_passed_filter_read_bp = lambda x: x['total_read_bp'] * x['passed_filter'],
                                                                                                                      total_aligned_read_bp = lambda x: x['total_read_bp'] * x['aligned'],
                                                                                                                      total_unaligned_read_bp = lambda x: x['total_read_bp'] * x['unaligned'],
                                                                                                                      total_human_read_bp = lambda x: x['total_read_bp'] * x['human_read'],
                                                                                                                      total_decoy_read_bp = lambda x: x['total_read_bp'] * x['decoy_read'],
                                                                                                                      total_microbe_read_bp = lambda x: x['total_read_bp'] * x['microbe_read'],
                                                                                                                     ).groupby(['dummy'], as_index=False).sum().drop(['dummy'], axis=1).copy()

    read_statistics.O.read_list['quality_score_bin'] = numpy.floor(read_statistics.O.read_list['average_quality'] / read_statistics.I.quality_score_bin_size) * read_statistics.I.quality_score_bin_size
    read_statistics.O.read_list['read_length_bin'] = numpy.floor(read_statistics.O.read_list['read_length'] / read_statistics.I.read_length_bin_size) * read_statistics.I.read_length_bin_size
    read_statistics.O.read_list['read_length_bin'] = read_statistics.O.read_list['read_length_bin'].astype('int')

    read_statistics.O.quality_score_bin = read_statistics.O.read_list[['quality_score_bin']].assign(read_count = lambda x: 1).groupby(['quality_score_bin'], as_index=False).sum().copy()
    read_statistics.O.read_length_bin = read_statistics.O.read_list[['read_length_bin']].assign(read_count = lambda x: 1).groupby(['read_length_bin'], as_index=False).sum().copy()

    read_statistics.O.passed_quality_score_bin = read_statistics.O.read_list.query('passed_filter == 1')[['quality_score_bin']].assign(read_count = lambda x: 1).groupby(['quality_score_bin'], as_index=False).sum().copy()
    read_statistics.O.passed_read_length_bin = read_statistics.O.read_list.query('passed_filter == 1')[['read_length_bin']].assign(read_count = lambda x: 1).groupby(['read_length_bin'], as_index=False).sum().copy()

    read_statistics.O.human_quality_score_bin = read_statistics.O.read_list.query('human_read == 1')[['quality_score_bin']].assign(read_count = lambda x: 1).groupby(['quality_score_bin'], as_index=False).sum().copy()
    read_statistics.O.human_read_length_bin = read_statistics.O.read_list.query('human_read == 1')[['read_length_bin']].assign(read_count = lambda x: 1).groupby(['read_length_bin'], as_index=False).sum().copy()

    read_statistics.O.decoy_quality_score_bin = read_statistics.O.read_list.query('decoy_read == 1')[['quality_score_bin']].assign(read_count = lambda x: 1).groupby(['quality_score_bin'], as_index=False).sum().copy()
    read_statistics.O.decoy_read_length_bin = read_statistics.O.read_list.query('decoy_read == 1')[['read_length_bin']].assign(read_count = lambda x: 1).groupby(['read_length_bin'], as_index=False).sum().copy()

    read_statistics.O.microbe_quality_score_bin = read_statistics.O.read_list.query('microbe_read == 1')[['quality_score_bin']].assign(read_count = lambda x: 1).groupby(['quality_score_bin'], as_index=False).sum().copy()
    read_statistics.O.microbe_read_length_bin = read_statistics.O.read_list.query('microbe_read == 1')[['read_length_bin']].assign(read_count = lambda x: 1).groupby(['read_length_bin'], as_index=False).sum().copy()

    read_statistics.O.aligned_quality_score_bin = read_statistics.O.read_list.query('aligned == 1')[['quality_score_bin']].assign(read_count = lambda x: 1).groupby(['quality_score_bin'], as_index=False).sum().copy()
    read_statistics.O.aligned_read_length_bin = read_statistics.O.read_list.query('aligned == 1')[['read_length_bin']].assign(read_count = lambda x: 1).groupby(['read_length_bin'], as_index=False).sum().copy()

    read_statistics.O.unaligned_quality_score_bin = read_statistics.O.read_list.query('aligned == 0')[['quality_score_bin']].assign(read_count = lambda x: 1).groupby(['quality_score_bin'], as_index=False).sum().copy()
    read_statistics.O.unaligned_read_length_bin = read_statistics.O.read_list.query('aligned == 0')[['read_length_bin']].assign(read_count = lambda x: 1).groupby(['read_length_bin'], as_index=False).sum().copy()


    megapath_nano.log.print('end')


def step_format_output(megapath_nano, options):

    megapath_nano.log.print('start')
    megapath_nano.log.print_time()

    file_prefix_with_path = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix)

    # Genome sets used

    megapath_nano.assembly_info = get_assembly_info(
                                               assembly_metadata=megapath_nano.assembly_metadata,
                                               db_conn=megapath_nano.taxonomy_db_conn,
                                               log=megapath_nano.log,
                                               assembly_list=pandas.concat([
                                                                            megapath_nano.human_assembly_list[['assembly_id']],
                                                                            megapath_nano.decoy_assembly_list[['assembly_id']],
                                                                            megapath_nano.species_id_assembly_list[['assembly_id']],
                                                                            megapath_nano.assembly_id_assembly_list[['assembly_id']],
                                                                            pandas.DataFrame(data={'assembly_id': [megapath_nano.global_options['human_repetitive_region_filter_assembly_id'],]})
                                                                           ], axis=0).sort_values(['assembly_id']).drop_duplicates()
                                              )
    megapath_nano.assembly_info = megapath_nano.assembly_info.merge(
                                                          right=megapath_nano.human_assembly_list.assign(human = lambda x: 1).set_index('assembly_id'),
                                                          how='left', 
                                                          left_on='assembly_id', 
                                                          right_index = True,
                                                          suffixes=['', '_y'],
                                                          validate='1:1',
                                                          )
    megapath_nano.assembly_info = megapath_nano.assembly_info.merge(
                                                          right=megapath_nano.decoy_assembly_list.assign(decoy = lambda x: 1).set_index('assembly_id'),
                                                          how='left', 
                                                          left_on='assembly_id', 
                                                          right_index = True,
                                                          suffixes=['', '_y'],
                                                          validate='1:1',
                                                         )
    megapath_nano.assembly_info = megapath_nano.assembly_info.merge(
                                                          right=megapath_nano.species_id_assembly_list.assign(species = lambda x: 1).set_index('assembly_id'),
                                                          how='left', 
                                                          left_on='assembly_id', 
                                                          right_index = True,
                                                          suffixes=['', '_y'],
                                                          validate='1:1',
                                                         )
    megapath_nano.assembly_info = megapath_nano.assembly_info.merge(
                                                          right=megapath_nano.assembly_id_assembly_list.assign(assembly = lambda x: 1).set_index('assembly_id'),
                                                          how='left', 
                                                          left_on='assembly_id', 
                                                          right_index = True,
                                                          suffixes=['', '_y'],
                                                          validate='1:1',
                                                         )
    megapath_nano.assembly_info = megapath_nano.assembly_info.fillna(0)
    megapath_nano.assembly_info['human'] = megapath_nano.assembly_info['human'].astype('int')
    megapath_nano.assembly_info['decoy'] = megapath_nano.assembly_info['decoy'].astype('int')
    megapath_nano.assembly_info['species'] = megapath_nano.assembly_info['species'].astype('int')
    megapath_nano.assembly_info['assembly'] = megapath_nano.assembly_info['assembly'].astype('int')

    if options.output_genome_set == True: 
        megapath_nano.assembly_info.sort_values(['assembly_id']).to_csv(path_or_buf=file_prefix_with_path + '.genome_set',
                                                                   sep='\t', header=True, index=False,
                                                                   columns=['human', 'decoy', 'species', 'assembly',] + assembly_info_col_name)

    # Per-read data

    if options.output_per_read_data == True:
        temp_per_read_dir_name = tempfile.mkdtemp(prefix='per_read.', dir=megapath_nano.temp_dir_name)
        temp_per_read_prefix_with_path = os.path.join(temp_per_read_dir_name, megapath_nano.output_prefix)

        megapath_nano.read_list.sort_values(['read_id']).to_csv(path_or_buf=temp_per_read_prefix_with_path + '.read_list', sep='\t', header=True, index=False, columns=consolidated_read_list_col_name)

        megapath_nano.human_best_align_list.sort_values(['read_id']).to_csv(path_or_buf=temp_per_read_prefix_with_path + '.human_list', sep='\t', header=True, index=False, columns=align_list_col_name)
        megapath_nano.decoy_best_align_list.sort_values(['read_id']).to_csv(path_or_buf=temp_per_read_prefix_with_path + '.decoy_list', sep='\t', header=True, index=False, columns=align_list_col_name)
        if options.output_id_signal == True:
            megapath_nano.id_best_align_list.sort_values(['read_id']).to_csv(path_or_buf=temp_per_read_prefix_with_path + '.id_list', sep='\t', header=True, index=False, columns=align_list_col_name)
        if options.output_raw_signal == True:
            megapath_nano.raw_best_align_list.sort_values(['read_id']).to_csv(path_or_buf=temp_per_read_prefix_with_path + '.raw_list', sep='\t', header=True, index=False, columns=align_list_col_name)
        megapath_nano.align_list.sort_values(['read_id']).to_csv(path_or_buf=temp_per_read_prefix_with_path + '.list', sep='\t', header=True, index=False, columns=align_list_col_name)
        megapath_nano.best_align_list.sort_values(['read_id']).to_csv(path_or_buf=temp_per_read_prefix_with_path + '.microbe_list', sep='\t', header=True, index=False, columns=align_list_col_name)
        if options.unique_alignment == True:
            megapath_nano.unique_align_list.sort_values(['read_id']).to_csv(path_or_buf=temp_per_read_prefix_with_path + '.unique_list', sep='\t', header=True, index=False, columns=align_list_col_name)
        if options.output_noise_stat == True:
            megapath_nano.noise_removal_best_align_list.sort_values(['read_id']).to_csv(path_or_buf=temp_per_read_prefix_with_path + '.noise_list', sep='\t', header=True, index=False, columns=align_list_col_name)

        shutil.make_archive(file_prefix_with_path + '.per_read', format=options.archive_format, root_dir=temp_per_read_dir_name)
        shutil.rmtree(temp_per_read_dir_name)

    # read statistics

    if options.output_quality_score_histogram == True:
        temp_quality_score_dir_name = tempfile.mkdtemp(prefix='quality_score.', dir=megapath_nano.temp_dir_name)
        temp_quality_score_prefix_with_path = os.path.join(temp_quality_score_dir_name, megapath_nano.output_prefix)

        megapath_nano.quality_score_bin.to_csv(path_or_buf=temp_quality_score_prefix_with_path + '.quality_score_histogram', sep='\t', header=True, index=False, columns=quality_score_histogram_col_name, float_format='%.2f')
        megapath_nano.passed_quality_score_bin.to_csv(path_or_buf=temp_quality_score_prefix_with_path + '.passed_quality_score_histogram', sep='\t', header=True, index=False, columns=quality_score_histogram_col_name, float_format='%.2f')
        megapath_nano.human_quality_score_bin.to_csv(path_or_buf=temp_quality_score_prefix_with_path + '.human_quality_score_histogram', sep='\t', header=True, index=False, columns=quality_score_histogram_col_name, float_format='%.2f')
        megapath_nano.decoy_quality_score_bin.to_csv(path_or_buf=temp_quality_score_prefix_with_path + '.decoy_quality_score_histogram', sep='\t', header=True, index=False, columns=quality_score_histogram_col_name, float_format='%.2f')
        megapath_nano.microbe_quality_score_bin.to_csv(path_or_buf=temp_quality_score_prefix_with_path + '.microbe_quality_score_histogram', sep='\t', header=True, index=False, columns=quality_score_histogram_col_name, float_format='%.2f')
        megapath_nano.aligned_quality_score_bin.to_csv(path_or_buf=temp_quality_score_prefix_with_path + '.aligned_quality_score_histogram', sep='\t', header=True, index=False, columns=quality_score_histogram_col_name, float_format='%.2f')
        megapath_nano.unaligned_quality_score_bin.to_csv(path_or_buf=temp_quality_score_prefix_with_path + '.unaligned_quality_score_histogram', sep='\t', header=True, index=False, columns=quality_score_histogram_col_name, float_format='%.2f')

        shutil.make_archive(file_prefix_with_path + '.quality_score', format=options.archive_format, root_dir=temp_quality_score_dir_name)
        shutil.rmtree(temp_quality_score_dir_name)

    if options.output_read_length_histogram == True:
        temp_read_length_dir_name = tempfile.mkdtemp(prefix='quality_score.', dir=megapath_nano.temp_dir_name)
        temp_read_length_prefix_with_path = os.path.join(temp_read_length_dir_name, megapath_nano.output_prefix)

        megapath_nano.read_length_bin.to_csv(path_or_buf=temp_read_length_prefix_with_path + '.read_length_histogram', sep='\t', header=True, index=False, columns=read_length_histogram_col_name)
        megapath_nano.passed_read_length_bin.to_csv(path_or_buf=temp_read_length_prefix_with_path + '.passed_read_length_histogram', sep='\t', header=True, index=False, columns=read_length_histogram_col_name)
        megapath_nano.human_read_length_bin.to_csv(path_or_buf=temp_read_length_prefix_with_path + '.human_read_length_histogram', sep='\t', header=True, index=False, columns=read_length_histogram_col_name)
        megapath_nano.decoy_read_length_bin.to_csv(path_or_buf=temp_read_length_prefix_with_path + '.decoy_read_length_histogram', sep='\t', header=True, index=False, columns=read_length_histogram_col_name)
        megapath_nano.microbe_read_length_bin.to_csv(path_or_buf=temp_read_length_prefix_with_path + '.microbe_read_length_histogram', sep='\t', header=True, index=False, columns=read_length_histogram_col_name)
        megapath_nano.aligned_read_length_bin.to_csv(path_or_buf=temp_read_length_prefix_with_path + '.aligned_read_length_histogram', sep='\t', header=True, index=False, columns=read_length_histogram_col_name)
        megapath_nano.unaligned_read_length_bin.to_csv(path_or_buf=temp_read_length_prefix_with_path + '.unaligned_read_length_histogram', sep='\t', header=True, index=False, columns=read_length_histogram_col_name)

        shutil.make_archive(file_prefix_with_path + '.read_length', format=options.archive_format, root_dir=temp_read_length_dir_name)
        shutil.rmtree(temp_read_length_dir_name)

    megapath_nano.read_stat.to_csv(path_or_buf=file_prefix_with_path + '.read_stat', sep='\t', header=True, index=False, columns=read_stat_col_name)


    # Analysis
    if not megapath_nano.human_assembly_list.empty:
        megapath_nano.human_sequence_list = megapath_nano.assembly_metadata.get_sequence_tax_id(assembly_list=megapath_nano.human_assembly_list)
        megapath_nano.human_sequence_name = get_sequence_name(db_conn=megapath_nano.taxonomy_db_conn, log=megapath_nano.log, sequence_list=megapath_nano.human_sequence_list)
        megapath_nano.human_best_align_stat = generate_align_stat_group_by_sequence_id(
                                                                                  assembly_metadata=megapath_nano.assembly_metadata,
                                                                                  log=megapath_nano.log,
                                                                                  align_list=megapath_nano.human_best_align_list,
                                                                                  sequence_info=megapath_nano.human_sequence_name,
                                                                                 )
    if not megapath_nano.decoy_assembly_list.empty:
        megapath_nano.decoy_sequence_list = megapath_nano.assembly_metadata.get_sequence_tax_id(assembly_list=megapath_nano.decoy_assembly_list)
        megapath_nano.decoy_sequence_name = get_sequence_name(db_conn=megapath_nano.taxonomy_db_conn, log=megapath_nano.log, sequence_list=megapath_nano.decoy_sequence_list)
        megapath_nano.decoy_best_align_stat = generate_align_stat_group_by_sequence_id(
                                                                                  assembly_metadata=megapath_nano.assembly_metadata,
                                                                                  log=megapath_nano.log,
                                                                                  align_list=megapath_nano.decoy_best_align_list,
                                                                                  sequence_info=megapath_nano.decoy_sequence_name,
                                                                                 )
    
    #test generate_align_stat_by_sequence_id
    #sequence_list = megapath_nano.assembly_metadata.get_sequence_tax_id(assembly_list=megapath_nano.assembly_list)
    #sequence_name = get_sequence_name(db_conn=megapath_nano.taxonomy_db_conn, log=megapath_nano.log, sequence_list=sequence_list)
    #best_align_stat = generate_align_stat_group_by_sequence_id(
    #                                                                    assembly_metadata=megapath_nano.assembly_metadata,
    #                                                                    log=megapath_nano.log,
    #                                                                    align_list=megapath_nano.best_align_list,
    #                                                                    sequence_info=sequence_name,
    #                                                                    noise_bed=megapath_nano.noise_bed,
    #                                                                    )
    #sequence_list.to_csv(path_or_buf=file_prefix_with_path + '.sequence_list', sep='\t', header=True, index=False)
    #sequence_name.to_csv(path_or_buf=file_prefix_with_path + '.sequence_name', sep='\t', header=True, index=False)
    #best_align_stat.to_csv(path_or_buf=file_prefix_with_path + '.microbe_stat_by_sequence_id', sep='\t', header=True, index=False)
    #end of testing

    if FLAGS.output_id_signal == True:
        megapath_nano.id_best_align_stat = generate_align_stat_group_by_assembly_id(
                                                                               assembly_metadata=megapath_nano.assembly_metadata,
                                                                               log=megapath_nano.log,
                                                                               align_list=megapath_nano.id_best_align_list,
                                                                               assembly_info=megapath_nano.assembly_info,
                                                                              )

    megapath_nano.raw_best_align_stat = generate_align_stat_group_by_assembly_id(
                                                                            assembly_metadata=megapath_nano.assembly_metadata,
                                                                            log=megapath_nano.log,
                                                                            align_list=megapath_nano.raw_best_align_list,
                                                                            assembly_info=megapath_nano.assembly_info,
                                                                           )

    megapath_nano.best_align_stat = generate_align_stat_group_by_assembly_id(
                                                                        assembly_metadata=megapath_nano.assembly_metadata,
                                                                        log=megapath_nano.log,
                                                                        align_list=megapath_nano.best_align_list,
                                                                        assembly_info=megapath_nano.assembly_info,
                                                                        noise_bed=megapath_nano.noise_bed,
                                                                       )
                                                                       

    raw_best_align_stat = megapath_nano.raw_best_align_stat.rename(columns={
                                                                       'total_number_of_read': 'pre_total_number_of_read',
                                                                       'total_read_bp': 'pre_total_read_bp',
                                                                       'average_read_length': 'pre_average_read_length',
                                                                       'total_aligned_bp': 'pre_total_aligned_bp',
                                                                       'average_depth': 'pre_average_depth',
                                                                       'covered_percent': 'pre_covered_percent',
                                                                       'average_identity': 'pre_average_identity',
                                                                       'average_edit_dist': 'pre_average_edit_dist',
                                                                       'average_alignment_score': 'pre_average_alignment_score',
                                                                      })
    raw_no_best_align_stat = raw_best_align_stat.merge(
                                                       right=megapath_nano.best_align_stat[['assembly_id']].assign(joined = lambda x: 1).set_index('assembly_id'),
                                                       how='left', 
                                                       left_on='assembly_id', 
                                                       right_index=True,
                                                       suffixes=['', '_y'],
                                                       validate='1:1',
                                                      ).fillna(0).query('joined == 0').drop(['joined'], axis=1)
    raw_no_best_align_stat = raw_no_best_align_stat[[
                                                     'assembly_id',
                                                     'assembly_length',
                                                     'tax_id',
                                                     'tax_name',
                                                     'species_tax_id',
                                                     'species_tax_name',
                                                     'genus_tax_id',
                                                     'genus_tax_name',
                                                     'pre_total_number_of_read',
                                                     'pre_total_read_bp',
                                                     'pre_average_read_length',
                                                     'pre_total_aligned_bp',
                                                     'pre_average_depth',
                                                     'pre_covered_percent',
                                                     'pre_average_identity',
                                                     'pre_average_edit_dist',
                                                     'pre_average_alignment_score',
                                                    ]]

    megapath_nano.best_align_stat = megapath_nano.best_align_stat.merge(
                                                              right=raw_best_align_stat[[
                                                                                         'assembly_id',
                                                                                         'pre_total_number_of_read',
                                                                                         'pre_total_read_bp',
                                                                                         'pre_average_read_length',
                                                                                         'pre_total_aligned_bp',
                                                                                         'pre_average_depth',
                                                                                         'pre_covered_percent',
                                                                                         'pre_average_identity',
                                                                                         'pre_average_edit_dist',
                                                                                         'pre_average_alignment_score',
                                                                                        ]].set_index('assembly_id'),
                                                              how='left', 
                                                              left_on='assembly_id',
                                                              right_index=True,
                                                              suffixes=['', '_y'],
                                                              validate='1:1',
                                                             ).fillna(0)

    num_assembly_id = megapath_nano.best_align_stat.shape[0]
    megapath_nano.best_align_stat = megapath_nano.best_align_stat.merge(
                                                              right=megapath_nano.max_adjusted_abundance.set_index('assembly_id'),
                                                              how='inner', 
                                                              left_on='assembly_id', 
                                                              right_index=True,
                                                              suffixes=['', '_y'],
                                                              validate='1:1',
                                                             )
    if megapath_nano.best_align_stat.shape[0] != num_assembly_id:
        raise RuntimeError('Assembly ID unmatched')

    megapath_nano.best_align_stat['adjusted_total_aligned_bp'] = megapath_nano.best_align_stat[['adjusted_total_aligned_bp', 'max_adjusted_total_aligned_bp']].min(axis=1)
    megapath_nano.best_align_stat['adjusted_average_depth'] = megapath_nano.best_align_stat[['adjusted_average_depth', 'max_adjusted_average_depth']].min(axis=1)
    megapath_nano.best_align_stat['adjusted_covered_percent'] = megapath_nano.best_align_stat[['adjusted_covered_percent', 'max_adjusted_covered_percent']].min(axis=1)
    megapath_nano.best_align_stat = megapath_nano.best_align_stat.drop(['max_adjusted_total_aligned_bp', 'max_adjusted_average_depth', 'max_adjusted_covered_percent'], axis=1)

    megapath_nano.best_align_stat = pandas.concat([megapath_nano.best_align_stat, raw_no_best_align_stat], axis=0,sort=True).fillna(0)

    megapath_nano.best_align_stat = megapath_nano.best_align_stat.merge(
                                                              right=megapath_nano.similar_species_marker.set_index('assembly_id'),
                                                              how='left', 
                                                              left_on='assembly_id', 
                                                              right_index=True,
                                                              suffixes=['', '_y'],
                                                              validate='1:1',
                                                             ).fillna(0)

    for key, value in align_stat_col_type_with_pre_noise_with_similar_species_marker_by_assembly_id.items():
        if value == int:
            megapath_nano.best_align_stat[key] = megapath_nano.best_align_stat[key].astype(value)


    if options.unique_alignment == True:
        megapath_nano.unique_align_stat = generate_align_stat_group_by_assembly_id(
                                                                              assembly_metadata=megapath_nano.assembly_metadata,
                                                                              log=megapath_nano.log,
                                                                              align_list=megapath_nano.unique_align_list,
                                                                              assembly_info=megapath_nano.assembly_info,
                                                                              noise_bed=megapath_nano.noise_bed,
                                                                             )
        raw_no_unique_align_stat = raw_best_align_stat.merge(
                                                             right=megapath_nano.unique_align_stat[['assembly_id']].assign(joined = lambda x: 1).set_index('assembly_id'),
                                                             how='left', 
                                                             left_on='assembly_id', 
                                                             right_index=True,
                                                             suffixes=['', '_y'],
                                                             validate='1:1',
                                                            ).fillna(0).query('joined == 0').drop(['joined'], axis=1)
        raw_no_unique_align_stat = raw_no_unique_align_stat[[
                                                             'assembly_id',
                                                             'assembly_length',
                                                             'tax_id',
                                                             'tax_name',
                                                             'species_tax_id',
                                                             'species_tax_name',
                                                             'genus_tax_id',
                                                             'genus_tax_name',
                                                             'pre_total_number_of_read',
                                                             'pre_total_read_bp',
                                                             'pre_average_read_length',
                                                             'pre_total_aligned_bp',
                                                             'pre_average_depth',
                                                             'pre_covered_percent',
                                                             'pre_average_identity',
                                                             'pre_average_edit_dist',
                                                             'pre_average_alignment_score',
                                                            ]]

        megapath_nano.unique_align_stat = megapath_nano.unique_align_stat.merge(
                                                                      right=raw_best_align_stat[[
                                                                                                 'assembly_id',
                                                                                                 'pre_total_number_of_read',
                                                                                                 'pre_total_read_bp',
                                                                                                 'pre_average_read_length',
                                                                                                 'pre_total_aligned_bp',
                                                                                                 'pre_average_depth',
                                                                                                 'pre_covered_percent',
                                                                                                 'pre_average_identity',
                                                                                                 'pre_average_edit_dist',
                                                                                                 'pre_average_alignment_score',
                                                                                                ]].set_index('assembly_id'),
                                                                      how='left', 
                                                                      left_on='assembly_id',
                                                                      right_index=True,
                                                                      suffixes=['', '_y'],
                                                                      validate='1:1',
                                                                     ).fillna(0)

        num_assembly_id = megapath_nano.unique_align_stat.shape[0]
        megapath_nano.unique_align_stat = megapath_nano.unique_align_stat.merge(
                                                                      right=megapath_nano.max_adjusted_abundance.set_index('assembly_id'),
                                                                      how='inner', 
                                                                      left_on='assembly_id', 
                                                                      right_index=True,
                                                                      suffixes=['', '_y'],
                                                                      validate='1:1',
                                                                     )
        if megapath_nano.unique_align_stat.shape[0] != num_assembly_id:
            raise RuntimeError('Assembly ID unmatched')

        megapath_nano.unique_align_stat['adjusted_total_aligned_bp'] = megapath_nano.unique_align_stat[['adjusted_total_aligned_bp', 'max_adjusted_total_aligned_bp']].min(axis=1)
        megapath_nano.unique_align_stat['adjusted_average_depth'] = megapath_nano.unique_align_stat[['adjusted_average_depth', 'max_adjusted_average_depth']].min(axis=1)
        megapath_nano.unique_align_stat['adjusted_covered_percent'] = megapath_nano.unique_align_stat[['adjusted_covered_percent', 'max_adjusted_covered_percent']].min(axis=1)
        megapath_nano.unique_align_stat = megapath_nano.unique_align_stat.drop(['max_adjusted_total_aligned_bp', 'max_adjusted_average_depth', 'max_adjusted_covered_percent'], axis=1)

        megapath_nano.unique_align_stat = pandas.concat([megapath_nano.unique_align_stat, raw_no_unique_align_stat], axis=0,sort=True).fillna(0)

        megapath_nano.unique_align_stat = megapath_nano.unique_align_stat.merge(
                                                                      right=megapath_nano.similar_species_marker.set_index('assembly_id'),
                                                                      how='left', 
                                                                      left_on='assembly_id', 
                                                                      right_index=True,
                                                                      suffixes=['', '_y'],
                                                                      validate='1:1',
                                                                     ).fillna(0)

        for key, value in align_stat_col_type_with_pre_noise_with_similar_species_marker_by_assembly_id.items():
            if value == int:
                megapath_nano.unique_align_stat[key] = megapath_nano.unique_align_stat[key].astype(value)



    if options.output_human_stat == True and not megapath_nano.human_assembly_list.empty:
        megapath_nano.human_best_align_stat.sort_values(['sequence_length'], ascending=[False]).to_csv(path_or_buf=file_prefix_with_path + '.human_stat', sep='\t', header=True, index=False, columns=align_stat_raw_col_name_by_sequence_id)
    if options.output_decoy_stat == True and not megapath_nano.decoy_assembly_list.empty:
        megapath_nano.decoy_best_align_stat.sort_values(['total_aligned_bp'], ascending=[False]).to_csv(path_or_buf=file_prefix_with_path + '.decoy_stat', sep='\t', header=True, index=False, columns=align_stat_raw_col_name_by_sequence_id)
    if options.output_id_signal == True:
        megapath_nano.id_best_align_stat.sort_values(['adjusted_total_aligned_bp'], ascending=[False]).to_csv(path_or_buf=file_prefix_with_path + '.id_stat', sep='\t', header=True, index=False, columns=align_stat_raw_col_name_by_assembly_id)
    if options.output_raw_signal == True:
        megapath_nano.raw_best_align_stat.sort_values(['adjusted_total_aligned_bp'], ascending=[False]).to_csv(path_or_buf=file_prefix_with_path + '.raw_stat', sep='\t', header=True, index=False, columns=align_stat_raw_col_name_by_assembly_id)
    megapath_nano.best_align_stat = megapath_nano.best_align_stat.sort_values(['adjusted_total_aligned_bp', 'pre_total_aligned_bp'], ascending=[False, False])
    megapath_nano.best_align_stat.to_csv(path_or_buf=file_prefix_with_path + '.microbe_stat', sep='\t', header=True, index=False, columns=align_stat_col_name_with_pre_noise_with_similar_species_marker_by_assembly_id)
    if options.unique_alignment == True:
        megapath_nano.unique_align_stat.sort_values(['adjusted_total_aligned_bp', 'pre_total_aligned_bp'], ascending=[False, False]).to_csv(path_or_buf=file_prefix_with_path + '.unique_stat', sep='\t', header=True, index=False, columns=align_stat_col_name_with_pre_noise_with_similar_species_marker_by_assembly_id)


    # for downstream processing 

    megapath_nano.best_align_stat.query('adjusted_total_aligned_bp > 0').sort_values(['adjusted_total_aligned_bp'], ascending=[False]).to_csv(path_or_buf=file_prefix_with_path + '.preport', sep='\t', header=True, index=False, columns=align_stat_col_name_dot_report)
    
    if FLAGS.reassignment==False:
        taxon_df=pandas.read_csv(f'{megapath_nano.global_options["db_folder"]}/sequence_name',sep='\t',header=None,names=['sequence_id','name'])
        if FLAGS.resolution=='species':
            taxon_df['name']=taxon_df['name'].apply(lambda x: " ".join(x.split(" ",2)[0:2]) if ' sp. ' not in x else " ".join(x.split(" ",3)[0:3]))
        align_list_species_name=megapath_nano.id_best_align_list.merge(right=taxon_df,on=['sequence_id'],how='left')
        #update name
        align_list_species_name=pandas.concat([align_list_species_name[align_list_species_name['name'].isnull()].assign(name=lambda x: x['sequence_id']),align_list_species_name[~align_list_species_name['name'].isnull()]])
    elif FLAGS.reassignment==True:
        align_list_species_name=megapath_nano.id_best_align_list
    sequence_id_groupby_count=align_list_species_name.groupby(['name']).count()['read_id'].sort_values(ascending=False)
    sequence_id_groupby_count.to_csv(path_or_buf=f'{file_prefix_with_path}.read_count_by_name',sep='\t')

    if FLAGS.reassignment==False:
        #testing for sequence id in best_align_stat
        megapath_nano.sequence_list = megapath_nano.assembly_metadata.get_sequence_tax_id(assembly_list=megapath_nano.assembly_list)
        megapath_nano.sequence_name = get_sequence_name(db_conn=megapath_nano.taxonomy_db_conn, log=megapath_nano.log, sequence_list=megapath_nano.sequence_list)
        megapath_nano.best_align_stat_by_sequence_id = generate_align_stat_group_by_sequence_id(
                                                                                            assembly_metadata=megapath_nano.assembly_metadata,
                                                                                            log=megapath_nano.log,
                                                                                            align_list=megapath_nano.best_align_list,
                                                                                            sequence_info=megapath_nano.sequence_name,
                                                                                            )
        megapath_nano.best_align_stat_by_sequence_id.sort_values(['total_aligned_bp'], ascending=[False]).to_csv(path_or_buf=file_prefix_with_path + '.microbe_stat_by_sequence_id', sep='\t', header=True, index=False, columns=align_stat_raw_col_name_by_sequence_id) 
        megapath_nano.sequence_info = megapath_nano.sequence_list[['assembly_id', 'sequence_id', 'sequence_length']].merge(
                                                                                                                right=megapath_nano.assembly_info.set_index('assembly_id'),
                                                                                                                how='left',
                                                                                                                left_on='assembly_id',
                                                                                                                right_index=True,
                                                                                                                validate='m:1',
                                                                                                                )

        megapath_nano.best_align_stat_by_sequence_id_with_assembly_info = megapath_nano.best_align_stat_by_sequence_id.merge(
                                                                                                            right=megapath_nano.sequence_info.set_index('sequence_id'),
                                                                                                            how='left',
                                                                                                            left_on='sequence_id',
                                                                                                            right_index = True,
                                                                                                            suffixes=['', '_y'],
                                                                                                            validate='m:1',
                                                                                                            )
        megapath_nano.best_align_stat_by_sequence_id_with_assembly_info = megapath_nano.best_align_stat_by_sequence_id_with_assembly_info.merge(
                                                                                                            right = megapath_nano.best_align_stat.set_index('assembly_id'),
                                                                                                            how = 'left',
                                                                                                            left_on = 'assembly_id',
                                                                                                            right_index = True,
                                                                                                            suffixes=['', '_y'],
                                                                                                            validate='m:1',
                                                                                                            ).rename(columns={'total_aligned_bp':'sequence_total_aligned_bp', 'average_depth': 'sequence_average_depth','covered_percent': 'sequence_covered_percent', 'adjusted_total_aligned_bp': 'assembly_adjusted_total_aligned_bp', 'adjusted_average_depth': 'assembly_adjusted_average_depth', 'adjusted_covered_percent': 'assembly_adjusted_covered_percent'})
        
        megapath_nano.best_align_stat_by_sequence_id_with_assembly_info.reindex(columns = align_stat_raw_col_name_by_sequence_id_with_assembly_id).to_csv(path_or_buf=file_prefix_with_path + '.microbe_stat_by_sequence_id_assembly_info', sep='\t', header=True, index=False)
    
    #Independent megapath_nano report for plasmids
    #megapath_nano.best_align_stat.to_csv(path_or_buf=file_prefix_with_path + '.micro.info', sep='\t', header=True, index=False)
    #megapath_nano.decoy_best_align_stat.to_csv(path_or_buf=file_prefix_with_path + '.plasmid.info', sep='\t', header=True, index=False)
    #end of Independent megapath_nano report for plasmids

    #end of testing

    bam_filter = good_align_list(
                                 align_list=megapath_nano.align_list,
                                 good_align_threshold=megapath_nano.global_options['good_alignment_threshold'],
                                )

    bam_filter.sort_values(['read_id', 'sequence_id', 'sequence_from', 'sequence_to']).to_csv(path_or_buf=file_prefix_with_path + '.bam_filter_good', sep='\t', header=True, index=False, columns=align_list_col_name_bam_filter)

    megapath_nano.best_align_list.sort_values(['read_id', 'sequence_id', 'sequence_from', 'sequence_to']).to_csv(path_or_buf=file_prefix_with_path + '.bam_filter', sep='\t', header=True, index=False, columns=align_list_col_name_bam_filter)

    megapath_nano.noise_bed.moveto(os.path.join(file_prefix_with_path + '.noise.bed'))

    megapath_nano.similar_noise_bed = merge_bed_with_assembly_id(bed_list=[megapath_nano.human_repetitive_region_noise_bed, megapath_nano.microbe_repetitive_region_noise_bed])

    if options.output_separate_noise_bed == True and (options.spike_filter == True or options.human_repetitive_region_filter == True or options.microbe_repetitive_region_filter == True):
        if options.spike_filter == True:
            megapath_nano.spike_noise_bed.moveto(file_prefix_with_path + '.spike_noise.bed')
        if options.human_repetitive_region_filter == True:
            megapath_nano.human_repetitive_region_noise_bed.moveto(file_prefix_with_path + '.human_noise.bed')
        if options.microbe_repetitive_region_filter == True:
            megapath_nano.microbe_repetitive_region_noise_bed.moveto(file_prefix_with_path + '.microbe_noise.bed')
        if options.closing_spike_filter == True:
            megapath_nano.closing_spike_noise_bed.moveto(file_prefix_with_path + '.closing_spike_noise.bed')


    # Noise analysis

    if options.output_noise_stat == True:
        temp_noise_analysis_dir_name = tempfile.mkdtemp(prefix='noise_analysis.', dir=megapath_nano.temp_dir_name)
        temp_noise_analysis_prefix_with_path = os.path.join(temp_noise_analysis_dir_name, megapath_nano.output_prefix)

        sort_order = megapath_nano.best_align_stat[['assembly_id']]

        megapath_nano.noise_stat = megapath_nano.noise_detection_stat.merge(
                                                                  right=megapath_nano.noise_removal_stat.set_index('assembly_id'),
                                                                  how='outer', 
                                                                  left_on='assembly_id', 
                                                                  right_index = True,
                                                                  suffixes=['', '_y'],
                                                                  validate='1:1',
                                                                 )
     
        megapath_nano.noise_stat_with_tax_name = megapath_nano.noise_stat.merge(
                                                                      right=megapath_nano.assembly_info.set_index('assembly_id'),
                                                                      how='left', 
                                                                      left_on='assembly_id', 
                                                                      right_index = True,
                                                                      suffixes=['', '_y'],
                                                                      validate='1:1',
                                                                     )
        megapath_nano.noise_stat_with_tax_name = megapath_nano.assembly_metadata.get_assembly_length(assembly_list=megapath_nano.noise_stat_with_tax_name, how='left').fillna(0)[noise_stat_col_name_with_tax_name]

        megapath_nano.noise_stat_with_tax_name = sort_order.merge(
                                                             right=megapath_nano.noise_stat_with_tax_name.set_index('assembly_id'),
                                                             how='inner', 
                                                             left_on='assembly_id', 
                                                             right_index = True,
                                                             suffixes=['', '_y'],
                                                             validate='1:1',
                                                            )

        for key, value in noise_stat_col_type_with_tax_name.items():
            if value == int:
                megapath_nano.noise_stat_with_tax_name[key] = megapath_nano.noise_stat_with_tax_name[key].astype(value)

        megapath_nano.noise_stat_with_tax_name.to_csv(path_or_buf=temp_noise_analysis_prefix_with_path + '.noise_stat', sep='\t', header=True, index=False, columns=noise_stat_col_name_with_tax_name)


        source_assembly_info = megapath_nano.assembly_info[['assembly_id', 'species_tax_id', 'species_tax_name']].rename(columns={'species_tax_id': 'source_species_tax_id', 'species_tax_name': 'source_species_tax_name'})
        target_assembly_info = megapath_nano.assembly_info[['assembly_id', 'species_tax_id', 'species_tax_name']].rename(columns={'species_tax_id': 'target_species_tax_id', 'species_tax_name': 'target_species_tax_name'})

        megapath_nano.noise_source_stat_with_tax_name = megapath_nano.noise_source_stat.merge(
                                                                                    right=source_assembly_info.set_index('assembly_id'),
                                                                                    how='left', 
                                                                                    left_on='source_assembly_id', 
                                                                                    right_index = True,
                                                                                    suffixes=['', '_y'],
                                                                                    validate='m:1',
                                                                                   ).merge(
                                                                                           right=target_assembly_info.set_index('assembly_id'),
                                                                                           how='left', 
                                                                                           left_on='target_assembly_id', 
                                                                                           right_index = True,
                                                                                           suffixes=['', '_y'],
                                                                                           validate='m:1',
                                                                                          )
        megapath_nano.noise_source_stat_with_tax_name = megapath_nano.assembly_metadata.get_assembly_length(assembly_list=megapath_nano.noise_source_stat_with_tax_name.assign(assembly_id = lambda x: x['source_assembly_id']), how='left').fillna(0)
        megapath_nano.noise_source_stat_with_tax_name = megapath_nano.noise_source_stat_with_tax_name.rename(columns={'assembly_length': 'source_assembly_length'}).drop(['assembly_id'], axis=1)
        megapath_nano.noise_source_stat_with_tax_name = megapath_nano.assembly_metadata.get_assembly_length(assembly_list=megapath_nano.noise_source_stat_with_tax_name.assign(assembly_id = lambda x: x['target_assembly_id']), how='left').fillna(0)
        megapath_nano.noise_source_stat_with_tax_name = megapath_nano.noise_source_stat_with_tax_name.rename(columns={'assembly_length': 'target_assembly_length'}).drop(['assembly_id'], axis=1)[noise_source_stat_col_name_with_tax_name]
        for key, value in noise_source_stat_col_type_with_tax_name.items():
            if value == int:
                megapath_nano.noise_source_stat_with_tax_name[key] = megapath_nano.noise_source_stat_with_tax_name[key].astype(value)

        megapath_nano.noise_source_stat_with_tax_name['source_species_tax_name'] = numpy.where(megapath_nano.noise_source_stat_with_tax_name['source_assembly_id'] == 'human', 'Human Genome Set', megapath_nano.noise_source_stat_with_tax_name['source_species_tax_name'])
        megapath_nano.noise_source_stat_with_tax_name['source_species_tax_name'] = numpy.where(megapath_nano.noise_source_stat_with_tax_name['source_assembly_id'] == 'decoy', 'Decoy Genome Set', megapath_nano.noise_source_stat_with_tax_name['source_species_tax_name'])
        megapath_nano.noise_source_stat_with_tax_name['source_species_tax_name'] = numpy.where(megapath_nano.noise_source_stat_with_tax_name['source_assembly_id'] == 'unidentified', 'Unidentified', megapath_nano.noise_source_stat_with_tax_name['source_species_tax_name'])
        megapath_nano.noise_source_stat_with_tax_name['source_assembly_id'] = numpy.where(megapath_nano.noise_source_stat_with_tax_name['source_assembly_id'] == 'human', 'N/A', megapath_nano.noise_source_stat_with_tax_name['source_assembly_id'])
        megapath_nano.noise_source_stat_with_tax_name['source_assembly_id'] = numpy.where(megapath_nano.noise_source_stat_with_tax_name['source_assembly_id'] == 'decoy', 'N/A', megapath_nano.noise_source_stat_with_tax_name['source_assembly_id'])
        megapath_nano.noise_source_stat_with_tax_name['source_assembly_id'] = numpy.where(megapath_nano.noise_source_stat_with_tax_name['source_assembly_id'] == 'unidentified', 'N/A', megapath_nano.noise_source_stat_with_tax_name['source_assembly_id'])

        megapath_nano.noise_source_stat_with_tax_name = megapath_nano.noise_source_stat_with_tax_name.sort_values(['target_assembly_id', 'source_assembly_id'])
        megapath_nano.noise_source_stat_with_tax_name = sort_order.rename(columns={'assembly_id': 'target_assembly_id'}).merge(
                                                                                                                          right=megapath_nano.noise_source_stat_with_tax_name.set_index('target_assembly_id'),
                                                                                                                          how='inner', 
                                                                                                                          left_on='target_assembly_id', 
                                                                                                                          right_index = True,
                                                                                                                          suffixes=['', '_y'],
                                                                                                                          validate='1:m',
                                                                                                                         )
        megapath_nano.noise_source_stat_with_tax_name.to_csv(path_or_buf=temp_noise_analysis_prefix_with_path + '.noise_source', sep='\t', header=True, index=False, columns=noise_source_stat_col_name_with_tax_name)


        megapath_nano.similar_noise_span_bp = pandas.DataFrame(columns=similar_noise_span_bp_col_name)

        if options.human_repetitive_region_filter == True:
            human_repetitive_region_noise_span_bp = megapath_nano.noise_detection_stat.query('human_span_bp > 0')[['assembly_id', 'human_span_bp']].rename(columns={'assembly_id': 'target_assembly_id', 'human_span_bp': 'microbe_span_bp'})
            human_repetitive_region_noise_span_bp = human_repetitive_region_noise_span_bp.assign(source_assembly_id = lambda x: megapath_nano.global_options['human_repetitive_region_filter_assembly_id'])
            megapath_nano.similar_noise_span_bp = pandas.concat([megapath_nano.similar_noise_span_bp, human_repetitive_region_noise_span_bp], axis=0,sort=True)

        if options.microbe_repetitive_region_filter == True:
            megapath_nano.similar_noise_span_bp = pandas.concat([megapath_nano.similar_noise_span_bp, megapath_nano.source_target_noise_span_bp], axis=0,sort=True)

        if options.human_repetitive_region_filter == True or options.microbe_repetitive_region_filter == True:
            megapath_nano.similar_noise_span_bp = megapath_nano.similar_noise_span_bp.merge(
                                                                                  right=bed_to_covered_bp_by_assembly_id(bed=megapath_nano.similar_noise_bed).set_index('assembly_id'),
                                                                                  how='left', 
                                                                                  left_on='target_assembly_id', 
                                                                                  right_index = True,
                                                                                  suffixes=['', '_y'],
                                                                                  validate='m:1',
                                                                                  ).rename(columns={'covered_bp': 'total_similar_noise_bp'}).fillna(0)

            megapath_nano.similar_noise_span_bp_with_tax_name = megapath_nano.similar_noise_span_bp.merge(
                                                                                                right=source_assembly_info.set_index('assembly_id'),
                                                                                                how='left', 
                                                                                                left_on='source_assembly_id', 
                                                                                                right_index = True,
                                                                                                suffixes=['', '_y'],
                                                                                                validate='m:1',
                                                                                               ).merge(
                                                                                                       right=target_assembly_info.set_index('assembly_id'),
                                                                                                       how='left', 
                                                                                                       left_on='target_assembly_id', 
                                                                                                       right_index = True,
                                                                                                       suffixes=['', '_y'],
                                                                                                       validate='m:1',
                                                                                                      )

            megapath_nano.similar_noise_span_bp_with_tax_name = megapath_nano.assembly_metadata.get_assembly_length(assembly_list=megapath_nano.similar_noise_span_bp_with_tax_name.assign(assembly_id = lambda x: x['source_assembly_id']), how='left').fillna(0)
            megapath_nano.similar_noise_span_bp_with_tax_name = megapath_nano.similar_noise_span_bp_with_tax_name.rename(columns={'assembly_length': 'source_assembly_length'}).drop(['assembly_id'], axis=1)
            megapath_nano.similar_noise_span_bp_with_tax_name = megapath_nano.assembly_metadata.get_assembly_length(assembly_list=megapath_nano.similar_noise_span_bp_with_tax_name.assign(assembly_id = lambda x: x['target_assembly_id']), how='left').fillna(0)
            megapath_nano.similar_noise_span_bp_with_tax_name = megapath_nano.similar_noise_span_bp_with_tax_name.rename(columns={'assembly_length': 'target_assembly_length'}).drop(['assembly_id'], axis=1)
            megapath_nano.similar_noise_span_bp_with_tax_name = megapath_nano.similar_noise_span_bp_with_tax_name.fillna(0).assign(total_similar_noise_percent = lambda x: x['total_similar_noise_bp'] / x['target_assembly_length'])[similar_noise_span_bp_col_name_with_tax_name]

            for key, value in similar_noise_span_bp_col_type_with_tax_name.items():
                if value == int:
                    megapath_nano.similar_noise_span_bp_with_tax_name[key] = megapath_nano.similar_noise_span_bp_with_tax_name[key].astype(value)

            megapath_nano.similar_noise_span_bp_with_tax_name = megapath_nano.similar_noise_span_bp_with_tax_name.sort_values(['target_assembly_id', 'source_assembly_id'])
            megapath_nano.similar_noise_span_bp_with_tax_name = sort_order.rename(columns={'assembly_id': 'target_assembly_id'}).merge(
                                                                                                                                  right=megapath_nano.similar_noise_span_bp_with_tax_name.set_index('target_assembly_id'),
                                                                                                                                  how='inner', 
                                                                                                                                  left_on='target_assembly_id', 
                                                                                                                                  right_index = True,
                                                                                                                                  suffixes=['', '_y'],
                                                                                                                                  validate='1:m',
                                                                                                                                 )
            megapath_nano.similar_noise_span_bp_with_tax_name.to_csv(path_or_buf=temp_noise_analysis_prefix_with_path + '.similar_region', sep='\t', header=True, index=False, columns=similar_noise_span_bp_col_name_with_tax_name)

        shutil.make_archive(file_prefix_with_path + '.noise', format=options.archive_format, root_dir=temp_noise_analysis_dir_name)
        shutil.rmtree(temp_noise_analysis_dir_name)


    megapath_nano.log.print('end')


def main():

    pandas.set_option('display.width', 120)
    pandas.set_option('display.max_colwidth', 100)

    megapath_nano = Placeholder()

    megapath_nano.global_options = {

                               # Experimental
                               'mapping_only':FLAGS.mapping_only,

                               # Non-user configurable parameters
                               'taxonomy_db':FLAGS.taxonomy_db,
                               'tool_folder':FLAGS.tool_folder,
                               'config_folder':FLAGS.config_folder,
                               'assembly_folder':FLAGS.assembly_folder,
                               'db_folder':FLAGS.db_folder,
                               'aligner':FLAGS.aligner,
                               'aligner_log':FLAGS.aligner_log,
                               'read_sim_log':FLAGS.read_sim_log,
                               'read_simulator':FLAGS.read_simulator,
                               'read_simulation_profiles':FLAGS.read_simulation_profiles,
                               'adaptor_trimming_log':FLAGS.adaptor_trimming_log,
                               'python':FLAGS.python,

                               'human_repetitive_region_filter_assembly_id':FLAGS.human_repetitive_region_filter_assembly_id,

                               'max_aligner_thread':FLAGS.max_aligner_thread,
                               'max_aligner_target_GBase_per_batch':FLAGS.max_aligner_target_GBase_per_batch,
                               'max_porechop_thread':FLAGS.max_porechop_thread,
                               'max_AMR_thread':FLAGS.max_AMR_thread,


                               'genus_height':FLAGS.genus_height,

                               # user configurable parameters
                               'human':FLAGS.human,
                               'decoy':FLAGS.decoy,
                               'species':FLAGS.species,
                               'assembly':FLAGS.assembly,

                               'min_alignment_score':FLAGS.min_alignment_score,
                               'head_crop':FLAGS.head_crop,
                               'tail_crop':FLAGS.tail_crop,
                               'min_read_length':FLAGS.min_read_length,
                               'min_read_quality':FLAGS.min_read_quality,
                               'human_filter_alignment_score_threshold':FLAGS.human_filter_alignment_score_threshold,
                               'human_filter_alignment_score_percent_threshold':FLAGS.human_filter_alignment_score_percent_threshold,
                               'decoy_filter_alignment_score_threshold':FLAGS.decoy_filter_alignment_score_threshold,
                               'decoy_filter_alignment_score_percent_threshold':FLAGS.decoy_filter_alignment_score_percent_threshold,
                               'variable_region_percent':FLAGS.variable_region_percent,
                               'expected_max_depth_stdev':FLAGS.expected_max_depth_stdev,
                               'closing_expected_max_depth_stdev':FLAGS.closing_expected_max_depth_stdev,
                               'species_id_min_aligned_bp':FLAGS.species_id_min_aligned_bp,
                               'assembly_id_min_average_depth':FLAGS.assembly_id_min_average_depth,
                               'max_alignment_noise_overlap':FLAGS.max_alignment_noise_overlap,
                               'min_alignment_length':FLAGS.min_alignment_length,
                               'unique_alignment_threshold':FLAGS.unique_alignment_threshold,
                               'number_of_genus_to_perform_noise_projection':FLAGS.number_of_genus_to_perform_noise_projection,
                               'min_percent_abundance_to_perform_noise_projection':FLAGS.min_percent_abundance_to_perform_noise_projection,
                               'noise_projection_simulated_read_length_bin_size':FLAGS.noise_projection_simulated_read_length_bin_size,
                               'noise_projection_simulated_read_length_multiplier':FLAGS.noise_projection_simulated_read_length_multiplier,
                               'noise_projection_simulated_read_error_profile':FLAGS.noise_projection_simulated_read_error_profile,
                               'noise_projection_num_read_to_simulate':FLAGS.noise_projection_num_read_to_simulate,
                               'similar_species_marker_num_genus':FLAGS.similar_species_marker_num_genus,
                               'similar_species_marker_alignment_similarity_1':FLAGS.similar_species_marker_alignment_similarity_1,
                               'similar_species_marker_aligned_region_threshold_1':FLAGS.similar_species_marker_aligned_region_threshold_1,
                               'similar_species_marker_alignment_similarity_2':FLAGS.similar_species_marker_alignment_similarity_2,
                               'similar_species_marker_aligned_region_threshold_2':FLAGS.similar_species_marker_aligned_region_threshold_2,
                               'similar_species_marker_similarity_combine_logic':FLAGS.similar_species_marker_similarity_combine_logic,
                               'microbe_repetitive_region_filter_abundance_threshold_80':FLAGS.microbe_repetitive_region_filter_abundance_threshold_80,
                               'microbe_repetitive_region_filter_abundance_threshold_90':FLAGS.microbe_repetitive_region_filter_abundance_threshold_90,
                               'microbe_repetitive_region_filter_abundance_threshold_95':FLAGS.microbe_repetitive_region_filter_abundance_threshold_95,
                               'microbe_repetitive_region_filter_abundance_threshold_98':FLAGS.microbe_repetitive_region_filter_abundance_threshold_98,
                               'microbe_repetitive_region_filter_abundance_threshold_99':FLAGS.microbe_repetitive_region_filter_abundance_threshold_99,
                               'microbe_repetitive_region_filter_abundance_threshold_99_2':FLAGS.microbe_repetitive_region_filter_abundance_threshold_99_2,
                               'microbe_repetitive_region_filter_targeted_max_span_percent':FLAGS.microbe_repetitive_region_filter_targeted_max_span_percent,
                               'microbe_repetitive_region_filter_allowed_max_span_percent':FLAGS.microbe_repetitive_region_filter_allowed_max_span_percent,
                               'microbe_repetitive_region_filter_min_average_depth':FLAGS.microbe_repetitive_region_filter_min_average_depth,
                               'microbe_repetitive_region_filter_max_span_percent_overall':FLAGS.microbe_repetitive_region_filter_max_span_percent_overall,
                               'good_alignment_threshold':FLAGS.good_alignment_threshold,

                               # system global variables
                               'debug':FLAGS.debug,
                              }

    megapath_nano.log = Log(os.sys.stderr, 'megapath_nano')

    megapath_nano.log.print('program start')
    megapath_nano.log.print_time()

    if os.path.exists(FLAGS.output_folder) == False or os.path.isdir(FLAGS.output_folder) == False:
        os.sys.exit('Output folder does not exist')
    if os.path.exists(os.path.join(megapath_nano.global_options['config_folder'], FLAGS.human)) == False or os.path.isfile(os.path.join(megapath_nano.global_options['config_folder'], FLAGS.human)) == False:
        os.sys.exit('Human genome set not found')
    if os.path.exists(os.path.join(megapath_nano.global_options['config_folder'], FLAGS.decoy)) == False or os.path.isfile(os.path.join(megapath_nano.global_options['config_folder'], FLAGS.decoy)) == False:
        os.sys.exit('Decoy genome set not found')
    if os.path.exists(os.path.join(megapath_nano.global_options['config_folder'], FLAGS.species)) == False or os.path.isfile(os.path.join(megapath_nano.global_options['config_folder'], FLAGS.species)) == False:
        os.sys.exit('Species ID genome set not found')
    if os.path.exists(os.path.join(megapath_nano.global_options['config_folder'], FLAGS.assembly)) == False or os.path.isfile(os.path.join(megapath_nano.global_options['config_folder'], FLAGS.assembly)) == False:
        os.sys.exit('Assembly ID genome set not found')


    megapath_nano.taxonomy_db_conn = sqlite3.connect(megapath_nano.global_options['taxonomy_db'])

    if FLAGS.temp_folder != '':
        megapath_nano.temp_dir_name = tempfile.mkdtemp(prefix='megapath_nano.', dir=FLAGS.temp_folder)
    else:
        megapath_nano.temp_dir_name = tempfile.mkdtemp(prefix='megapath_nano.')

    if FLAGS.RAM_folder != '':
        megapath_nano.RAM_dir_name = tempfile.mkdtemp(prefix='megapath_nano.', dir=FLAGS.RAM_folder)
    else:
        megapath_nano.RAM_dir_name = tempfile.mkdtemp(prefix='megapath_nano.')

    if megapath_nano.global_options['debug'] == False:
        atexit.register(exit_with_cleanup, temp_root=megapath_nano.temp_dir_name, RAM_root=megapath_nano.RAM_dir_name)
    else:
        pybedtools.KEEP_TEMPFILES=True

    #pybedtools.helpers.set_bedtools_path(os.path.join(megapath_nano.global_options['tool_folder'], 'bedtools'))
    pybedtools.helpers.set_tempdir(megapath_nano.RAM_dir_name)
    if FLAGS.debug == True:
        pybedtools.debug_mode(True)


    num_core = psutil.cpu_count(logical=True)
    ram_available = psutil.virtual_memory().available
    if megapath_nano.global_options['max_aligner_target_GBase_per_batch']=='auto':
        megapath_nano.global_options['max_aligner_target_GBase_per_batch'] = str(ram_available // 1024 // 1024 // 1024 // 8 )

    megapath_nano.global_options['alignerThreadOption'] = '-t ' + str(min(num_core, megapath_nano.global_options['max_aligner_thread'])) + ' -I ' + megapath_nano.global_options['max_aligner_target_GBase_per_batch'] + 'G'
    megapath_nano.global_options['porechopThreadOption'] = '-t ' + str(min(num_core, megapath_nano.global_options['max_porechop_thread']))
    megapath_nano.global_options['AMRThreadOption'] = str(min(num_core, megapath_nano.global_options['max_AMR_thread']))

    megapath_nano.log.print('Loading assembly metadata')

    megapath_nano.assembly_metadata = AssemblyMetadata(assembly_folder=megapath_nano.global_options['assembly_folder'])

    megapath_nano.log.print('Number of assemblies: {num_assembly}'.format(num_assembly=megapath_nano.assembly_metadata.get_num_assembly()))
    megapath_nano.log.print_time()

    if megapath_nano.global_options['aligner_log'] == '':
        megapath_nano.aligner_log = os.sys.stderr
    else:
        megapath_nano.aligner_log = os.open(megapath_nano.global_options['aligner_log'], flags=os.O_CREAT | os.O_WRONLY, mode=0o644)
    if megapath_nano.global_options['read_sim_log'] == '':
        megapath_nano.read_sim_log = os.sys.stderr
    else:
        if FLAGS.noise_projection == True:
            megapath_nano.read_sim_log = os.open(megapath_nano.global_options['read_sim_log'], flags=os.O_CREAT | os.O_WRONLY, mode=0o644)

    megapath_nano.query_filename_list = pandas.DataFrame(FLAGS.query.split(','), columns=['path'])

    megapath_nano.output_folder = FLAGS.output_folder
    if FLAGS.output_prefix == '':
        temp_query_file_name = os.path.split(megapath_nano.query_filename_list['path'][0])[1]
        last_dot_pos = temp_query_file_name.rfind('.')
        if temp_query_file_name[last_dot_pos + 1:] == 'gz':
            last_dot_pos = temp_query_file_name[:last_dot_pos].rfind('.')
        if last_dot_pos == 0:
            megapath_nano.output_prefix = temp_query_file_name
        else:
            megapath_nano.output_prefix = temp_query_file_name[:last_dot_pos]
    else:
        megapath_nano.output_prefix = FLAGS.output_prefix

    for file in megapath_nano.query_filename_list['path']:
        if os.path.exists(file) == True:
            megapath_nano.log.print('input query file: ' + file)
        else:
            megapath_nano.log.print('input query file: ' + file + ' not found')
            os.sys.exit()




    # step 0: 'adaptor_trimming': trim adaptors from reads

    adaptor_trimming = IO_Container()

    adaptor_trimming.I.query_filename_list = megapath_nano.query_filename_list

    if FLAGS.adaptor_trimming == True:
        step_adaptor_trimming(megapath_nano, adaptor_trimming)
        if FLAGS.output_adaptor_trimmed_query == True:
            for filename in adaptor_trimming.O.query_filename_list['path']:
                shutil.copyfile(filename, os.path.join(FLAGS.output_folder, os.path.split(filename)[1]))
    else:
        adaptor_trimming.O.query_filename_list = adaptor_trimming.I.query_filename_list

    # Output
    # adaptor_trimming.O.query_filename_list    	format: ['path'] (files reside in megapath_nano.temp_dir_name)


    # step 1: 'read_trimming_and_filter': trim and filter reads

    read_trimming_and_filter = IO_Container()

    read_trimming_and_filter.I.query_filename_list = adaptor_trimming.O.query_filename_list
    read_trimming_and_filter.I.head_crop = megapath_nano.global_options['head_crop']
    read_trimming_and_filter.I.tail_crop = megapath_nano.global_options['tail_crop']
    read_trimming_and_filter.I.min_read_length = megapath_nano.global_options['min_read_length']
    read_trimming_and_filter.I.min_read_quality = megapath_nano.global_options['min_read_quality']
    read_trimming_and_filter.I.reassign_read_id = FLAGS.reassign_read_id

    if FLAGS.read_filter == False:
        read_trimming_and_filter.I.min_read_length = 0
        read_trimming_and_filter.I.min_read_quality = 0

    if FLAGS.read_trimming == False:
        read_trimming_and_filter.I.head_crop = 0
        read_trimming_and_filter.I.tail_crop = 0

    step_read_trimming_and_filter(megapath_nano, read_trimming_and_filter)

    if FLAGS.output_trimmed_and_filtered_query == True:
        for filename in read_trimming_and_filter.O.query_filename_list['path']:
            shutil.copyfile(filename, os.path.join(FLAGS.output_folder, os.path.split(filename)[1]))

    # Output
    # read_trimming_and_filter.O.query_filename_list         format: ['path'] (files reside in megapath_nano.temp_dir_name)
    # read_trimming_and_filter.O.read_info_filename_list     format: ['path'] (files reside in megapath_nano.temp_dir_name)
    # read_trimming_and_filter.O.read_info                   format: read_info_col_name
    # read_trimming_and_filter.O.passed_read_id_list
    # read_trimming_and_filter.O.num_read
    # read_trimming_and_filter.O.num_read_passed_filter

    if read_trimming_and_filter.O.num_read_passed_filter <= 0:
        os.sys.exit('No read remained after read filter')


    # step 2: 'human_and_decoy_filter': align all reads to human and decoy; reads that can be aligned with a alignment score >= threshold are classified as human or decoy reads

    human_and_decoy_filter = IO_Container()

    human_and_decoy_filter.I.query_filename_list = read_trimming_and_filter.O.query_filename_list
    human_and_decoy_filter.I.read_id_list = read_trimming_and_filter.O.passed_read_id_list
    human_and_decoy_filter.I.human_min_alignment_score = megapath_nano.global_options['human_filter_alignment_score_threshold']
    human_and_decoy_filter.I.human_min_alignment_score_percent = megapath_nano.global_options['human_filter_alignment_score_percent_threshold']
    human_and_decoy_filter.I.decoy_min_alignment_score = megapath_nano.global_options['decoy_filter_alignment_score_threshold']
    human_and_decoy_filter.I.decoy_min_alignment_score_percent = megapath_nano.global_options['decoy_filter_alignment_score_percent_threshold']

    if FLAGS.human_filter == True and megapath_nano.global_options['human'] != '':
        megapath_nano.log.print('Human filter genome set: ' + megapath_nano.global_options['human'])
        human_and_decoy_filter.I.human_assembly_list = read_genome_set(global_options=megapath_nano.global_options, genome_set_name=megapath_nano.global_options['human'])
    else:
        human_and_decoy_filter.I.human_assembly_list = pandas.DataFrame(columns=['assembly_id'])

    if FLAGS.decoy_filter == True and megapath_nano.global_options['decoy'] != '':
        megapath_nano.log.print('Decoy filter genome set: ' + megapath_nano.global_options['decoy'])
        human_and_decoy_filter.I.decoy_assembly_list = read_genome_set(global_options=megapath_nano.global_options, genome_set_name=megapath_nano.global_options['decoy'])
    else:
        human_and_decoy_filter.I.decoy_assembly_list = pandas.DataFrame(columns=['assembly_id'])


    if human_and_decoy_filter.I.human_assembly_list.empty == False or human_and_decoy_filter.I.decoy_assembly_list.empty == False:
        step_human_and_decoy_filter(megapath_nano, human_and_decoy_filter)
        if FLAGS.output_human_and_decoy_filtered_query == True:
            for filename in human_and_decoy_filter.O.query_filename_list['path']:
                shutil.copyfile(filename, os.path.join(FLAGS.output_folder, os.path.split(filename)[1])) #filename ended with .human_and_decoy_filtered
    else:
        human_and_decoy_filter.O.query_filename_list = human_and_decoy_filter.I.query_filename_list
        human_and_decoy_filter.O.human_best_align_list = pandas.DataFrame(columns=align_list_col_name)
        for key, value in align_list_col_type.items():
            human_and_decoy_filter.O.human_best_align_list[key] = human_and_decoy_filter.O.human_best_align_list[key].astype(value)

        human_and_decoy_filter.O.decoy_best_align_list = pandas.DataFrame(columns=align_list_col_name)
        for key, value in align_list_col_type.items():
            human_and_decoy_filter.O.decoy_best_align_list[key] = human_and_decoy_filter.O.decoy_best_align_list[key].astype(value)
        #human_and_decoy_filter.O.microbe_best_align_list to be filled in the next step
        human_and_decoy_filter.O.human_read_id_list = pandas.DataFrame(columns=['read_id', 'read_length'])
        human_and_decoy_filter.O.human_read_id_list['read_length'] = human_and_decoy_filter.O.human_read_id_list['read_length'].astype('int')
        human_and_decoy_filter.O.decoy_read_id_list = pandas.DataFrame(columns=['read_id', 'read_length'])
        human_and_decoy_filter.O.decoy_read_id_list['read_length'] = human_and_decoy_filter.O.decoy_read_id_list['read_length'].astype('int')
        human_and_decoy_filter.O.microbe_read_id_list = read_trimming_and_filter.O.passed_read_id_list
        #read_id

    # Output
    # human_and_decoy_filter.O.query_filename_list      format: ['path'] (files reside in megapath_nano.temp_dir_name)
    # human_and_decoy_filter.O.human_best_align_list    format: align_list_col_name
    # human_and_decoy_filter.O.decoy_best_align_list    format: align_list_col_name
    # human_and_decoy_filter.O.microbe_best_align_list    format: align_list_col_name
    # human_and_decoy_filter.O.human_read_id_list       format: ['read_id', 'read_length']
    # human_and_decoy_filter.O.decoy_read_id_list       format: ['read_id', 'read_length']
    # human_and_decoy_filter.O.microbe_read_id_list       format: ['read_id', 'read_length']
    
    
    if FLAGS.all_taxon_steps == False:
        print('Finished all taxon steps')
        os.sys.exit()
    
    
    
    # step 3: 'placement_to_species': align filtered reads to all genome sets; species with number of reads best aligned to it >= threshold will be included

    placement_to_species = IO_Container()

    placement_to_species.I.query_filename_list = human_and_decoy_filter.O.query_filename_list
    placement_to_species.I.species_ID_min_aligned_bp = megapath_nano.global_options['species_id_min_aligned_bp']
    placement_to_species.I.target_assembly_list = read_genome_set(global_options=megapath_nano.global_options, genome_set_name=megapath_nano.global_options['species'])

    if FLAGS.output_PAF == True:
        placement_to_species.I.paf_path_and_prefix = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix + '.species')
    else:
        placement_to_species.I.paf_path_and_prefix = None

    megapath_nano.log.print('Species identification genome set: ' + megapath_nano.global_options['species'])

    step_placement_to_species(megapath_nano, placement_to_species)

    # Output
    # placement_to_species.O.best_align_list            format: align_list_col_name
    # placement_to_species.O.align_list                 format: align_list_col_name
    # placement_to_species.O.aligned_species_list       format: ['species_tax_id', 'aligned_bp']
    # placement_to_species.O.selected_species_list      format: ['species_tax_id', 'aligned_bp']
    # placement_to_species.O.read_id_species_id         format: ['read_id', 'species_tax_id']

    # When no human or decoy has been specified
    if human_and_decoy_filter.I.human_assembly_list.empty == True and human_and_decoy_filter.I.decoy_assembly_list.empty == True:
        human_and_decoy_filter.O.microbe_best_align_list = placement_to_species.O.best_align_list


    if FLAGS.assembly_selection == True:
        # step X: 'placement_to_assembly': align reads that best align to a species to all assembly of the species and then choose the best assembly



        placement_to_assembly.I.query_filename_list = human_and_decoy_filter.O.query_filename_list
        placement_to_assembly.I.species_list = placement_to_species.O.selected_species_list
        placement_to_assembly.I.read_id_species_id = placement_to_species.O.read_id_species_id
        placement_to_assembly.I.target_assembly_list = read_genome_set(global_options=megapath_nano.global_options, genome_set_name=megapath_nano.global_options['assembly'])
        placement_to_assembly.I.species_id_assembly_id = placement_to_species.I.target_assembly_list[['assembly_id']]

        megapath_nano.log.print('Assembly identification genome set: ' + megapath_nano.global_options['assembly'])

        step_placement_to_assembly(megapath_nano, placement_to_assembly)

        # Output
        # placement_to_assembly.O.align_list           		format: align_list_col_name
        # placement_to_assembly.O.num_assembly_candidate


        # step X: 'assembly_selection': align reads that best align to a species to all assembly of the species and then choose the best assembly

        assembly_selection = IO_Container()

        assembly_selection.I.species_list = placement_to_species.O.selected_species_list
        assembly_selection.I.assembly_align_list = placement_to_assembly.O.align_list
        assembly_selection.I.species_align_list = placement_to_species.O.align_list
        assembly_selection.I.assembly_ID_min_average_depth = megapath_nano.global_options['assembly_id_min_average_depth']
        assembly_selection.I.good_align_threshold = megapath_nano.global_options['good_alignment_threshold']
        assembly_selection.I.read_id_species_id = placement_to_species.O.read_id_species_id

        step_assembly_selection(megapath_nano, assembly_selection)

        # Output
        # assembly_selection.O.align_stat       	format: align_stat_col_name
        # assembly_selection.O.align_list       	format: align_stat_col_name
        # assembly_selection.O.best_align_list  	format: align_stat_col_name
        # assembly_selection.O.good_align_list  	format: align_stat_col_name
        # assembly_selection.O.assembly_list    	format: align_stat_col_name


        # step X: 'align_assembly_set': align filtered reads to assembly_set

        align_assembly_set = IO_Container()

        align_assembly_set.I.query_filename_list = human_and_decoy_filter.O.query_filename_list
        align_assembly_set.I.assembly_list = assembly_selection.O.assembly_list[['assembly_id']]
        align_assembly_set.I.species_id_assembly_id = placement_to_species.I.target_assembly_list[['assembly_id']]
        align_assembly_set.I.species_align_list = placement_to_species.O.align_list
        
        if FLAGS.output_PAF == True:
            align_assembly_set.I.paf_path_and_prefix = os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix + '.assembly')
        else:
            align_assembly_set.I.paf_path_and_prefix = None

        step_align_assembly_set(megapath_nano, align_assembly_set)

        # Output
        # align_assembly_set.O.align_list       	format: align_list_col_name
        # align_assembly_set.O.best_align_list  	format: align_list_col_name
    else:
        placement_to_assembly = IO_Container()
        placement_to_assembly.I.target_assembly_list = read_genome_set(global_options=megapath_nano.global_options, genome_set_name=megapath_nano.global_options['assembly'])
        
        assembly_selection = IO_Container()
        species_align_stat = align_list_to_align_stat_by_assembly_id(
                                                                 assembly_metadata=megapath_nano.assembly_metadata,
                                                                 log=megapath_nano.log,
                                                                 align_list=placement_to_species.O.align_list
                                                                )
        species_align_stat = species_align_stat.sort_values(['species_tax_id', 'adjusted_average_depth', 'alignment_score_tiebreaker']).drop_duplicates(subset=['species_tax_id'], keep='last')
        species_align_stat = species_align_stat.merge(
                                                      right=placement_to_species.O.selected_species_list[['species_tax_id']].set_index('species_tax_id'),
                                                      how='inner',
                                                      left_on='species_tax_id',
                                                      right_index = True,
                                                      suffixes=['', '_y'],
                                                      validate='1:1',
                                                     )

        assembly_selection.O.best_align_list = placement_to_species.O.best_align_list
        assembly_selection.O.assembly_list = species_align_stat.sort_values(['species_tax_id', 'adjusted_average_depth', 'alignment_score_tiebreaker']).drop_duplicates(subset=['species_tax_id'], keep='last').copy()

        align_assembly_set = IO_Container()
        align_assembly_set.O.align_list=placement_to_species.O.align_list
        align_assembly_set.O.best_align_list=placement_to_species.O.best_align_list


    # step 7: 'raw_stat': assembly raw best align list

    raw_stat = IO_Container()

    raw_stat.I.best_align_list = align_assembly_set.O.best_align_list
    raw_stat.I.huamn_and_decoy_best_align_list = human_and_decoy_filter.O.microbe_best_align_list   # microbe_best_align_list contains alignment to human and decoy

    step_raw_stat(megapath_nano, raw_stat)

    # Output
    # raw_stat.O.best_align_list          format: align_list_col_name


    # step X: 'variable_region': find variable regions of each aligned assembly
    
    variable_region = IO_Container()

    variable_region.I.selected_assembly_list = assembly_selection.O.assembly_list
    variable_region.I.target_assembly_list = placement_to_assembly.I.target_assembly_list
    variable_region.I.variable_region_percent = megapath_nano.global_options['variable_region_percent']

    if FLAGS.variable_region_adjustment == True:
        step_variable_region(megapath_nano, variable_region)
        variable_region.O.noise_stat.to_csv(path_or_buf=os.path.join(megapath_nano.output_folder, megapath_nano.output_prefix + '.variable_region'), sep='\t', header=True, index=False, columns=['assembly_id', 'variable_span_bp', 'variable_span_percent'])
    else:
        variable_region.O.variable_region_bed = BedTool('', from_string=True)
        variable_region.O.noise_stat = assembly_selection.O.assembly_list[['assembly_id']].assign(variable_span_bp = lambda x: 0, variable_span_percent = lambda x: 0).copy()
        variable_region.O.noise_summary = variable_region.O.noise_stat['variable_span_percent'].describe()

    # Output
    # variable_region.O.noise_bed
    # variable_region.O.noise_stat         format: ['assembly_id', 'variable_span_bp', 'variable_span_percent']
    # variable_region.O.noise_summary


    # step 8: 'spike_filter': calculate expected max_depth using best align_list and mark regions with depth > max_depth as noise region
    
    spike_filter = IO_Container()

    spike_filter.I.align_list = align_assembly_set.O.best_align_list                     # use best alignment to detect spike
    spike_filter.I.expected_max_depth_stdev = megapath_nano.global_options['expected_max_depth_stdev']

    if FLAGS.spike_filter == True:
        step_spike_filter(megapath_nano, spike_filter)
    else:
        spike_filter.O.noise_bed = BedTool('', from_string=True)
        spike_filter.O.noise_stat = assembly_selection.O.assembly_list[['assembly_id']].assign(spike_span_bp = lambda x: 0, spike_span_percent = lambda x: 0).copy()
        spike_filter.O.noise_summary = spike_filter.O.noise_stat['spike_span_percent'].describe()

    # Output
    # spike_filter.O.noise_bed
    # spike_filter.O.noise_stat         format: ['assembly_id', 'spike_span_bp', 'spike_span_percent']
    # spike_filter.O.noise_summary


    # step 9: 'human_repetitive_region_filter': mark regions similar to human as noise region

    human_repetitive_region_filter = IO_Container()

    human_repetitive_region_filter.I.assembly_list = assembly_selection.O.assembly_list
    
    if FLAGS.human_repetitive_region_filter == True:
        step_human_repetitive_region_filter(megapath_nano, human_repetitive_region_filter)
    else:
        human_repetitive_region_filter.O.noise_bed = BedTool('', from_string=True)
        human_repetitive_region_filter.O.noise_stat = human_repetitive_region_filter.I.assembly_list[['assembly_id']].assign(human_span_bp = lambda x: 0, human_span_percent = lambda x: 0).copy()

    # Output
    # human_repetitive_region_filter.O.noise_bed
    # human_repetitive_region_filter.O.noise_stat          format: ['assembly_id', 'human_span_bp', 'human_span_percent']
    # human_repetitive_region_filter.O.noise_summary


    # step 10: 'approx_abundance': remove noise and calculate abundance

    approx_abundance = IO_Container()

    approx_abundance.I.align_list = align_assembly_set.O.align_list
    approx_abundance.I.noise_bed = merge_bed_with_assembly_id(bed_list=[spike_filter.O.noise_bed, human_repetitive_region_filter.O.noise_bed])
    approx_abundance.I.max_align_noise_overlap = megapath_nano.global_options['max_alignment_noise_overlap']

    step_approx_abundance(megapath_nano, approx_abundance)

    # Output
    # approx_abundance.O.align_stat         format: align_stat_col_name


    # step 11: 'microbe_repetitive_region_filter': align high abundant species to low abundant species and mark aligned regions as noise regions

    microbe_repetitive_region_filter = IO_Container()

    microbe_repetitive_region_filter.I.align_stat = approx_abundance.O.align_stat
    microbe_repetitive_region_filter.I.assembly_list = assembly_selection.O.assembly_list
    microbe_repetitive_region_filter.I.abundance_threshold_80 = megapath_nano.global_options['microbe_repetitive_region_filter_abundance_threshold_80']
    microbe_repetitive_region_filter.I.abundance_threshold_90 = megapath_nano.global_options['microbe_repetitive_region_filter_abundance_threshold_90']
    microbe_repetitive_region_filter.I.abundance_threshold_95 = megapath_nano.global_options['microbe_repetitive_region_filter_abundance_threshold_95']
    microbe_repetitive_region_filter.I.abundance_threshold_98 = megapath_nano.global_options['microbe_repetitive_region_filter_abundance_threshold_98']
    microbe_repetitive_region_filter.I.abundance_threshold_99 = megapath_nano.global_options['microbe_repetitive_region_filter_abundance_threshold_99']
    microbe_repetitive_region_filter.I.abundance_threshold_99_2 = megapath_nano.global_options['microbe_repetitive_region_filter_abundance_threshold_99_2']
    microbe_repetitive_region_filter.I.targeted_max_span_percent = megapath_nano.global_options['microbe_repetitive_region_filter_targeted_max_span_percent']
    microbe_repetitive_region_filter.I.allowed_max_span_percent = megapath_nano.global_options['microbe_repetitive_region_filter_allowed_max_span_percent']
    microbe_repetitive_region_filter.I.min_average_depth = megapath_nano.global_options['microbe_repetitive_region_filter_min_average_depth']
    microbe_repetitive_region_filter.I.max_span_percent_overall = megapath_nano.global_options['microbe_repetitive_region_filter_max_span_percent_overall']
    microbe_repetitive_region_filter.I.genus_height = megapath_nano.global_options['genus_height']

    if FLAGS.microbe_repetitive_region_filter == True:
        step_microbe_repetitive_region_filter(megapath_nano, microbe_repetitive_region_filter)
    else:
        microbe_repetitive_region_filter.O.noise_bed = BedTool('', from_string=True)
        microbe_repetitive_region_filter.O.noise_stat = assembly_selection.O.assembly_list[['assembly_id']].assign(microbe_span_bp = lambda x: 0, microbe_span_percent = lambda x: 0).copy()
        microbe_repetitive_region_filter.O.source_target_noise_span_bp = pandas.DataFrame(columns=similar_noise_span_bp_col_name)
        for key, value in similar_noise_span_bp_col_type.items():
            microbe_repetitive_region_filter.O.source_target_noise_span_bp[key] = microbe_repetitive_region_filter.O.source_target_noise_span_bp[key].astype(value)

    # Output
    # microbe_repetitive_region_filter.O.noise_bed
    # microbe_repetitive_region_filter.O.noise_stat                    format: ['assembly_id', 'microbe_span_bp', 'microbe_span_percent']
    # microbe_repetitive_region_filter.O.noise_summary
    # microbe_repetitive_region_filter.O.source_target_noise_span_bp   format: similar_noise_span_bp_col_name


    # step 12: 'noise_removal': remove alignments that have overlap with noise region over max_noise_overlap

    noise_removal = IO_Container()

    noise_removal.I.align_list = align_assembly_set.O.align_list
    noise_removal.I.noise_bed = merge_bed_with_assembly_id(bed_list=[spike_filter.O.noise_bed, human_repetitive_region_filter.O.noise_bed, microbe_repetitive_region_filter.O.noise_bed])
    noise_removal.I.max_align_noise_overlap = megapath_nano.global_options['max_alignment_noise_overlap']
    noise_removal.I.non_zero_approx_abundence_assembly_id = approx_abundance.O.align_stat.query('adjusted_average_depth > 0')[['assembly_id']].merge(
                                                                                                                                                     right=raw_stat.O.best_align_list[['assembly_id']].drop_duplicates().set_index(['assembly_id']),
                                                                                                                                                     how='inner', 
                                                                                                                                                     left_on='assembly_id', 
                                                                                                                                                     right_index = True,
                                                                                                                                                     suffixes=['', '_y'],
                                                                                                                                                     validate='1:1',
                                                                                                                                                    )

    step_noise_removal(megapath_nano, noise_removal)

    # Output
    # noise_removal.O.align_list                    format: align_list_col_name
    # noise_removal.O.num_align_before
    # noise_removal.O.num_align_after
    # noise_removal.O.noise_bed


    # step 13: 'short_alignment_removal'

    short_alignment_removal = IO_Container()

    short_alignment_removal.I.align_list = noise_removal.O.align_list
    short_alignment_removal.I.min_align_length = megapath_nano.global_options['min_alignment_length']

    if FLAGS.short_alignment_filter == True:
        step_short_alignment_removal(megapath_nano, short_alignment_removal)
    else:
        short_alignment_removal.O.align_list = short_alignment_removal.I.align_list

    # Output
    # short_alignment_removal.O.align_list          format: align_list_col_name
    # short_alignment_removal.O.num_read_before
    # short_alignment_removal.O.num_read_after


    # step 14: 'closing_spike_filter': a more round of spike filter where noise reads (instead of noise alignments) are removed
    
    closing_spike_filter = IO_Container()

    closing_spike_filter.I.align_list = short_alignment_removal.O.align_list
    closing_spike_filter.I.align_list_with_short_alignment = noise_removal.O.align_list
    closing_spike_filter.I.noise_bed = noise_removal.O.noise_bed
    closing_spike_filter.I.expected_max_depth_stdev = megapath_nano.global_options['closing_expected_max_depth_stdev']
    closing_spike_filter.I.max_align_noise_overlap = megapath_nano.global_options['max_alignment_noise_overlap']

    if FLAGS.closing_spike_filter == True:
        step_closing_spike_filter(megapath_nano, closing_spike_filter)
    else:
        closing_spike_filter.O.closing_spike_noise_bed = BedTool('', from_string=True)
        closing_spike_filter.O.noise_bed = closing_spike_filter.I.noise_bed
        closing_spike_filter.O.noise_stat = assembly_selection.O.assembly_list[['assembly_id']].assign(closing_spike_span_bp = lambda x: 0, closing_spike_span_percent = lambda x: 0).copy()
        closing_spike_filter.O.noise_summary = closing_spike_filter.O.noise_stat['closing_spike_span_percent'].describe()
        closing_spike_filter.O.align_list = closing_spike_filter.I.align_list

    # Output
    # closing_spike_filter.O.closing_spike_noise_bed
    # closing_spike_filter.O.noise_bed
    # closing_spike_filter.O.noise_stat               format: ['assembly_id', 'closing_spike_span_bp', 'closing_spike_span_percent']
    # closing_spike_filter.O.noise_summary
    # closing_spike_filter.O.align_list


    # step 15: 'combine_with_human_and_decoy'

    combine_with_human_and_decoy = IO_Container() 

    combine_with_human_and_decoy.I.align_list = closing_spike_filter.O.align_list
    combine_with_human_and_decoy.I.huamn_and_decoy_best_align_list = human_and_decoy_filter.O.microbe_best_align_list   # microbe_best_align_list contains alignment to human and decoy

    step_combine_with_human_and_decoy(megapath_nano, combine_with_human_and_decoy)

    # Output
    # combine_with_human_and_decoy.O.align_list             format: align_list_col_name


    # step 16: 'best_alignment'

    best_alignment = IO_Container()

    best_alignment.I.align_list = combine_with_human_and_decoy.O.align_list

    step_best_alignment(megapath_nano, best_alignment)

    # Output
    # best_alignment.O.best_align_list              format: align_list_col_name


    # step 17: 'separate_human_and_decoy'

    separate_human_and_decoy = IO_Container()

    separate_human_and_decoy.I.best_align_list = best_alignment.O.best_align_list
    separate_human_and_decoy.I.human_best_align_list = human_and_decoy_filter.O.human_best_align_list
    separate_human_and_decoy.I.decoy_best_align_list = human_and_decoy_filter.O.decoy_best_align_list
    separate_human_and_decoy.I.human_read_id_list = human_and_decoy_filter.O.human_read_id_list
    separate_human_and_decoy.I.decoy_read_id_list = human_and_decoy_filter.O.decoy_read_id_list
    separate_human_and_decoy.I.read_id_list = read_trimming_and_filter.O.passed_read_id_list
    separate_human_and_decoy.I.human_assembly_list = human_and_decoy_filter.I.human_assembly_list
    separate_human_and_decoy.I.decoy_assembly_list = human_and_decoy_filter.I.decoy_assembly_list

    step_separate_human_and_decoy(megapath_nano, separate_human_and_decoy)

    # Output
    # separate_human_and_decoy.O.best_align_list        format: align_list_col_name
    # separate_human_and_decoy.O.human_best_align_list  format: align_list_col_name
    # separate_human_and_decoy.O.decoy_best_align_list  format: align_list_col_name
    # separate_human_and_decoy.O.read_list              format: aligned_read_list_col_name
    # separate_human_and_decoy.O.assembly_list          format: ['assembly_id']


    # step 18: 'unique_alignment'

    unique_alignment = IO_Container()

    unique_alignment.I.best_align_list = separate_human_and_decoy.O.best_align_list
    unique_alignment.I.human_best_align_list = separate_human_and_decoy.O.human_best_align_list
    unique_alignment.I.decoy_best_align_list = separate_human_and_decoy.O.decoy_best_align_list
    unique_alignment.I.align_list = short_alignment_removal.O.align_list
    unique_alignment.I.unique_align_threshold = megapath_nano.global_options['unique_alignment_threshold']

    if FLAGS.unique_alignment == True:
        step_unique_alignment(megapath_nano, unique_alignment)
    else:
        unique_alignment.O.best_align_list = unique_alignment.I.best_align_list

    # Output
    # unique_alignment.O.best_align_list                format: align_list_col_name
    # unique_alignment.O.num_align_before
    # unique_alignment.O.num_align_after


    # step X: 'noise_projection'

    noise_projection = IO_Container()

    noise_projection.I.best_align_list = separate_human_and_decoy.O.best_align_list
    noise_projection.I.noise_bed = noise_removal.O.noise_bed
    noise_projection.I.num_genus = megapath_nano.global_options['number_of_genus_to_perform_noise_projection']
    noise_projection.I.min_percent_abundance = megapath_nano.global_options['min_percent_abundance_to_perform_noise_projection']
    noise_projection.I.read_length_bin_size = megapath_nano.global_options['noise_projection_simulated_read_length_bin_size']
    noise_projection.I.read_length_multiplier = megapath_nano.global_options['noise_projection_simulated_read_length_multiplier']
    noise_projection.I.error_profile = megapath_nano.global_options['noise_projection_simulated_read_error_profile']
    noise_projection.I.num_read_to_simulate = megapath_nano.global_options['noise_projection_num_read_to_simulate']
    noise_projection.I.min_read_length = megapath_nano.global_options['min_read_length']
    noise_projection.I.genus_height = megapath_nano.global_options['genus_height']

    if FLAGS.noise_projection == True:
        step_noise_projection(megapath_nano, noise_projection)
    else:
        noise_projection.O.noise_aligned_bp = pandas.DataFrame(columns=noise_projection_aligned_bp_name)
        for key, value in noise_projection_aligned_bp_type.items():
            noise_projection.O.noise_aligned_bp[key] = noise_projection.O.noise_aligned_bp[key].astype(value)

    # Output
    # noise_projection.O.noise_aligned_bp                format: assembly_id, as_noise_source, as_noise_target, projected_noise_aligned_bp


    # step 19: 'similar_species_marker'

    similar_species_marker = IO_Container()

    similar_species_marker.I.best_align_list = separate_human_and_decoy.O.best_align_list
    similar_species_marker.I.noise_bed = noise_removal.O.noise_bed
    similar_species_marker.I.num_genus = megapath_nano.global_options['similar_species_marker_num_genus']
    similar_species_marker.I.alignment_similarity_1 = megapath_nano.global_options['similar_species_marker_alignment_similarity_1']
    similar_species_marker.I.aligned_region_threshold_1 = megapath_nano.global_options['similar_species_marker_aligned_region_threshold_1']
    similar_species_marker.I.alignment_similarity_2 = megapath_nano.global_options['similar_species_marker_alignment_similarity_2']
    similar_species_marker.I.aligned_region_threshold_2 = megapath_nano.global_options['similar_species_marker_aligned_region_threshold_2']
    similar_species_marker.I.similarity_combine_logic = megapath_nano.global_options['similar_species_marker_similarity_combine_logic']
    similar_species_marker.I.genus_height = megapath_nano.global_options['genus_height']

    if FLAGS.similar_species_marker == True:
        step_similar_species_marker(megapath_nano, similar_species_marker)
    else:
        similar_species_marker.O.similar_species_marker = pandas.DataFrame(columns=similar_species_marker_name_with_assembly_id)
        for key, value in similar_species_marker_type_with_assembly_id.items():
            similar_species_marker.O.similar_species_marker[key] = similar_species_marker.O.similar_species_marker[key].astype(value)

    # Output
    # similar_species_marker.O.similar_species_marker                format: similar_species_marker_name_with_assembly_id


    # step 20: 'noise_detection_statistics'

    noise_detection_statistics = IO_Container()

    noise_detection_statistics.I.assembly_list = assembly_selection.O.assembly_list
    noise_detection_statistics.I.spike_noise_stat = spike_filter.O.noise_stat
    noise_detection_statistics.I.human_noise_stat = human_repetitive_region_filter.O.noise_stat
    noise_detection_statistics.I.microbe_noise_stat = microbe_repetitive_region_filter.O.noise_stat
    noise_detection_statistics.I.noise_bed = noise_removal.O.noise_bed
    noise_detection_statistics.I.closing_spike_noise_stat = closing_spike_filter.O.noise_stat

    if FLAGS.output_noise_stat == True:
        step_noise_detection_statistics(megapath_nano, noise_detection_statistics)
    else:
        noise_detection_statistics.O.noise_stat = pandas.DataFrame(columns=noise_detection_stat_col_name_by_assembly_id)
        for key, value in noise_detection_stat_col_type_by_assembly_id.items():
            noise_detection_statistics.O.noise_stat[key] = noise_detection_statistics.O.noise_stat[key].astype(value)        

    # Output
    # noise_detection_statistics.O.noise_stat          format: noise_detection_stat_col_name_by_assembly_id


    # step 21: 'noise_removal_statistics'

    noise_removal_statistics = IO_Container()

    noise_removal_statistics.I.assembly_list = assembly_selection.O.assembly_list
    noise_removal_statistics.I.best_align_list = raw_stat.O.best_align_list          # take the best_align_list before noise removal
    noise_removal_statistics.I.spike_noise_bed = spike_filter.O.noise_bed
    noise_removal_statistics.I.human_noise_bed = human_repetitive_region_filter.O.noise_bed
    noise_removal_statistics.I.microbe_noise_bed = microbe_repetitive_region_filter.O.noise_bed
    noise_removal_statistics.I.closing_spike_noise_bed = closing_spike_filter.O.closing_spike_noise_bed
    noise_removal_statistics.I.noise_bed = closing_spike_filter.O.noise_bed
    noise_removal_statistics.I.max_align_noise_overlap = megapath_nano.global_options['max_alignment_noise_overlap']
    if FLAGS.short_alignment_filter == True:
        noise_removal_statistics.I.min_align_length = megapath_nano.global_options['min_alignment_length']
    else:
        noise_removal_statistics.I.min_align_length = 0

    if FLAGS.output_noise_stat == True:
        step_noise_removal_statistics(megapath_nano, noise_removal_statistics)
    else:
        noise_removal_statistics.O.noise_best_align_list = pandas.DataFrame(columns=align_list_col_name)
        for key, value in align_list_col_type.items():
            noise_removal_statistics.O.noise_best_align_list[key] = noise_removal_statistics.O.noise_best_align_list[key].astype(value)        


    # Output
    # noise_removal_statistics.O.noise_best_align_list      format: align_list_col_name
    # noise_removal_statistics.O.noise_stat                 format: noise_removal_stat_col_name_by_assembly_id


    # step 22: 'noise_source_statistics'

    noise_source_statistics = IO_Container()

    noise_source_statistics.I.noise_best_align_list = noise_removal_statistics.O.noise_best_align_list
    noise_source_statistics.I.read_list = separate_human_and_decoy.O.read_list
    noise_source_statistics.I.best_align_list = separate_human_and_decoy.O.best_align_list

    if FLAGS.output_noise_stat == True:
        step_noise_source_statistics(megapath_nano, noise_source_statistics)
    else:
        noise_source_statistics.O.noise_source_stat = pandas.DataFrame(columns=noise_source_stat_col_name)
        for key, value in noise_source_stat_col_type.items():
            noise_source_statistics.O.noise_source_stat[key] = noise_source_statistics.O.noise_source_stat[key].astype(value)        
    

    # Output

    # noise_source_statistics.O.noise_source_stat               format: noise_source_stat_col_name


    # step 23: 'max_adjusted_abundance'

    max_adjusted_abundance = IO_Container()

    max_adjusted_abundance.I.noise_best_align_list = noise_removal_statistics.O.noise_best_align_list
    max_adjusted_abundance.I.best_align_list = separate_human_and_decoy.O.best_align_list

    step_max_adjusted_abundance(megapath_nano, max_adjusted_abundance)
  

    # Output

    # max_adjusted_abundance.O.abundance_stat               format: ['assembly_id', 'max_adjusted_total_aligned_bp', 'max_adjusted_average_depth', 'max_adjusted_covered_percent']


    # step 24: 'read_statistics'

    read_statistics = IO_Container()

    read_statistics.I.read_list = separate_human_and_decoy.O.read_list
    read_statistics.I.read_info = read_trimming_and_filter.O.read_info
    read_statistics.I.quality_score_bin_size = FLAGS.quality_score_bin_size
    read_statistics.I.read_length_bin_size = FLAGS.read_length_bin_size

    step_read_statistics(megapath_nano, read_statistics)

    # Output
    # read_statistics.O.read_list
    # read_statistics.O.read_stat
    # read_statistics.O.quality_score_bin
    # read_statistics.O.read_length_bin
    # read_statistics.O.passed_quality_score_bin
    # read_statistics.O.passed_read_length_bin
    # read_statistics.O.human_quality_score_bin
    # read_statistics.O.human_read_length_bin
    # read_statistics.O.decoy_quality_score_bin
    # read_statistics.O.decoy_read_length_bin
    # read_statistics.O.microbe_quality_score_bin
    # read_statistics.O.microbe_read_length_bin
    # read_statistics.O.aligned_quality_score_bin
    # read_statistics.O.aligned_read_length_bin
    # read_statistics.O.unaligned_quality_score_bin
    # read_statistics.O.unaligned_read_length_bin


    # step 25: 'format_output'

    # List of pathogens

    megapath_nano.assembly_list = assembly_selection.O.assembly_list[['assembly_id']]

    # Genome sets used

    megapath_nano.human_genome_set = megapath_nano.global_options['human']
    megapath_nano.decoy_genome_set = megapath_nano.global_options['decoy']
    megapath_nano.species_id_genome_set = megapath_nano.global_options['species']
    megapath_nano.assembly_id_genome_set = megapath_nano.global_options['assembly']

    megapath_nano.human_assembly_list = human_and_decoy_filter.I.human_assembly_list[['assembly_id']]
    megapath_nano.decoy_assembly_list = human_and_decoy_filter.I.decoy_assembly_list[['assembly_id']]
    megapath_nano.species_id_assembly_list = placement_to_species.I.target_assembly_list[['assembly_id']]
    megapath_nano.assembly_id_assembly_list = placement_to_assembly.I.target_assembly_list[['assembly_id']]

    # Per-read data

    megapath_nano.read_list = read_statistics.O.read_list
    megapath_nano.human_best_align_list = separate_human_and_decoy.O.human_best_align_list
    megapath_nano.decoy_best_align_list = separate_human_and_decoy.O.decoy_best_align_list
    megapath_nano.align_list = combine_with_human_and_decoy.O.align_list
    megapath_nano.best_align_list = separate_human_and_decoy.O.best_align_list
    megapath_nano.unique_align_list = unique_alignment.O.best_align_list
    megapath_nano.max_adjusted_abundance = max_adjusted_abundance.O.abundance_stat
    megapath_nano.id_best_align_list = assembly_selection.O.best_align_list
    megapath_nano.raw_best_align_list = raw_stat.O.best_align_list
    megapath_nano.noise_removal_best_align_list = noise_removal_statistics.O.noise_best_align_list

    # Similar species marker

    megapath_nano.similar_species_marker = similar_species_marker.O.similar_species_marker

    # Read statistics

    megapath_nano.read_stat = read_statistics.O.read_stat
    megapath_nano.quality_score_bin = read_statistics.O.quality_score_bin
    megapath_nano.read_length_bin = read_statistics.O.read_length_bin
    megapath_nano.passed_quality_score_bin = read_statistics.O.passed_quality_score_bin
    megapath_nano.passed_read_length_bin = read_statistics.O.passed_read_length_bin
    megapath_nano.human_quality_score_bin = read_statistics.O.human_quality_score_bin
    megapath_nano.human_read_length_bin = read_statistics.O.human_read_length_bin
    megapath_nano.decoy_quality_score_bin = read_statistics.O.decoy_quality_score_bin
    megapath_nano.decoy_read_length_bin = read_statistics.O.decoy_read_length_bin
    megapath_nano.microbe_quality_score_bin = read_statistics.O.microbe_quality_score_bin
    megapath_nano.microbe_read_length_bin = read_statistics.O.microbe_read_length_bin
    megapath_nano.aligned_quality_score_bin = read_statistics.O.aligned_quality_score_bin
    megapath_nano.aligned_read_length_bin = read_statistics.O.aligned_read_length_bin
    megapath_nano.unaligned_quality_score_bin = read_statistics.O.unaligned_quality_score_bin
    megapath_nano.unaligned_read_length_bin = read_statistics.O.unaligned_read_length_bin


    # Noise analysis

    megapath_nano.noise_bed = closing_spike_filter.O.noise_bed
    megapath_nano.spike_noise_bed = spike_filter.O.noise_bed
    megapath_nano.human_repetitive_region_noise_bed = human_repetitive_region_filter.O.noise_bed
    megapath_nano.microbe_repetitive_region_noise_bed = microbe_repetitive_region_filter.O.noise_bed
    megapath_nano.closing_spike_noise_bed = closing_spike_filter.O.closing_spike_noise_bed

    megapath_nano.noise_detection_stat = noise_detection_statistics.O.noise_stat                                 # format: noise_detection_stat_col_name_by_assembly_id
    megapath_nano.noise_removal_stat = noise_removal_statistics.O.noise_align_stat                               # format: noise_removal_stat_col_name_by_assembly_id
    megapath_nano.noise_source_stat = noise_source_statistics.O.noise_source_stat                                # format: noise_source_stat_col_name
    megapath_nano.source_target_noise_span_bp = microbe_repetitive_region_filter.O.source_target_noise_span_bp             # format: similar_noise_span_bp_col_name


    step_format_output(megapath_nano, FLAGS)


    megapath_nano.taxonomy_db_conn.close()

    if FLAGS.debug == False:
        try:
            shutil.rmtree(megapath_nano.temp_dir_name)
        except(FileNotFoundError):
            pass
        try:
            shutil.rmtree(megapath_nano.RAM_dir_name)
        except(FileNotFoundError):
            pass


    megapath_nano.log.print('program end')
    megapath_nano.log.print_time()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='MegaPath-Nano: Compositional Analysis')

    parser.add_argument('--debug', dest='debug', default=False, action='store_true')

    # Set up for processing options
    group_adaptor_trimming = parser.add_mutually_exclusive_group(required=False)
    group_adaptor_trimming.add_argument('--adaptor_trimming', dest='adaptor_trimming', action='store_true')
    group_adaptor_trimming.add_argument('--no-adaptor_trimming', dest='adaptor_trimming', action='store_false')

    group_read_trimming = parser.add_mutually_exclusive_group(required=False)
    group_read_trimming.add_argument('--read_trimming', dest='read_trimmingr', action='store_true')
    group_read_trimming.add_argument('--no-read_trimming', dest='read_trimming', action='store_false')

    group_read_filter = parser.add_mutually_exclusive_group(required=False)
    group_read_filter.add_argument('--read_filter', dest='read_filter', action='store_true')
    group_read_filter.add_argument('--no-read_filter', dest='read_filter', action='store_false')

    group_human_filter = parser.add_mutually_exclusive_group(required=False)
    group_human_filter.add_argument('--human_filter', dest='human_filter', action='store_true')
    group_human_filter.add_argument('--no-human_filter', dest='human_filter', action='store_false')

    group_decoy_filter = parser.add_mutually_exclusive_group(required=False)
    group_decoy_filter.add_argument('--decoy_filter', dest='decoy_filter', action='store_true')
    group_decoy_filter.add_argument('--no-decoy_filter', dest='decoy_filter', action='store_false')

    group_variable_region_adjustment = parser.add_mutually_exclusive_group(required=False)
    group_variable_region_adjustment.add_argument('--variable_region_adjustment', dest='variable_region_adjustment', action='store_true')
    group_variable_region_adjustment.add_argument('--no-variable_region_adjustment', dest='variable_region_adjustment', action='store_false')

    group_spike_filter = parser.add_mutually_exclusive_group(required=False)
    group_spike_filter.add_argument('--spike_filter', dest='spike_filter', action='store_true', help='beta version')
    group_spike_filter.add_argument('--no-spike_filter', dest='spike_filter', action='store_false')

    group_closing_spike_filter = parser.add_mutually_exclusive_group(required=False)
    group_closing_spike_filter.add_argument('--closing_spike_filter', dest='closing_spike_filter', action='store_true', help='beta version')
    group_closing_spike_filter.add_argument('--no-closing_spike_filter', dest='closing_spike_filter', action='store_false')

    group_human_repetitive_region_filter = parser.add_mutually_exclusive_group(required=False)
    group_human_repetitive_region_filter.add_argument('--human_repetitive_region_filter', dest='human_repetitive_region_filter', action='store_true', help='beta version')
    group_human_repetitive_region_filter.add_argument('--no-human_repetitive_region_filter', dest='human_repetitive_region_filter', action='store_false')

    group_microbe_repetitive_region_filter = parser.add_mutually_exclusive_group(required=False)
    group_microbe_repetitive_region_filter.add_argument('--microbe_repetitive_region_filter', dest='microbe_repetitive_region_filter', action='store_true', help='beta version')
    group_microbe_repetitive_region_filter.add_argument('--no-microbe_repetitive_region_filter', dest='microbe_repetitive_region_filter', action='store_false')

    group_short_alignment_filter = parser.add_mutually_exclusive_group(required=False)
    group_short_alignment_filter.add_argument('--short_alignment_filter', dest='short_alignment_filter', action='store_true')
    group_short_alignment_filter.add_argument('--no-short_alignment_filter', dest='short_alignment_filter', action='store_false')

    group_unique_alignment = parser.add_mutually_exclusive_group(required=False)
    group_unique_alignment.add_argument('--unique_alignment', dest='unique_alignment', action='store_true')
    group_unique_alignment.add_argument('--no-unique_alignment', dest='unique_alignment', action='store_false')

    group_noise_projection = parser.add_mutually_exclusive_group(required=False)
    group_noise_projection.add_argument('--noise_projection', dest='noise_projection', action='store_true')
    group_noise_projection.add_argument('--no-noise_projection', dest='noise_projection', action='store_false')

    group_similar_species_marker = parser.add_mutually_exclusive_group(required=False)
    group_similar_species_marker.add_argument('--similar_species_marker', dest='similar_species_marker', action='store_true')
    group_similar_species_marker.add_argument('--no-similar_species_marker', dest='similar_species_marker', action='store_false')

    group_mapping_only = parser.add_mutually_exclusive_group(required=False)
    group_mapping_only.add_argument('--mapping_only', dest='mapping_only', action='store_true')
    group_mapping_only.add_argument('--no-mapping_only', dest='mapping_only', action='store_false')

    group_reassign_read_id = parser.add_mutually_exclusive_group(required=False)
    group_reassign_read_id.add_argument('--reassign_read_id', dest='reassign_read_id', action='store_true')
    group_reassign_read_id.add_argument('--no-reassign_read_id', dest='reassign_read_id', action='store_false')
    
    group_filter_fq_only = parser.add_mutually_exclusive_group(required=False)
    group_filter_fq_only.add_argument('--all_taxon_steps', dest='all_taxon_steps', action='store_true')
    group_filter_fq_only.add_argument('--filter_fq_only', dest='all_taxon_steps', action='store_false')

    group_taxon_and_AMR_module_option = parser.add_mutually_exclusive_group(required=False)
    group_taxon_and_AMR_module_option.add_argument('--taxon_and_AMR_module', dest='taxon_and_AMR_module_option', action='store_const',const='taxon_and_AMR_module')
    group_taxon_and_AMR_module_option.add_argument('--taxon_module_only', dest='taxon_and_AMR_module_option', action='store_const',const='taxon_module_only')
    group_taxon_and_AMR_module_option.add_argument('--AMR_module_only', dest='taxon_and_AMR_module_option', action='store_const',const='AMR_module_only')
    
    group_reassignment = parser.add_mutually_exclusive_group(required=False)
    group_reassignment.add_argument('--reassignment', dest='reassignment', action='store_true')
    group_reassignment.add_argument('--no-reassignment', dest='reassignment', action='store_false')

    group_reassignment = parser.add_mutually_exclusive_group(required=False)
    group_reassignment.add_argument('--assembly_selection', dest='assembly_selection', action='store_true')
    group_reassignment.add_argument('--no-assembly_selection', dest='assembly_selection', action='store_false')

    group_reassignment = parser.add_mutually_exclusive_group(required=False)
    group_reassignment.add_argument('--species_level', dest='resolution', action='store_const',const='species')
    group_reassignment.add_argument('--strain_level', dest='resolution', action='store_const',const='strain', help='beta version')
    # Set up for output options

    group_output_adaptor_trimmed_query = parser.add_mutually_exclusive_group(required=False)
    group_output_adaptor_trimmed_query.add_argument('--output_adaptor_trimmed_query', dest='output_adaptor_trimmed_query', action='store_true')
    group_output_adaptor_trimmed_query.add_argument('--no-output_adaptor_trimmed_query', dest='output_adaptor_trimmed_query', action='store_false')

    group_output_trimmed_and_filtered_query = parser.add_mutually_exclusive_group(required=False)
    group_output_trimmed_and_filtered_query.add_argument('--output_trimmed_and_filtered_query', dest='output_trimmed_and_filtered_query', action='store_true')
    group_output_trimmed_and_filtered_query.add_argument('--no-output_trimmed_and_filtered_query', dest='output_trimmed_and_filtered_query', action='store_false')

    group_output_human_and_decoy_filtered_query = parser.add_mutually_exclusive_group(required=False)
    group_output_human_and_decoy_filtered_query.add_argument('--output_human_and_decoy_filtered_query', dest='output_human_and_decoy_filtered_query', action='store_true')
    group_output_human_and_decoy_filtered_query.add_argument('--no-output_human_and_decoy_filtered_query', dest='output_human_and_decoy_filtered_query', action='store_false')

    group_output_PAF = parser.add_mutually_exclusive_group(required=False)
    group_output_PAF.add_argument('--output_PAF', dest='output_PAF', action='store_true')
    group_output_PAF.add_argument('--no-output_PAF', dest='output_PAF', action='store_false')

    group_output_noise_stat = parser.add_mutually_exclusive_group(required=False)
    group_output_noise_stat.add_argument('--output_noise_stat', dest='output_noise_stat', action='store_true')
    group_output_noise_stat.add_argument('--no-output_noise_stat', dest='output_noise_stat', action='store_false')

    group_output_separate_noise_bed = parser.add_mutually_exclusive_group(required=False)
    group_output_separate_noise_bed.add_argument('--output_separate_noise_bed', dest='output_separate_noise_bed', action='store_true')
    group_output_separate_noise_bed.add_argument('--no-output_separate_noise_bed', dest='output_separate_noise_bed', action='store_false')

    group_output_raw_signal = parser.add_mutually_exclusive_group(required=False)
    group_output_raw_signal.add_argument('--output_raw_signal', dest='output_raw_signal', action='store_true')
    group_output_raw_signal.add_argument('--no-output_raw_signal', dest='output_raw_signal', action='store_false')

    group_output_id_signal = parser.add_mutually_exclusive_group(required=False)
    group_output_id_signal.add_argument('--output_id_signal', dest='output_id_signal', action='store_true')
    group_output_id_signal.add_argument('--no-output_id_signal', dest='output_id_signal', action='store_false')

    group_output_per_read_data = parser.add_mutually_exclusive_group(required=False)
    group_output_per_read_data.add_argument('--output_per_read_data', dest='output_per_read_data', action='store_true')
    group_output_per_read_data.add_argument('--no-output_per_read_data', dest='output_per_read_data', action='store_false')

    group_output_quality_score_histogram = parser.add_mutually_exclusive_group(required=False)
    group_output_quality_score_histogram.add_argument('--output_quality_score_histogram', dest='output_quality_score_histogram', action='store_true')
    group_output_quality_score_histogram.add_argument('--no-output_quality_score_histogram', dest='output_quality_score_histogram', action='store_false')

    group_output_read_length_histogram = parser.add_mutually_exclusive_group(required=False)
    group_output_read_length_histogram.add_argument('--output_read_length_histogram', dest='output_read_length_histogram', action='store_true')
    group_output_read_length_histogram.add_argument('--no-output_read_length_histogram', dest='output_read_length_histogram', action='store_false')

    group_output_human_stat = parser.add_mutually_exclusive_group(required=False)
    group_output_human_stat.add_argument('--output_human_stat', dest='output_human_stat', action='store_true')
    group_output_human_stat.add_argument('--no-output_human_stat', dest='output_human_stat', action='store_false')

    group_output_decoy_stat = parser.add_mutually_exclusive_group(required=False)
    group_output_decoy_stat.add_argument('--output_decoy_stat', dest='output_decoy_stat', action='store_true')
    group_output_decoy_stat.add_argument('--no-output_decoy_stat', dest='output_decoy_stat', action='store_false')

    group_output_genome_set = parser.add_mutually_exclusive_group(required=False)
    group_output_genome_set.add_argument('--output_genome_set', dest='output_genome_set', action='store_true')
    group_output_genome_set.add_argument('--no-output_genome_set', dest='output_genome_set', action='store_false')

    # processing options

    parser.set_defaults(adaptor_trimming=True)
    parser.set_defaults(read_trimming=True)
    parser.set_defaults(read_filter=True)
    parser.set_defaults(human_filter=False)
    parser.set_defaults(decoy_filter=False)
    parser.set_defaults(assembly_selection=False)
    parser.set_defaults(variable_region_adjustment=False)
    parser.set_defaults(spike_filter=False)
    parser.set_defaults(human_repetitive_region_filter=False)
    parser.set_defaults(microbe_repetitive_region_filter=False)
    parser.set_defaults(closing_spike_filter=False)
    parser.set_defaults(short_alignment_filter=False)
    parser.set_defaults(unique_alignment=False)
    parser.set_defaults(noise_projection=False)
    parser.set_defaults(similar_species_marker=False)
    parser.set_defaults(reassign_read_id=False)
    parser.set_defaults(all_taxon_steps=True)
    parser.set_defaults(taxon_and_AMR_module_option='taxon_and_AMR_module')
    parser.set_defaults(reassignment=False)
    parser.set_defaults(resolution='species')

    # experimental
    parser.set_defaults(mapping_only=False)

    # options to retain filtered query files
    parser.set_defaults(output_adaptor_trimmed_query=False)
    parser.set_defaults(output_trimmed_and_filtered_query=True)
    parser.set_defaults(output_human_and_decoy_filtered_query=True)

    # output options

    parser.set_defaults(output_PAF=True)
    parser.set_defaults(output_noise_stat=True)
    parser.set_defaults(output_separate_noise_bed=True)
    parser.set_defaults(output_human_stat=True)
    parser.set_defaults(output_decoy_stat=True)
    parser.set_defaults(output_id_signal=True)
    parser.set_defaults(output_raw_signal=True)
    parser.set_defaults(output_per_read_data=True)
    parser.set_defaults(output_quality_score_histogram=True)
    parser.set_defaults(output_read_length_histogram=True)
    parser.set_defaults(output_genome_set=True)


    # Non-user configurable parameters
    
    
    
    CWD=os.path.dirname(os.path.realpath(__file__))
    NANO_DIR=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    parser.add_argument('--temp_folder', help='temporary folder', default='')
    parser.add_argument('--RAM_folder', help='temporary folder in RAM', default='/run/shm')
    parser.add_argument('--taxonomy_db', help='taxonomy database', default=f'{NANO_DIR}/db/ncbi_taxonomy.db')
    parser.add_argument('--tool_folder', help='Tool folder', default=f'{CWD}/tools')
    parser.add_argument('--config_folder', help='Config file folder', default=f'{NANO_DIR}/config')
    parser.add_argument('--assembly_folder', help='Assembly folder', default=f'{NANO_DIR}/genomes')
    parser.add_argument('--db_folder', help='Db folder', default=f'{NANO_DIR}/db')
    parser.add_argument('--aligner', help='Aligner program', default='minimap2')
    parser.add_argument('--read_simulator', help='Read simulation program', default=f'{CWD}/tools/nanosim/simulator.py')
    parser.add_argument('--read_simulation_profiles', help='Read simulation profiles', default=f'{CWD}/tools/nanosim/nanosim_profiles')
    parser.add_argument('--aligner_log', help='Log for stderr output from aligner program', default='minimap2.log')
    parser.add_argument('--read_sim_log', help='Log for stderr output from read simulator', default='read_sim.log')
    parser.add_argument('--adaptor_trimming_log', help='Log for stdout output from adaptor trimming program', default='adaptor_trimming.log')
    parser.add_argument('--python', help='Python entry point', default='python3')

    parser.add_argument('--human_repetitive_region_filter_assembly_id', help='Assembly ID for human similar region filter', default='GCF_000001405.39')

    parser.add_argument('--max_aligner_thread', help='Maximum number of threads used by aligner', type=int, default=64)
    parser.add_argument('--max_aligner_target_GBase_per_batch', help='Maximum size of reference loaded in memory per batch by aligner', default='auto')
    parser.add_argument('--max_porechop_thread', help='Maximum number of threads used by porechop', type=int, default=64)
    parser.add_argument('--max_AMR_thread', help='Maximum number of threads used by AMR module', type=int, default=64)

    parser.add_argument('--genus_height', help='Height in taxonomy to be considered as genus. Full rank info in db_preparation/genAssemblyMetadata.py', type=int, default=11)

    # for all alignments
    parser.add_argument('--min_alignment_score', help='Min alignment score', type=int, default=0)

    # read_trimming_and_filter
    parser.add_argument('--head_crop', help='Number of nucleotides to crop at start of reads', type=int, default=0)
    parser.add_argument('--tail_crop', help='Number of nucleotides to crop at end of reads', type=int, default=0)
    parser.add_argument('--min_read_length', help='Minimum read length', type=int, default=0)
    parser.add_argument('--min_read_quality', help='Minimum average base quality of read', type=float, default=7.0)

    # human_filter, decoy_filter
    parser.add_argument('--human_filter_alignment_score_threshold', help='Alignment score threshold for flagging a read as a human read', type=int, default=1000)
    parser.add_argument('--human_filter_alignment_score_percent_threshold', help='Alignment score (normalized by read length) threshold (in percent) for flagging a read as a human read', type=int, default=100)
    parser.add_argument('--decoy_filter_alignment_score_threshold', help='Alignment score threshold for flagging a read as a decoy read', type=int, default=1000)
    parser.add_argument('--decoy_filter_alignment_score_percent_threshold', help='Alignment score (normalized by read length) threshold (in percent) for flagging a read as a decoy read', type=int, default=100)

    # species_ID
    parser.add_argument('--species_id_min_aligned_bp', help='Min aligned base pairs to include a species for analysis', type=int, default=0)

    # assembly_ID
    parser.add_argument('--good_alignment_threshold', help='Alignment score threshold in percentage of best alignment score', type=int, default=80)
    parser.add_argument('--assembly_id_min_average_depth', help='Min average depth to perform assembly selection (default skip assembly selection)', type=float, default=0.5)

    # variable_region_adjustment
    parser.add_argument('--variable_region_percent', help='Maximum percentage of strands aligned for a region to be labeled as variable', type=int, default=50)

    # spike_filter
    parser.add_argument('--expected_max_depth_stdev', help='Number of standard deviations for calculating expected max depth', type=int, default=6)

    # similar_region_filter
    parser.add_argument('--microbe_repetitive_region_filter_abundance_threshold_80', help='Difference (no. of times) in apparent abundance to trigger similar region filter with 80% similarity', type=float, default=160)
    parser.add_argument('--microbe_repetitive_region_filter_abundance_threshold_90', help='Difference (no. of times) in apparent abundance to trigger similar region filter with 90% similarity', type=float, default=80)
    parser.add_argument('--microbe_repetitive_region_filter_abundance_threshold_95', help='Difference (no. of times) in apparent abundance to trigger similar region filter with 95% similarity', type=float, default=40)
    parser.add_argument('--microbe_repetitive_region_filter_abundance_threshold_98', help='Difference (no. of times) in apparent abundance to trigger similar region filter with 98% similarity', type=float, default=16)
    parser.add_argument('--microbe_repetitive_region_filter_abundance_threshold_99', help='Difference (no. of times) in apparent abundance to trigger similar region filter with 99% similarity', type=float, default=8)
    parser.add_argument('--microbe_repetitive_region_filter_abundance_threshold_99_2', help='Difference (no. of times) in apparent abundance to trigger similar region filter with 99.2% similarity', type=float, default=6.4)
    parser.add_argument('--microbe_repetitive_region_filter_targeted_max_span_percent', help='Maximum percent of regions (targeted) to be marked as similar region', type=int, default=90)
    parser.add_argument('--microbe_repetitive_region_filter_allowed_max_span_percent', help='Maximum percent of regions (allowed) to be marked as similar region', type=int, default=97)
    parser.add_argument('--microbe_repetitive_region_filter_min_average_depth', help='Minimum average depth to be considered as source of noise', type=float, default=0.2)
    # microbe_repetitive_region_filter_max_span_percent_overall is not used currently
    parser.add_argument('--microbe_repetitive_region_filter_max_span_percent_overall', help='Maximum percent of regions to be marked as similar region (overall)', type=int, default=97)

    # noise_removal
    parser.add_argument('--max_alignment_noise_overlap', help='The maximum percent for an alignment to overlap with noise regions without being removed', type=int, default=50)

    # short_alignment_filter
    parser.add_argument('--min_alignment_length', help='Minimum alignment length to be considered as evidence', type=int, default=0)

    # closing_spike_filter
    parser.add_argument('--closing_expected_max_depth_stdev', help='Number of standard deviations for calculating expected max depth for closing spike filter', type=int, default=9)

    # unique_alignment
    parser.add_argument('--unique_alignment_threshold', help='Unique alignments shall have no alignments with alignment score within this percent', type=int, default=80)

    # noise projection
    parser.add_argument('--number_of_genus_to_perform_noise_projection', help='Number of genus to perform noise projection', type=int, default=3)
    parser.add_argument('--min_percent_abundance_to_perform_noise_projection', help='Minimum percent of abundance relative to the most abundant species in a genus to perform noise projection', type=int, default=25)
    parser.add_argument('--noise_projection_simulated_read_length_bin_size', help='Read length bin size for generating simulated reads', type=int, default=1000)
    parser.add_argument('--noise_projection_simulated_read_length_multiplier', help='Multiplier over average read length to obtain maximum read length', type=float, default=0.5)
    parser.add_argument('--noise_projection_simulated_read_error_profile', help='Error profile for generating simulated reads', default='ecoli_R91D')
    parser.add_argument('--noise_projection_num_read_to_simulate', help='Number of simulated reads to generate', type=int, default=10000)

    # similar species marker
    parser.add_argument('--similar_species_marker_num_genus', help='Number of top most abundant species (1 per genus) to be considered as possible source of noise', type=int, default=3)
    parser.add_argument('--similar_species_marker_alignment_similarity_1', help='Similarity cutoff (1) used for alignment', type=int, choices=[99, 98, 95, 90, 80], default=98)
    parser.add_argument('--similar_species_marker_aligned_region_threshold_1', help='Percentage of aligned region (1) to be considered as highly similar', type=int, default=50)
    parser.add_argument('--similar_species_marker_alignment_similarity_2', help='Similarity cutoff (2) used for alignment', type=int, choices=[99, 98, 95, 90, 80], default=95)
    parser.add_argument('--similar_species_marker_aligned_region_threshold_2', help='Percentage of aligned region (2) to be considered as highly similar', type=int, default=75)
    parser.add_argument('--similar_species_marker_similarity_combine_logic', help='Logic for combining criteria 1 and 2 (and / or)', choices=['and', 'or'], default='or')

    # read statistics
    parser.add_argument('--quality_score_bin_size', help='Bin size for quality score histogram', type=float, default=0.2)
    parser.add_argument('--read_length_bin_size', help='Bin size for read length histogram', type=int, default=100)

    # default genome set
    parser.add_argument('--human', help='Human genome set', default='human.genome_set')
    parser.add_argument('--decoy', help='Decoy genome set', default='plasmid.genome_set')
    parser.add_argument('--species', help='Genome set for species identification', default='species_id.genome_set')
    parser.add_argument('--assembly', help='Genome set for assembly identification', default='assembly_id.genome_set')

    # input and output
    parser.add_argument('--query', help='Query file (fastq or fasta)', required=True, default='')
    parser.add_argument('--output_prefix', help='Output prefix', default='')     # by default query file name will be used as output_prefix
    parser.add_argument('--output_folder', help='Output folder', default='./')
    parser.add_argument('--archive_format', help='Format used for archive file', choices=['zip', 'tar', 'gztar', 'bztar'], default='gztar') 

    FLAGS = parser.parse_args()

    main()
