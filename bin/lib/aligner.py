# Calling aligners

import os
import sqlite3
import pandas
import argparse
import tempfile
import subprocess
import shlex
import psutil
import shutil
import hashlib
import random
import pybedtools
from pybedtools.bedtool import BedTool
from bioconvert.sam2paf import SAM2PAF


FLAGS = None
UNPARSED = None


align_list_col_name_no_assembly_id = ['read_id', 'read_length', 'read_from', 'read_to', 'strand', 'sequence_id', 'sequence_length', 'sequence_from', 'sequence_to',
                                      'match', 'mapq', 'edit_dist', 'alignment_score']
align_list_col_type_no_assembly_id = {'read_id': str, 'read_length': int, 'read_from': int, 'read_to': int, 'strand': str, 'sequence_id': str, 'sequence_length': int, 'sequence_from' : int, 'sequence_to': int,
                                      'match': int, 'mapq': int, 'edit_dist': int, 'alignment_score': int}
align_list_col_name = ['read_id', 'read_length', 'read_from', 'read_to', 'strand', 'sequence_id', 'sequence_length', 'sequence_from', 'sequence_to',
                       'match', 'mapq', 'edit_dist', 'alignment_score',
                       'assembly_id', 'tax_id', 'species_tax_id', 'genus_tax_id', 'alignment_score_tiebreaker']
align_list_col_type = {'read_id': str, 'read_length': int, 'read_from': int, 'read_to': int, 'strand': str, 'sequence_id': str, 'sequence_length': int, 'sequence_from' : int, 'sequence_to': int,
                       'match': int, 'mapq': int, 'edit_dist': int, 'alignment_score': int,
                       'assembly_id': str, 'tax_id': int, 'species_tax_id': int, 'genus_tax_id': int, 'alignment_score_tiebreaker': float}
align_stat_col_name_no_assembly_id = ['sequence_id', 'count', 'alignment_score']
align_stat_col_type_no_assembly_id = {'sequence_id': str, 'count': int, 'alignment_score': int}
align_stat_col_name = ['assembly_id', 'count', 'alignment_score']
align_stat_col_type = {'assembly_id': str, 'count': int, 'alignment_score': int}



    
def align_list_to_bed(*,
                      align_list):

    temp_align_list = align_list.assign(assembly_id_sequence_id = lambda x: f"{x['assembly_id']},{x['sequence_id']}")
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


# query_assembly_list and target_assembly_list are namedtuple with key 'assembly_id'
# query and target are list of filenames of fasta or fastq file with full path

def Align(*, 
          assembly_metadata, 
          global_options, 
          temp_dir_name,
          log_file,
          query_filename_list=None, 
          query_assembly_list=None, 
          target_filename_list=None, 
          target_assembly_list=None, 
          aligner_options=None,
          paf_path_and_prefix=None,
          mapping_only=False,
          module_option='',
          AMR_output_folder='',
          align_concat_fa=False,
          ):
    NANO_DIR=os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
    local_temp_dir_name = tempfile.mkdtemp(prefix='Align.', dir=temp_dir_name)

    if align_concat_fa==False:
        num_target_specification = 0
        if target_filename_list is not None and target_filename_list['path'].shape[0] > 0:
            num_target_specification += 1
        if target_assembly_list is not None and target_assembly_list['assembly_id'].shape[0] > 0:
            num_target_specification += 1
        if num_target_specification != 1:
            os.sys.exit('Exactly one of target_filename_list and target_assembly_list must be specified')

        #target_filename_list has only contig-level genomes, target_assembly_* has all
        if target_assembly_list is not None and len(target_assembly_list['assembly_id']) > 0:
            target_filename_list = assembly_metadata.get_assembly_path(assembly_list=target_assembly_list)
            if target_filename_list is None:
                os.sys.exit('Target assembly_id not found')
            target_filename_list['path'] = target_filename_list['path'].map(lambda x: os.path.join(global_options['assembly_folder'], x))
        
        if target_assembly_list is not None and len(target_assembly_list['assembly_id']) > 0:
            target_assembly_length = assembly_metadata.get_assembly_length(assembly_list=target_assembly_list)
            if target_assembly_length is None: 
                os.sys.exit('Target assembly_length not found')
        else:
            target_assembly_length = target_filename_list.assign(assembly_length = lambda x: 1)
    
    
    
        num_target = target_filename_list.shape[0]
        if target_assembly_length.shape[0] != num_target:
            os.sys.exit('Number of target_assembly_length does not match number of target_filename_list') 

        
        temp_pipe_target_fasta_name = os.path.join(local_temp_dir_name, 'temp_pipe_target_fasta')
        os.mkfifo(temp_pipe_target_fasta_name)

    num_query_specification = 0
    if query_filename_list is not None and query_filename_list['path'].shape[0] > 0:
        num_query_specification += 1
    if query_assembly_list is not None and query_assembly_list['assembly_id'].shape[0] > 0:
        num_query_specification += 1

    if num_query_specification != 1:
        os.sys.exit('Exactly one of query_filename_list and query_assembly_list must be specified')
    if query_assembly_list is not None and query_assembly_list['assembly_id'].shape[0] > 0:
        query_filename_list = assembly_metadata.get_assembly_path(assembly_list=query_assembly_list)
        if query_filename_list is None or query_filename_list.shape[0] != query_assembly_list['assembly_id'].shape[0]:
            os.sys.exit('Query assembly_id not found')
        query_filename_list['path'] = query_filename_list['path'].map(lambda x: os.path.join(global_options['assembly_folder'], x))

    random_hash_string = ''
    
    for query in (query_filename_list['path']):
        if os.path.isfile(query) == False:
            os.sys.exit('Query file ' + query + ' not exists')
        random_hash_string = random_hash_string + os.path.split(query)[1]
    
    random_hash_object = hashlib.md5(random_hash_string.encode())
    random.seed(random_hash_object.hexdigest())



    read_length = None
    align_list = pandas.DataFrame(columns=align_list_col_name_no_assembly_id)
    for key, value in align_list_col_type_no_assembly_id.items():
        align_list[key] = align_list[key].astype(value)

    if paf_path_and_prefix is None or paf_path_and_prefix == '':
        output_PAF = False
    else:
        output_PAF = True
    
    if output_PAF == True:
        paf_filename = f'{paf_path_and_prefix}.paf'
        sam_filename = f'{paf_path_and_prefix}.sam'

    #  align query to target (input from named pipe)
    aligner_command = [global_options['aligner'],]
    if mapping_only == False:
        aligner_command.extend(('-c',))
    
    if output_PAF == True:
        aligner_command.extend(('-a',))
        #aligner_command.extend(['--split-prefix','tmp'])
    aligner_command.extend(aligner_options)

    if align_concat_fa==False:    
        aligner_command.append(temp_pipe_target_fasta_name)
        aligner_command.extend(query_filename_list['path'])
        aligner_command.extend(['--split-prefix','tmp'])
        if output_PAF == True:
            sam_file = os.open(sam_filename, flags=os.O_CREAT |  os.O_WRONLY|os.O_TRUNC, mode=0o644)
            aligner_stdout=sam_file
        else:
            aligner_stdout=subprocess.PIPE
        
        aligner_process = subprocess.Popen(aligner_command, close_fds=True, stdout=aligner_stdout, stderr=log_file)
        
        #  cat target to pipe
        temp_pipe_target_fasta = os.open(temp_pipe_target_fasta_name, flags=os.O_WRONLY)
        
        #chunk
        chunk_size=1000
        for i in range(0,len(target_filename_list['path']),chunk_size):
            cat_command = ['cat',]
            cat_command.extend(target_filename_list['path'][i:i+chunk_size].tolist())
            cat_process = subprocess.Popen(cat_command, close_fds=True, stdout=temp_pipe_target_fasta)
            cat_process.wait()

    elif align_concat_fa==True:
        if module_option!='amplicon_filter_module':
            aligner_command.append(f'{NANO_DIR}/genomes/refseq/refseq.fna.gz')
        else:
            aligner_command.append(f'{target_filename_list["path"][0]}')
        aligner_command.extend(query_filename_list['path'])
        if output_PAF == True:# and module_option!='amplicon_filter_module':
            sam_file = os.open(sam_filename, flags=os.O_CREAT |  os.O_WRONLY|os.O_TRUNC, mode=0o644)
            aligner_stdout=sam_file
        else:
            aligner_stdout=subprocess.PIPE
        
        aligner_process = subprocess.Popen(aligner_command, close_fds=True, stdout=aligner_stdout, stderr=log_file)
        #if module_option=='amplicon_filter_module':
        #    bam_output = subprocess.check_output(('samtools', 'sort','-o',f'{paf_path_and_prefix}.bam'), stdin=aligner_process.stdout)
        #    aligner_process.wait()
        #    index_process = subprocess.Popen(('samtools','index',f'{paf_path_and_prefix}.bam'))
        #    index_process.wait()
        #    os.sys.exit('Finished amplicon filter module.')

    if output_PAF == True:
        if align_concat_fa==False:
            cat_process.wait()
            os.close(temp_pipe_target_fasta)
        aligner_process.wait()
        print('Finished species alignment step.')
        bam_filename=f'{paf_path_and_prefix}.bam'
        exclude_flag='1796'
        if module_option=='amplicon_filter_module':
            exclude_flag='4'
        bam_operation_command = f'samtools view -F{exclude_flag} -b {sam_filename}|samtools sort -o {bam_filename};samtools index {bam_filename};'
        bam_operation_process = subprocess.Popen(bam_operation_command, shell=True, stderr=subprocess.DEVNULL)
        if module_option in ['taxon_and_AMR_module','AMR_module_only']:
            bam_operation_process.wait()
            run_amr_command=f'python {NANO_DIR}/megapath_nano_amr.py --query_bam {bam_filename} --output_folder {AMR_output_folder} --threads {global_options["AMRThreadOption"]}'
            run_amr_process=subprocess.Popen(run_amr_command, shell=True)
        if module_option in ['AMR_module_only','amplicon_filter_module']:
            bam_operation_process.wait()
            os.sys.exit()
        convert=SAM2PAF(sam_filename,paf_filename)
        convert(pri_only=False)
        awk_stdin= os.open(paf_filename, flags=os.O_RDONLY)
    elif output_PAF == False:
        awk_stdin = aligner_process.stdout

    # step 3: retain the required columns only (connect from step 2 through stdout-stdin)
    
    aligner_output_filename = os.path.join(local_temp_dir_name, 'align_list')
    aligner_output = os.open(aligner_output_filename, flags=os.O_CREAT | os.O_EXCL | os.O_WRONLY, mode=0o644)
    
    awk_command = ['awk', ]
    awk_command.append('{OFS="\t"};{gsub("NM:i:","")};{gsub("AS:i:","")};{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$15}')
    awk_process = subprocess.Popen(awk_command, close_fds=True, stdin=awk_stdin, stdout=aligner_output)


    if output_PAF == False:
        cat_process.wait()
        os.close(temp_pipe_target_fasta)
        
        aligner_process.wait()

    
    awk_process.wait()
    os.close(aligner_output)
    if output_PAF == True:
        os.close(sam_file)


    # step 3: Filter for unaligned alignment

    list_col = ('read_id', 'read_length', 'read_from', 'read_to', 'strand', 'sequence_id', 'sequence_length', 'sequence_from', 'sequence_to',
                'match', 'alignment_block_length', 'mapq', 'edit_dist', 'alignment_score')
    col_type = {'read_id': str, 'read_length': int, 'read_from': int, 'read_to': int, 'strand': str, 'sequence_id': str, 'sequence_length': int,
                'sequence_from' : int, 'sequence_to': int, 'match': int, 'alignment_block_length': int, 'mapq': int, 'edit_dist': int, 'alignment_score': int}

    try:
        prefilter_align_list = pandas.read_csv(aligner_output_filename, 
                                               sep='\t', names=list_col, index_col=False, dtype=col_type, header=None
                                              )
    except pandas.errors.ParserError:
        prefilter_align_list = pandas.DataFrame(columns=list_col)
        for key, value in col_type.items():
            prefilter_align_list[key] = prefilter_align_list[key].astype(value)

    if global_options['debug'] == False:
        os.remove(aligner_output_filename)


    # step 4: Filter alignments

    min_alignment_score = global_options['min_alignment_score']
    align_list = prefilter_align_list.query('sequence_length > 0 and alignment_score >= @min_alignment_score')


    # step 5: Join align_list with sequence_assembly_tax_id

    if target_assembly_list is not None and len(target_assembly_list['assembly_id']) > 0:
        sequence_assembly_tax_id = assembly_metadata.get_sequence_tax_id(assembly_list=target_assembly_list).set_index(['sequence_id'])[['assembly_id', 'tax_id', 'species_tax_id', 'genus_tax_id']]

        num_align = align_list.shape[0]
        align_list = align_list.merge(
                                      right=sequence_assembly_tax_id,
                                      how='inner',
                                      left_on='sequence_id',
                                      right_index = True,
                                      suffixes=['', '_y'],
                                      validate='m:1'
                                     )

        
        if align_list.shape[0] != num_align:
            print('Some sequence id cannot be matched', file=os.sys.stderr)

    align_list = align_list.assign(alignment_score_tiebreaker = lambda x: 0)
    align_list['alignment_score_tiebreaker'] = align_list['alignment_score_tiebreaker'].apply(lambda x: random.random())

    # step 6: Return align_list and read_info

    if global_options['debug'] == False:
        shutil.rmtree(local_temp_dir_name)

    return align_list



def main():

    pandas.set_option('display.width', 120)
    pandas.set_option('display.max_colwidth', 100)

    num_target_specified = sum(map(lambda x: 0 if x == '' else 1, (FLAGS.target, FLAGS.target_list, FLAGS.target_assembly_id, )))
    num_query_specified = sum(map(lambda x: 0 if x == '' else 1, (FLAGS.query, FLAGS.query_list, FLAGS.query_assembly_id, )))

    if num_target_specified != 1:
        os.sys.exit('Exactly one of --target, --target_list, and --target_assembly_id must be specified')

    if num_query_specified != 1:
        os.sys.exit('Exactly one of --query, --query_list, and --query_assembly_id must be specified')


    global_options = {
                      'tool_folder':FLAGS.tool_folder,
                      'assembly_folder':FLAGS.assembly_folder,
                      'aligner':FLAGS.aligner,
                      'aligner_log':FLAGS.aligner_log,
                      'max_target_GBase_per_batch':FLAGS.max_target_GBase_per_batch,
                      'min_alignment_score':FLAGS.min_alignment_score,
                      'debug':FLAGS.debug,
                     }

    target_filename_list = None
    if FLAGS.target:
        target_filename_list = pandas.DataFrame(FLAGS.target.split(','), columns=['path'])

    target_assembly_list = None
    if FLAGS.target_assembly_id:
        target_assembly_list = pandas.DataFrame(FLAGS.target_assembly_id.split(','), columns=['assembly_id'])
    if FLAGS.target_list:
        target_assembly_list = pandas.read_csv(FLAGS.target_list, sep='\t', names=('assembly_id',), usecols=(0,), header=None)

    query_filename_list = None
    if FLAGS.query:
        query_filename_list = pandas.DataFrame(FLAGS.query.split(','), columns=['path'])

    query_assembly_list = None
    if FLAGS.query_assembly_id:
        query_assembly_list = pandas.DataFrame(FLAGS.query_assembly_id.split(','), columns=['assembly_id'])
    if FLAGS.query_list:
        query_assembly_list = pandas.read_csv(FLAGS.query_list, sep='\t', names=('assembly_id',), usecols=(0,), header=None)


    assembly_metadata = AssemblyMetadata(assembly_folder=FLAGS.assembly_folder)  

    pybedtools.helpers.set_bedtools_path(os.path.join(FLAGS.tool_folder, 'bedtools'))

    if FLAGS.temp_folder != '':
        temp_dir_name = tempfile.mkdtemp(prefix='aligner.', dir=FLAGS.temp_folder)
    else:
        temp_dir_name = tempfile.mkdtemp(prefix='aligner.')

    if FLAGS.RAM_folder != '':
        pybedtools.helpers.set_tempdir(FLAGS.RAM_folder)
 
    if global_options['aligner_log'] == '':
        aligner_log = os.sys.stderr
    else:
        aligner_log = os.open(global_options['aligner_log'], flags=os.O_CREAT | os.O_WRONLY, mode=0o644)

    align_list = Align(
                       assembly_metadata=assembly_metadata, 
                       global_options=global_options,
                       temp_dir_name=temp_dir_name,
                       log_file=aligner_log,
                       query_filename_list=query_filename_list,
                       query_assembly_list=query_assembly_list,
                       target_filename_list=target_filename_list,
                       target_assembly_list=target_assembly_list,
                       aligner_options=UNPARSED,
                       paf_path_and_prefix=FLAGS.paf_prefix,
                       mapping_only=FLAGS.mapping_only,
                      )

    if FLAGS.output_format == 'best' or FLAGS.output_format == 'stat':
        best_align_list = align_list.sort_values(['read_id', 'alignment_score', 'alignment_score_tiebreaker']).drop_duplicates(subset=['read_id'], keep='last')
        if FLAGS.output_format == 'best':
            if FLAGS.target:
                best_align_list.to_csv(path_or_buf=os.sys.stdout, sep='\t', header=False, index=False, columns=align_list_col_name_no_assembly_id)
            else:
                best_align_list.to_csv(path_or_buf=os.sys.stdout, sep='\t', header=False, index=False, columns=align_list_col_name)
        else:
            if FLAGS.target:
                best_align_stat = best_align_list.assign(count = lambda x: 1).groupby(['sequence_id'], as_index=False).sum()
                best_align_stat.to_csv(path_or_buf=os.sys.stdout, sep='\t', na_rep='NaN', header=False, index=False, columns=align_stat_col_name_no_assembly_id)
            else:
                best_align_stat = best_align_list.assign(count = lambda x: 1).groupby(['assembly_id'], as_index=False).sum()
                best_align_stat.to_csv(path_or_buf=os.sys.stdout, sep='\t', na_rep='NaN', header=False, index=False, columns=align_stat_col_name)
    elif FLAGS.output_format == 'bed' or FLAGS.output_format == 'bp':
        if FLAGS.target:
            align_list.assign(assembly_id = lambda x: '')
        bed = align_list_to_bed(align_list=align_list)
        if FLAGS.output_format == 'bed':
            print(bed)
        else:
            if FLAGS.target:
                covered_bp = bed_to_covered_bp_by_sequence_id(bed=bed)
                covered_bp.to_csv(path_or_buf=os.sys.stdout, sep='\t', na_rep='NaN', header=False, index=False, columns=['sequence_id', 'covered_bp'])
            else:
                covered_bp = bed_to_covered_bp_by_assembly_id(bed=bed)
                covered_bp.to_csv(path_or_buf=os.sys.stdout, sep='\t', na_rep='NaN', header=False, index=False, columns=['assembly_id', 'covered_bp'])
        os.remove(bed.fn)
    elif FLAGS.output_format == 'all':
        if FLAGS.target:
            align_list.to_csv(path_or_buf=os.sys.stdout, sep='\t', na_rep='NaN', header=False, index=False, columns=align_list_col_name_no_assembly_id)
        else:
            align_list.to_csv(path_or_buf=os.sys.stdout, sep='\t', na_rep='NaN', header=False, index=False, columns=align_list_col_name)
    else:
        pass
    

    if FLAGS.debug == False:
        shutil.rmtree(temp_dir_name)


if __name__ == '__main__':
    NANO_DIR=os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
    parser = argparse.ArgumentParser(description='Aligner')

    parser.add_argument('--debug', dest='debug', default=False, action='store_true')

    parser.add_argument('--mapping_only', dest='mapping_only', default=False, action='store_true')

    parser.add_argument('--temp_folder', help='temporary folder', default='')
    parser.add_argument('--RAM_folder', help='temporary folder in RAM', default='/run/shm')

    parser.add_argument('--tool_folder', help='Tool folder', default=f'{NANO_DIR}/tools')
    parser.add_argument('--assembly_folder', help='Assembly folder', default=f'{NANO_DIR}/genomes')
    parser.add_argument('--aligner', help='Aligner program', default='minimap2')
    parser.add_argument('--aligner_log', help='Log for stderr output from aligner program', default='')

    group_output_format = parser.add_mutually_exclusive_group(required=False)
    group_output_format.add_argument('--best', dest='output_format', action='store_const', const='best')
    group_output_format.add_argument('--all', dest='output_format', action='store_const', const='all')
    group_output_format.add_argument('--bed', dest='output_format', action='store_const', const='bed')
    group_output_format.add_argument('--bp', dest='output_format', action='store_const', const='bp')
    group_output_format.add_argument('--stat', dest='output_format', action='store_const', const='stat')

    parser.set_defaults(output_format='all')

    parser.add_argument('--max_target_GBase_per_batch', help='Max number of target GBases per batch', type=int, default=8)
    parser.add_argument('--min_alignment_score', help='Min alignment score', type=int, default=0)

    parser.add_argument('--target', help='Target file (fasta)', default='')
    parser.add_argument('--target_list', help='Target file (list of assembly ID)', default='')
    parser.add_argument('-t', '--target_assembly_id', help='Target assembly ID (comma-delimited)', default='')

    parser.add_argument('--query', help='Query file (fastq or fasta)', default='')
    parser.add_argument('--query_list', help='Query file (list of assembly ID)', default='')
    parser.add_argument('-q', '--query_assembly_id', help='Query assembly ID (comma-delimited)', default='')

    parser.add_argument('--paf_prefix', help='Prefix of PAF file output', default='')
   

    FLAGS, UNPARSED = parser.parse_known_args()

    os.sys.path.append(os.getcwd())
    from assembly_metadata import AssemblyMetadata


    main()
