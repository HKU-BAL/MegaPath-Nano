import sys
import re
import os
import json
import random
import shlex
from os.path import isfile
from argparse import ArgumentParser
from sys import stderr
from subprocess import PIPE, Popen
from collections import defaultdict

import fast_align_reads2ref as fast_align

#parameters
cigarRe = r"(\d+)([MIDNSHP=X])"
random.seed(0)
min_mapping_quality = 20
min_base_quality = 20
calling_min_mapping_quality = 0
calling_min_base_quality = 0
expandReferenceRegion = 1000000
no_of_positions = 33
realign_windows = 200

class Position(object):
    def __init__(self, pos, ref_base=None, alt_base=None, genotype=None, phase_set=None):
        self.pos = pos
        self.ref_base = ref_base
        self.alt_base = alt_base
        self.read_info = {}
        self.name_base_dict = {}
        self.genotype = genotype
        self.phase_set = phase_set
        self.read_name_seq = defaultdict(str)
        self.ref_seq = None
        self.alt_seq = None

    def set_name_base_dict(self, base_list, read_name_list, base_quality_list=None):
        self.base_list = base_list
        self.read_name_list = read_name_list
        self.read_name_dict = dict(zip(read_name_list, base_list))
        if base_quality_list:
            self.bq_dict = dict(zip(read_name_list, base_quality_list))


def byte(x):
    return bytes(x, encoding="utf8")

def get_reference_seq(sequence, start ,end, reference_start_0_based):
    #0 base
    if end < start:
        end, start = start, end
    return sequence[start - reference_start_0_based: end - reference_start_0_based]

def subprocess_popen(args, stdin=None, stdout=PIPE, stderr=stderr, bufsize=8388608):
    return Popen(args, stdin=stdin, stdout=stdout, stderr=stderr, bufsize=bufsize, universal_newlines=True)


def region_from(ctg_name, ctg_start=None, ctg_end=None):
    """
    1-based region string [start, end]
    """
    if ctg_name is None:
        return ""
    if (ctg_start is None) != (ctg_end is None):
        return ""

    if ctg_start is None and ctg_end is None:
        return "{}".format(ctg_name)
    return "{}:{}-{}".format(ctg_name, ctg_start, ctg_end)

def get_nearby_reads(center_pos, pileup_dict):
    all_nearby_read_name_dict = defaultdict(str)
    all_nearby_bq_dict = defaultdict(str)
    read_name_start = defaultdict()
    start_pos, end_pos = center_pos - realign_windows, center_pos + realign_windows + 1
    for p in range(start_pos, end_pos):
        if p not in pileup_dict:
            if p > center_pos:
                break
            continue
        nearby_read_name = pileup_dict[p].read_name_list
        for read_name in nearby_read_name:
            if read_name not in all_nearby_read_name_dict:
                read_name_start[read_name] = p
            base, indel = pileup_dict[p].read_name_dict[read_name]
            if base in '#*': # skip Nn
                continue
            if base.upper() in "ACGT":
                all_nearby_read_name_dict[read_name] += base.upper()
                if len(pileup_dict[p].bq_dict[read_name]) != 1:
                    print (1)
                all_nearby_bq_dict[read_name] += pileup_dict[p].bq_dict[read_name]
            if indel != '' and indel[0] == '+':
                all_nearby_read_name_dict[read_name] += indel[1:].upper()
                all_nearby_bq_dict[read_name] += "A" * len(indel[1:])
    for read_name in all_nearby_bq_dict:
        if len(all_nearby_bq_dict[read_name]) != len(all_nearby_read_name_dict[read_name]):
            print ('not match')
    return all_nearby_read_name_dict, read_name_start, all_nearby_bq_dict

def decode_pileup_bases(pileup_bases):
    base_idx = 0
    base_list = []
    base_offset = 0
    while base_idx < len(pileup_bases):
        base = pileup_bases[base_idx]
        if base =='+' or base == '-':
            base_idx += 1
            advance = 0
            while True:
                num = pileup_bases[base_idx]
                if num.isdigit():
                    advance = advance * 10 + int(num)
                    base_idx += 1
                else:
                    break
            base_list[-1][1] = base + pileup_bases[base_idx: base_idx + advance]
            base_idx += advance - 1

        elif base in "ACGTNacgtn#*":
            base_list.append([base, ""])
            base_offset += 1
        elif base =='^': # start of read, next base is mq
            base_idx += 1
        base_idx += 1
    return base_list

def reference_sequence_from(samtools_execute_command, fasta_file_path, regions):
    refernce_sequences = []
    region_value_for_faidx = " ".join(regions)

    samtools_faidx_process = subprocess_popen(
        shlex.split("{} faidx {} {}".format(samtools_execute_command, fasta_file_path, region_value_for_faidx))
    )
    while True:
        row = samtools_faidx_process.stdout.readline()
        is_finish_reading_output = row == '' and samtools_faidx_process.poll() is not None
        if is_finish_reading_output:
            break
        if row:
            refernce_sequences.append(row.rstrip())

    # first line is reference name ">xxxx", need to be ignored
    reference_sequence = "".join(refernce_sequences[1:])

    # uppercase for masked sequences
    reference_sequence = reference_sequence.upper()

    samtools_faidx_process.stdout.close()
    samtools_faidx_process.wait()
    if samtools_faidx_process.returncode != 0:
        return None

    return reference_sequence


def count_align_length(cigar, match=2,mismatch=6, gap_open=8, gap_extend=2):
    length = 0
    new_cigar = ""
    soft_clip_start, soft_clip_end, found_start = 0, 0 , False
    for index, m in enumerate(re.finditer(cigarRe, cigar)):
        l, op = int(m.group(1)), m.group(2)
        if op in 'MX=':
            length += int(l)
        elif op == 'I':
            length += int(l)
        elif op == 'S':
            if index == 0:
                soft_clip_start = int(l)
            else:
                soft_clip_end = int(l)
        if op != "S":
            new_cigar += str(l) + op
    return length, soft_clip_start, soft_clip_end, new_cigar

def local_realignment(args):
    ctg_start = args.ctgStart
    ctg_end = args.ctgEnd
    position = args.pos
    bed_file_path = args.bed_fn
    fasta_file_path = args.ref_fn
    ctg_name = args.ctgName
    samtools_execute_command = args.samtools
    ilmn_realign_folder = args.ilmn_realign_folder
    bam_file_path = args.bam_fn
    platform = args.platform
    is_illumina_pltform = platform == 'illumina'
    pos = int(position)
    test_pos = pos
    need_phasing_pos_set = set([pos])
    ctg_start = pos - realign_windows
    ctg_end = pos + realign_windows + 1

    if is_illumina_pltform and ilmn_realign_folder is not None:
        # replace bam with realigned bam fn
        fn = bam_file_path.split('/')
        directry, bam_fn = '/'.join(fn[:-1]), fn[-1]
        bam_file_path = os.path.join(ilmn_realign_folder, bam_fn)
        if not os.path.exists(bam_file_path) or not os.path.exists(bam_file_path + '.bai'):
            exit('[ERROR] Illumina realigned folder not found')
    # if args.test_illumina:
    #     _, pos = position.split('^')
    #     ctg_name = 'NC_000962.3'
    #     pos = int(pos)
    #     test_pos = pos
    #     need_phasing_pos_set = set([pos])
    #     ctg_start = pos - realign_windows
    #     ctg_end = pos + realign_windows + 1
    # elif position:
    #     bam, pos = position.split('^')
    #     ctg_name = 'NC_000962.3'
    #     pos = int(pos)
    #     test_pos = pos
    #     need_phasing_pos_set = set([pos])
    #     ctg_start = pos - realign_windows
    #     ctg_end = pos + realign_windows + 1

    if not isfile("{}.fai".format(fasta_file_path)):
        print("Fasta index {}.fai doesn't exist.".format(fasta_file_path), file=sys.stderr)
        sys.exit(1)

    # 1-based regions [start, end] (start and end inclusive)
    ref_regions = []
    reads_regions = []

    is_ctg_name_given = ctg_name is not None
    is_ctg_range_given = is_ctg_name_given and ctg_start is not None and ctg_end is not None

    if is_ctg_range_given:
        reads_regions.append(region_from(ctg_name=ctg_name, ctg_start=ctg_start - no_of_positions, ctg_end=ctg_end + no_of_positions))
        reference_start, reference_end = ctg_start - expandReferenceRegion, ctg_end + expandReferenceRegion
        reference_start = 1 if reference_start < 1 else reference_start
        ref_regions.append(region_from(ctg_name=ctg_name, ctg_start=reference_start, ctg_end=reference_end))
    elif is_ctg_name_given:
        reads_regions.append(region_from(ctg_name=ctg_name))
        ref_regions.append(region_from(ctg_name=ctg_name))
        reference_start = 1

    reference_sequence = reference_sequence_from(
        samtools_execute_command=samtools_execute_command,
        fasta_file_path=fasta_file_path,
        regions=ref_regions
    )
    if reference_sequence is None or len(reference_sequence) == 0:
        print("[ERROR] Failed to load reference seqeunce from file ({}).".format(fasta_file_path), file=sys.stderr)
        sys.exit(1)
    if bed_file_path:
        samtools_mpileup_process = subprocess_popen(
        shlex.split(
            "{} mpileup  {} -r {} -l {} --reverse-del --min-MQ {} --min-BQ {} --output-QNAME".format(samtools_execute_command, bam_file_path,
                                                                                " ".join(reads_regions),
                                                                                bed_file_path,
                                                                                calling_min_mapping_quality,
                                                                            calling_min_base_quality)))
    else:
        samtools_mpileup_process = subprocess_popen(
        shlex.split(
            "{} mpileup  {} -r {} --reverse-del --min-MQ {} --min-BQ {} --output-QNAME".format(samtools_execute_command, bam_file_path,
                                                                                " ".join(reads_regions),
                                                                                calling_min_mapping_quality,
                                                                                calling_min_base_quality)))  # without -F to use more reads for local-reassembly


    output_bam = args.output_bam
    pileup_dict = defaultdict(str)
    for row in samtools_mpileup_process.stdout: # chr postion N depth seq BQ read_name
        columns = row.strip().split('\t')
        pos = int(columns[1])
        pileup_bases = columns[4]
        read_name_list = columns[6].split(',')
        base_list = decode_pileup_bases(pileup_bases)
        base_quality_list = columns[5]
        if len(read_name_list) != len(base_list):
            continue
        if len(base_quality_list) != len(base_list):
            continue
        pileup_dict[pos] = Position(pos=pos)
        pileup_dict[pos].set_name_base_dict(base_list, read_name_list, base_quality_list)


    if is_illumina_pltform and test_pos not in pileup_dict:
        print (['ERROR pos {} has no support after realignment'.format(test_pos)])

    elif is_illumina_pltform:
        alt_dict = defaultdict(int)
        base_list = pileup_dict[test_pos].base_list
        depth = 0
        reference_base = reference_sequence[test_pos - reference_start].upper()
        for base, indel in base_list:
            if base in "#*":
                depth += 1
                continue
            depth += 1
            if indel != '':
                if indel[0] == '+':
                    alt_dict['I' + indel[1:].upper()] += 1
                else:  # del
                    alt_dict['D' + indel.upper()] += 1
            elif base.upper() != reference_base:
                alt_dict[base.upper()] += 1
        output = open(args.output_fn, 'w')
        output.write(
            '\t'.join([ctg_name, str(test_pos), str(depth), json.dumps(alt_dict), position]) + '\n')
        output.close()
        return


    header = subprocess_popen(
        shlex.split(
            "{} view -H {}".format(samtools_execute_command,bam_file_path))).stdout.read()
    output = open(args.output_fn, 'w')
    for pos_idx, pos in enumerate(sorted(list(need_phasing_pos_set))):
        all_read_name_dict, read_name_start, all_nearby_bq_dict = get_nearby_reads(pos, pileup_dict)
        all_read_name = list(all_read_name_dict.keys())
        all_seq_list = [all_read_name_dict[item] for item in all_read_name]
        ref_start = pos - realign_windows - 100
        ref_end = pos + realign_windows + 100 + 1
        ref_seq = reference_sequence[ref_start - reference_start: ref_end-reference_start].upper()
        my_align = fast_align.FastPassAligner()
        my_align.reference = ref_seq
        my_align.reference_start = ref_start
        my_align.reference_end = ref_end
        my_align.consensus = all_seq_list
        my_align.read_name_list = all_read_name
        align_info = sorted(my_align.align_reads(), key= lambda x: x[2])
        alt_info = defaultdict(int)
        if len(all_read_name) != len(align_info):
            print ('Read name length not equal, skip')
            continue
        output_bam_fn = "{}/{}.bam".format(os.path.dirname(args.output_fn), pos)
        if output_bam:
            save_file_fp = subprocess_popen(shlex.split("{} view -bh - -o {}".format(samtools_execute_command, output_bam_fn)), stdin=PIPE,stdout=PIPE)
            save_file_fp.stdin.write(header)

        output_read = []
        for SEQ, CIGAR, reference_position, read_name in align_info:
            if read_name == "0":
                continue
            base_quality = all_nearby_bq_dict[read_name]
            cigar_length,soft_clip_start,soft_clip_end, new_cigar = count_align_length(CIGAR)
            SEQ1 = SEQ

            SEQ = SEQ[:-soft_clip_end] if soft_clip_end > 0 else SEQ
            base_quality = base_quality[:-soft_clip_end] if soft_clip_end > 0 else base_quality
            if len(base_quality) != len(SEQ):
                print (read_name, len(base_quality), len(SEQ), len(SEQ1), CIGAR, new_cigar)
            new_rp = reference_position
            if output_bam:
                if cigar_length != len(SEQ):
                    print(CIGAR, len(SEQ), len(SEQ1), cigar_length, soft_clip_start,soft_clip_end)
                    continue
                if soft_clip_start > 0:
                    new_cigar = str(soft_clip_start) + 'S' + new_cigar
                    SEQ += "N" * soft_clip_start
                    base_quality += "N" * soft_clip_start
                output_read.append([read_name, "0", "NC_000962.3", str(new_rp),
                                      str(60), new_cigar, "*", "1", "1",
                                      SEQ,
                                      base_quality])
            advance = 0
            query_position = 0

            for c in str(new_cigar):

                if c.isdigit():
                    advance = advance * 10 + int(c)
                    continue
                # soft clip
                if c == "S":
                    query_position += advance
                # match / mismatch
                if c in 'M=X':
                    if pos >= reference_position and reference_position + advance > pos:

                        offset = pos - reference_position
                        if query_position + offset < len(SEQ):
                            alt_info[SEQ[query_position + offset]] += 1
                    reference_position += advance
                    query_position += advance
                # insertion
                if c == "I":
                    if pos == reference_position - 1:
                        ins_seq = SEQ[query_position:query_position + advance]
                        alt_info['I'+ins_seq.upper()] += 1
                    query_position += advance

                # deletion
                if c == "D":
                    if pos == reference_position - 1:
                        del_seq = reference_sequence[
                                  reference_position - reference_start: reference_position - reference_start + advance]
                        alt_info['D'+del_seq] += 1
                    reference_position += advance
                # reset advance
                advance = 0
        reference_base = reference_sequence[pos - reference_start].upper()
        depth = 0
        for alt_type, count in alt_info.items():
            if alt_type[0] not in 'ID':
                depth += count
        if reference_base in alt_info:
            del alt_info[reference_base]

        output.write( '\t'.join([ctg_name, str(pos), str(depth), json.dumps(alt_info), position]) + '\n')

        if output_bam:
            output_read = sorted(output_read, key=lambda x: int(x[3]))

            read_str = ['\t'.join(item) for item in output_read]
            save_file_fp.stdin.write('\n'.join(read_str) + '\n')

            save_file_fp.stdin.close()
            save_file_fp.wait()

            save_index_fp = subprocess_popen(shlex.split("samtools index {}".format(output_bam_fn)),
                stdin=PIPE, stdout=PIPE)
            save_index_fp.stdin.close()
            save_index_fp.wait()

    output.close()

    samtools_mpileup_process.stdout.close()
    samtools_mpileup_process.wait()



def main():
    parser = ArgumentParser(description="Generate 1-based variant candidates using alignments")

    parser.add_argument('--bam_fn', type=str, default="input.bam",
                        help="Sorted bam file input, default: %(default)s")

    parser.add_argument('--ref_fn', type=str, default="ref.fa",
                        help="Reference fasta file input, default: %(default)s")

    parser.add_argument('--bed_fn', type=str, default=None,
                        help="Call variant only in these regions, works in intersection with ctgName, ctgStart and ctgEnd, optional, default: as defined by ctgName, ctgStart and ctgEnd")

    parser.add_argument('--output_fn', type=str, default=None,
                        help="Candidate sites VCF file input, if provided, will choose candidate +/- 1 or +/- 2. Use together with gen4Training. default: %(default)s")

    parser.add_argument('--ctgName', type=str, default="NC_000962.3",
                        help="The name of sequence to be processed, default: %(default)s")

    parser.add_argument('--ctgStart', type=int, default=None,
                        help="The 1-based starting position of the sequence to be processed")

    parser.add_argument('--ctgEnd', type=int, default=None,
                        help="The 1-based inclusive ending position of the sequence to be processed")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', default: %(default)s")

    parser.add_argument('--pos', type=str, default=None,
                        help="Tensor output, use PIPE for standard output, default: %(default)s")

    parser.add_argument('--platform', type=str, default='ont',
                        help="Sequencing platform of the input. Options: 'ont,illumina', default: %(default)s")

    parser.add_argument('--output_bam', action='store_true',
                        help="test in haploid mode in illumina data, works for TB")

    parser.add_argument('--ilmn_realign_folder', type=str, default=None,
                        help="The folder to store illumina realigned bam files, only works for illumina platform, default: %(default)s")


    args = parser.parse_args()

    local_realignment(args)

if __name__ == "__main__":
    main()
