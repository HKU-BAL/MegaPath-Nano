import sys
import shlex
import signal
import gc
from subprocess import PIPE
from argparse import ArgumentParser
from collections import namedtuple

import shared.param as param
from shared.utils import subprocess_popen, IUPAC_base_to_num_dict as BASE2NUM

import random
is_pypy = '__pypy__' in sys.builtin_module_names

ReferenceResult = namedtuple('ReferenceResult', ['name', 'start', 'end', 'sequence', 'is_faidx_process_have_error'])

def phredscore2raw_score(qual):
    return ord(qual) - 33

def PypyGCCollect(signum, frame):
    gc.collect()
    signal.alarm(60)


BASES = set(list(BASE2NUM.keys()) + ["-"])

no_of_positions = 2 * param.flankingBaseNum + 1
matrix_row = param.matrixRow
matrix_num = param.matrixNum

matrix_depth = 100
max_iter = 100
#TODO
base_quality_cut_off = 7 # or int from 0-100

def generate_tensor(ctg_name, alignments, center, reference_sequence, reference_start_0_based, minimum_coverage, tensor_fp):

    # select center base higher than threhold
    if len(alignments) > matrix_depth:
        high_base_quality_alignments_idx_set = set()
        for read_idx, alignment in zip(range(len(alignments)), alignments):
            for reference_position, queryAdv, reference_base, query_base, STRAND, base_quality in alignment:
                if center != reference_position + 1: continue
                if reference_base == '-': continue # insert base
                base_quality = phredscore2raw_score(base_quality) if base_quality != "" else 100
                if base_quality >= base_quality_cut_off: high_base_quality_alignments_idx_set.add(read_idx)

        if len(high_base_quality_alignments_idx_set) > matrix_depth:
            new_alignments = []
            for read_idx in high_base_quality_alignments_idx_set:
                new_alignments.append(alignments[read_idx])
            alignments = new_alignments

    max_depth = len(alignments)

    itervals = max_depth // matrix_depth if not max_depth % matrix_depth else max_depth // matrix_depth + 1
    itervals = min(max_iter, itervals)
    iter_alignment = []
    random_sample = random.sample(range(max_depth), max_depth)
    for idx in range(itervals):
        if itervals == 1:
            iter_alignment.append(alignments)
            break
        if idx < itervals - 1:
            temp_aligment = []
            for index in random_sample[idx * matrix_depth:(idx + 1) * matrix_depth]:
                temp_aligment.append(alignments[index])
            iter_alignment.append(temp_aligment)
        else:
            temp_aligment = []
            for index in random_sample[-matrix_depth:]:
                temp_aligment.append(alignments[index])
            iter_alignment.append(temp_aligment)

    for alignments in iter_alignment:
        flanking_base_num = param.flankingBaseNum
        tensor = [[[0] * matrix_num for _ in range(matrix_row)] for _ in range(no_of_positions)]
        depth = [0] * no_of_positions

        for alignment in alignments:
            for reference_position, queryAdv, reference_base, query_base, STRAND,base_quality in alignment:
                # base_quality = phredscore2raw_score(base_quality) if base_quality != "" else 100
                # if base_quality_cut_off is not None:
                #     if base_quality < base_quality_cut_off:
                #         continue
                if str(reference_base) not in BASES or str(query_base) not in BASES:
                    continue
                position_index = reference_position - center + (flanking_base_num + 1)
                if not (0 <= position_index < no_of_positions):
                    continue

                strand_offset = 4 if STRAND else 0
                if query_base != "-" and reference_base != "-":
                    depth[position_index] = depth[position_index] + 1
                    tensor[position_index][BASE2NUM[reference_base] + strand_offset][0] += 1
                    tensor[position_index][BASE2NUM[query_base] + strand_offset][1] += 1
                    tensor[position_index][BASE2NUM[reference_base] + strand_offset][2] += 1
                    tensor[position_index][BASE2NUM[query_base] + strand_offset][3] += 1
                elif query_base != "-" and reference_base == "-":
                    position_index = min(position_index + queryAdv, no_of_positions - 1)
                    tensor[position_index][BASE2NUM[query_base] + strand_offset][1] += 1
                elif query_base == "-" and reference_base != "-":
                    tensor[position_index][BASE2NUM[reference_base] + strand_offset][2] += 1
                else:
                    print("Should not reach here: %s, %s" % (reference_base, query_base), file=sys.stderr)

        new_reference_position = center - reference_start_0_based
        if new_reference_position - (flanking_base_num+1) < 0 or depth[flanking_base_num] < minimum_coverage:
            return
        l =  "%s %d %s %s" % (
            ctg_name,
            center,
            reference_sequence[new_reference_position-(flanking_base_num+1):new_reference_position + flanking_base_num],
            " ".join((" ".join(" ".join("%d" % x for x in innerlist) for innerlist in outerlist)) for outerlist in tensor)
        )
        if l != None:
            tensor_fp.stdin.write(l)
            tensor_fp.stdin.write("\n")


def candidate_position_generator_from(
    candidate_file_path,
    ctg_start,
    ctg_end,
    is_consider_left_edge,
    flanking_base_num,
    begin_to_end
):
    is_read_file_from_standard_input = candidate_file_path == "PIPE"
    if is_read_file_from_standard_input:
        candidate_file_path_output = sys.stdin
    else:
        candidate_file_path_process = subprocess_popen(shlex.split("gzip -fdc %s" % (candidate_file_path)))
        candidate_file_path_output = candidate_file_path_process.stdout

    is_ctg_region_provided = ctg_start is not None and ctg_end is not None

    for row in candidate_file_path_output:
        row = row.split()
        position = int(row[1])  # 1-based position

        if is_ctg_region_provided and not (ctg_start <= position <= ctg_end):
            continue

        if is_consider_left_edge:
            # i is 0-based
            for i in range(position - (flanking_base_num + 1), position + (flanking_base_num + 1)):
                if i not in begin_to_end:
                    begin_to_end[i] = [(position + (flanking_base_num + 1), position)]
                else:
                    begin_to_end[i].append((position + (flanking_base_num + 1), position))
        else:
            begin_to_end[position - (flanking_base_num + 1)] = [(position + (flanking_base_num + 1), position)]

        yield position

    if not is_read_file_from_standard_input:
        candidate_file_path_output.close()
        candidate_file_path_process.wait()
    yield -1


class TensorStdout(object):
    def __init__(self, handle):
        self.stdin = handle

    def __del__(self):
        self.stdin.close()


def reference_result_from(
    ctg_name,
    ctg_start,
    ctg_end,
    samtools,
    reference_file_path,
    expand_reference_region
):
    region_str = ""
    reference_start, reference_end = None, None
    have_start_and_end_positions = ctg_start != None and ctg_end != None
    if have_start_and_end_positions:
        reference_start, reference_end = ctg_start - expand_reference_region, ctg_end + expand_reference_region
        reference_start = 1 if reference_start < 1 else reference_start
        region_str = "%s:%d-%d" % (ctg_name, reference_start, reference_end)
    else:
        region_str = ctg_name

    faidx_process = subprocess_popen(shlex.split("%s faidx %s %s" % (samtools, reference_file_path, region_str)),)
    if faidx_process is None:
        return None

    reference_name = None
    reference_sequences = []
    for row in faidx_process.stdout:
        if reference_name is None:
            reference_name = row.rstrip().lstrip(">") or ""
        else:
            reference_sequences.append(row.rstrip())
    reference_sequence = "".join(reference_sequences)

    faidx_process.stdout.close()
    faidx_process.wait()

    return ReferenceResult(
        name=reference_name,
        start=reference_start,
        end=reference_end,
        sequence=reference_sequence,
        is_faidx_process_have_error=faidx_process.returncode != 0,
    )


def samtools_view_process_from(
    ctg_name,
    ctg_start,
    ctg_end,
    samtools,
    bam_file_path
):
    have_start_and_end_position = ctg_start != None and ctg_end != None
    region_str = ("%s:%d-%d" % (ctg_name, ctg_start, ctg_end)) if have_start_and_end_position else ctg_name

    return subprocess_popen(
        shlex.split("%s view -F %d %s %s" % (samtools, param.SAMTOOLS_VIEW_FILTER_FLAG, bam_file_path, region_str))
    )


def OutputAlnTensor(args):
    available_slots = 10000000
    samtools = args.samtools
    tensor_file_path = args.tensor_fn
    bam_file_path = args.bam_fn
    reference_file_path = args.ref_fn
    candidate_file_path = args.can_fn
    dcov = args.dcov
    is_consider_left_edge = not args.stop_consider_left_edge
    min_coverage = args.minCoverage
    minimum_mapping_quality = args.minMQ
    ctg_name = args.ctgName
    ctg_start = args.ctgStart
    ctg_end = args.ctgEnd

    reference_result = reference_result_from(
        ctg_name=ctg_name,
        ctg_start=ctg_start,
        ctg_end=ctg_end,
        samtools=samtools,
        reference_file_path=reference_file_path,
        expand_reference_region=param.expandReferenceRegion,
    )

    reference_sequence = reference_result.sequence if reference_result is not None else ""
    is_faidx_process_have_error = reference_result is None or reference_result.is_faidx_process_have_error
    have_reference_sequence = reference_result is not None and len(reference_sequence) > 0

    if reference_result is None or is_faidx_process_have_error or not have_reference_sequence:
        print("Failed to load reference seqeunce. Please check if the provided reference fasta %s and the ctgName %s are correct." % (
            reference_file_path,
            ctg_name
        ), file=sys.stderr)
        sys.exit(1)

    reference_start = reference_result.start
    reference_start_0_based = 0 if reference_start is None else (reference_start - 1)
    begin_to_end = {}
    candidate_position = 0
    candidate_position_generator = candidate_position_generator_from(
        candidate_file_path=candidate_file_path,
        ctg_start=ctg_start,
        ctg_end=ctg_end,
        is_consider_left_edge=is_consider_left_edge,
        flanking_base_num=param.flankingBaseNum,
        begin_to_end=begin_to_end
    )

    samtools_view_process = samtools_view_process_from(
        ctg_name=ctg_name,
        ctg_start=ctg_start,
        ctg_end=ctg_end,
        samtools=samtools,
        bam_file_path=bam_file_path
    )

    center_to_alignment = {}

    if tensor_file_path != "PIPE":
        tensor_fpo = open(tensor_file_path, "wb")
        tensor_fp = subprocess_popen(shlex.split("gzip -c"), stdin=PIPE, stdout=tensor_fpo)
    else:
        tensor_fp = TensorStdout(sys.stdout)

    previous_position = 0
    depthCap = 0
    for l in samtools_view_process.stdout:
        l = l.split()
        if l[0][0] == "@":
            continue

        FLAG = int(l[1])
        POS = int(l[3]) - 1  # switch from 1-base to 0-base to match sequence index
        MQ = int(l[4])
        CIGAR = l[5]
        SEQ = l[9]
        reference_position = POS
        query_position = 0
        STRAND = (16 == (FLAG & 16))
        QUAL = l[10]

        if MQ < minimum_mapping_quality:
            continue

        end_to_center = {}
        active_set = set()

        while candidate_position != -1 and candidate_position < (POS + len(SEQ) + 100000):
            candidate_position = next(candidate_position_generator)

        if previous_position != POS:
            previous_position = POS
            depthCap = 0
        else:
            depthCap += 1
            if depthCap >= dcov:
                #print >> sys.stderr, "Bypassing POS %d at depth %d\n" % (POS, depthCap)
                continue

        advance = 0
        for c in str(CIGAR):
            if available_slots <= 0:
                break

            if c.isdigit():
                advance = advance * 10 + int(c)
                continue

            # soft clip
            if c == "S":
                query_position += advance

            # match / mismatch
            if c == "M" or c == "=" or c == "X":
                for _ in range(advance):
                    if reference_position in begin_to_end:
                        for rEnd, rCenter in begin_to_end[reference_position]:
                            if rCenter in active_set:
                                continue
                            end_to_center[rEnd] = rCenter
                            active_set.add(rCenter)
                            center_to_alignment.setdefault(rCenter, [])
                            center_to_alignment[rCenter].append([])
                    for center in list(active_set):
                        if available_slots <= 0:
                            break
                        available_slots -= 1

                        center_to_alignment[center][-1].append((
                            reference_position,
                            0,
                            reference_sequence[reference_position - reference_start_0_based],
                            SEQ[query_position],
                            STRAND,
                            QUAL[query_position]
                        ))
                    if reference_position in end_to_center:
                        center = end_to_center[reference_position]
                        active_set.remove(center)
                    reference_position += 1
                    query_position += 1

            # insertion
            if c == "I":
                for queryAdv in range(advance):
                    for center in list(active_set):
                        if available_slots <= 0:
                            break
                        available_slots -= 1

                        center_to_alignment[center][-1].append((
                            reference_position,
                            queryAdv,
                            "-",
                            SEQ[query_position],
                            STRAND,
                            QUAL[query_position]
                        ))
                    query_position += 1

            # deletion
            if c == "D":
                for _ in range(advance):
                    for center in list(active_set):
                        if available_slots <= 0:
                            break
                        available_slots -= 1

                        center_to_alignment[center][-1].append((
                            reference_position,
                            0,
                            reference_sequence[reference_position - reference_start_0_based],
                            "-",
                            STRAND,
                            ""
                        ))
                    if reference_position in begin_to_end:
                        for rEnd, rCenter in begin_to_end[reference_position]:
                            if rCenter in active_set:
                                continue
                            end_to_center[rEnd] = rCenter
                            active_set.add(rCenter)
                            center_to_alignment.setdefault(rCenter, [])
                            center_to_alignment[rCenter].append([])
                    if reference_position in end_to_center:
                        center = end_to_center[reference_position]
                        active_set.remove(center)
                    reference_position += 1

            # reset advance
            advance = 0

        if depthCap == 0:
            for center in list(center_to_alignment.keys()):
                if center + (param.flankingBaseNum + 1) >= POS:
                    continue
                generate_tensor(
                    ctg_name, center_to_alignment[center], center, reference_sequence, reference_start_0_based, min_coverage, tensor_fp
                )
                # if l != None:
                #     tensor_fp.stdin.write(l)
                #     tensor_fp.stdin.write("\n")
                available_slots += sum(len(i) for i in center_to_alignment[center])
                #print >> sys.stderr, "POS %d: remaining slots %d" % (center, available_slots)
                del center_to_alignment[center]

    for center in center_to_alignment.keys():
        generate_tensor(
            ctg_name, center_to_alignment[center], center, reference_sequence, reference_start_0_based, min_coverage, tensor_fp
        )
        # if l != None:
        #     tensor_fp.stdin.write(l)
        #     tensor_fp.stdin.write("\n")

    samtools_view_process.stdout.close()
    samtools_view_process.wait()
    if tensor_file_path != "PIPE":
        tensor_fp.stdin.close()
        tensor_fp.wait()
        tensor_fpo.close()


def main():
    parser = ArgumentParser(
        description="Generate tensors summarizing local alignments from a BAM file and a list of candidate locations")

    parser.add_argument('--bam_fn', type=str, default="input.bam",
                        help="Sorted bam file input, default: %(default)s")

    parser.add_argument('--ref_fn', type=str, default="ref.fa",
                        help="Reference fasta file input, default: %(default)s")

    parser.add_argument('--can_fn', type=str, default="PIPE",
                        help="Variant candidate list generated by ExtractVariantCandidates.py or true variant list generated by GetTruth.py, use PIPE for standard input, default: %(default)s")

    parser.add_argument('--tensor_fn', type=str, default="PIPE",
                        help="Tensor output, use PIPE for standard output, default: %(default)s")

    parser.add_argument('--minMQ', type=int, default=10,
                        help="Minimum Mapping Quality. Mapping quality lower than the setting will be filtered, default: %(default)d")

    parser.add_argument('--ctgName', type=str, default="chr17",
                        help="The name of sequence to be processed, default: %(default)s")

    parser.add_argument('--ctgStart', type=int, default=None,
                        help="The 1-based starting position of the sequence to be processed")

    parser.add_argument('--ctgEnd', type=int, default=None,
                        help="The 1-based inclusive ending position of the sequence to be processed")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', default: %(default)s")

    parser.add_argument('--stop_consider_left_edge', action='store_true',
                        help="If not set, would consider left edge only. That is, count the left-most base-pairs of a read for coverage even if the starting position of a read is after the starting position of a tensor")

    parser.add_argument('--dcov', type=int, default=1000,
                        help="Cap depth per position at %(default)d")

    parser.add_argument('--minCoverage', type=int, default=0,
                        help="Minimum coverage required to generate a tensor, default: %(default)d")

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    OutputAlnTensor(args)


if __name__ == "__main__":
    main()
