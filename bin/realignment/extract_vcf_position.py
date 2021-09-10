import sys
import shlex
from argparse import ArgumentParser
from sys import stderr
from subprocess import PIPE, Popen
import json
import logging
import re

def subprocess_popen(args, stdin=None, stdout=PIPE, stderr=stderr, bufsize=8388608):
    return Popen(args, stdin=stdin, stdout=stdout, stderr=stderr, bufsize=bufsize, universal_newlines=True)

def decode_alt(fn):
    alt_dict = {}
    f = open(fn, 'r').read().rstrip()
    lines = f.split('\n')
    for line in lines:
        line = line.rstrip().split('\t')
        try:
            contig, pos, depth = line[0], int(line[1]), int(line[2])
        except IndexError:
            return None
        try:
            alt_info = json.loads(line[3])
        except:
            seqs = line[3].split(' ')
            alt_info = dict(zip(seqs[::2], [int(item) for item in seqs[1::2]])) if len(seqs) else {}
        alt_dict[contig + ' ' + str(pos)] = [depth, alt_info]
    return alt_dict

def find_AF(depth, alt_info, ref_base, alt_base):
    count = 0
    if len(ref_base) == len(alt_base) and len(ref_base) == 1: #SNP
        if alt_base in alt_info:
            count = int(alt_info[alt_base])

    elif len(ref_base) < len(alt_base): # Insertion
        ins_base = 'I' + alt_base[1:]
        if ins_base in alt_info:
            count = int(alt_info[ins_base])
    elif len(ref_base) > len(alt_base): #Deletion
        del_base = 'D' + ref_base[1:]
        if del_base in alt_info:
            count = int(alt_info[del_base])
    if count > 0:
        af = (count) / float(depth)
        return af
    return None

def ExtractVcfPosition(args):
    vcf_fn = args.vcf_fn
    output_fn = args.output_fn
    unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (vcf_fn)))

    alt_dict = decode_alt(args.alt_fn)
    if alt_dict is None:
        return logging.error('[ERROR] Index error')
    header = []
    vcf_info = []
    for row in unzip_process.stdout:
        if type(row) == bytes:
            row = row.decode()
        row = row.rstrip()
        if row[0] == '#':
            header.append(row)
            continue
        columns = row.strip().split('\t')
        ctg_name = columns[0]
        pos = int(columns[1])
        key = ctg_name + ' ' + str(pos)
        if key in alt_dict:
            # update allele frequency and depth in clair vcf output
            ref_base = columns[3]
            alt_base = columns[4]
            try:
                GT, GQ, DP, AF = columns[-1].split(':')
                DP, alt_info = alt_dict[key]
                NEW_AF = find_AF(DP, alt_info, ref_base, alt_base)
                if NEW_AF and NEW_AF > 0:
                    print('[INFO] {}:{} updates allele frequency to {}'.format(ctg_name, pos, round(NEW_AF,4)))
                    last_column = "%s:%s:%d:%.4f" % (GT, GQ, DP, NEW_AF)
                    columns = '\t'.join(columns[:-1] + [last_column])
                    vcf_info.append(columns)
                else:
                    vcf_info.append(row)
            except:
                AF= float(re.findall("AF=\d+\.\d+", row)[0][3:])
                DP, alt_info = alt_dict[key]
                NEW_AF = find_AF(DP, alt_info, ref_base, alt_base)
                if NEW_AF and NEW_AF > 0:
                    print('[INFO] {}:{} updates allele frequency to {}'.format(ctg_name, pos, round(NEW_AF,4)))
                    last_column = "%d:%.4f" % (DP, NEW_AF)
                    columns = '\t'.join(columns[:] + [last_column])
                    vcf_info.append(columns)
                else:
                    vcf_info.append(row)

    output = header + vcf_info
    with open(output_fn, 'w') as output_file:
        output_file.write('\n'.join(output))
    unzip_process.stdout.close()
    unzip_process.wait()

def main():
    parser = ArgumentParser(description="Generate 1-based variant candidates using alignments")

    parser.add_argument(
        '--output_fn', type=str, default=None,
        help="Path to directory that stores small bins. (default: %(default)s)"
    )
    parser.add_argument(
        '--vcf_fn', type=str, default=None,
        help="Path of the output folder. (default: %(default)s)"
    )
    parser.add_argument(
        '--alt_fn', type=str, default=None,
        help="Path of the output folder. (default: %(default)s)"
    )
    parser.add_argument(
        '--realign_region_size', type=int, default=500,
        help="Name of the large bin. (default: %(default)s)"
    )
    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of sequence to be processed, default: %(default)s")

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    ExtractVcfPosition(args)


if __name__ == "__main__":
    main()
