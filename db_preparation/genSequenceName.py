import argparse
import os
import sys
import pandas
import gzip
from Bio import SeqIO

FLAGS = None


def summary(db_dir, sequenceName):
    genomes = ['archaea', "bacteria", "vertebrate_mammalian", "viral", "fungi", "protozoa"]
    sequenceNameWriter = open(sequenceName, 'w') if sequenceName is not None else sys.stdout

    for genome in genomes:
        assembly_summary = pandas.read_csv(os.path.join(db_dir, genome + '/assembly_summary.txt'), dtype=str, sep='\t',
                                       header=1)

        for i in range(assembly_summary.shape[0]):
            assemblyAccession = assembly_summary['# assembly_accession'][i]
            prefix = assembly_summary['ftp_path'][i][assembly_summary['ftp_path'][i].find(assemblyAccession):]
            path = '%s/%s/%s/%s_genomic.fna.gz'%(db_dir, genome, assemblyAccession, prefix)
            try:
                with gzip.open(path, 'rt') as f:
                    for record in SeqIO.parse(f, 'fasta'):
                        sequenceNameWriter.write("%s\t%s\n"%(record.id, ' '.join(record.description.split(' ')[1:])))
            except (EOFError,FileNotFoundError) as e:
                print(path, e,'\n', 'The RefSeq data were incompletely downloaded.')
        
    if sequenceName is not None:
        sequenceNameWriter.close()


def plasmidSummary(db_dir, num, sequenceName):
    sequenceNameWriter = open(sequenceName, 'w') if sequenceName is not None else sys.stdout
    for i in range(1, num+1):
        path='%s/plasmid/PLA_00000000%d.1/plasmid.%d.1.genomic.fna.gz'%(db_dir, i,i)
        with gzip.open(path, 'rt') as f:
            for record in SeqIO.parse(f, 'fasta'):
                sequenceNameWriter.write("%s\t%s\n"%(record.id, ' '.join(record.description.split(' ')[1:])))
    
    if sequenceName is not None:
        sequenceNameWriter.close()
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='gen sequence_name for db update')

    parser.add_argument('--sequenceName', default=None, required=False)
    parser.add_argument('--db_dir', default='./')
    parser.add_argument('--num', type=int, default=8)
    parser.add_argument('--function', type=int, help='1-abhvfp, 2-plasmid', required=True)

    FLAGS, UNPARSED = parser.parse_known_args()
    if FLAGS.function == 1:
        summary(sequenceName = FLAGS.sequenceName, db_dir=FLAGS.db_dir)
    elif FLAGS.function == 2:
        plasmidSummary(sequenceName = FLAGS.sequenceName, db_dir=FLAGS.db_dir, num = FLAGS.num)
