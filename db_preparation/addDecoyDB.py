import argparse
import os
from Bio import SeqIO

FLAGS = None

def main():
    decoy_name=os.path.basename(os.path.splitext(FLAGS.decoy_fasta))
    with open('%s/plasmid.genome_set' %(FLAGS.config_folder),'a') as config_plasmid:
        config_plasmid.write(decoy_name+'\n')
    assemblyLengthWriter = open('%s/assembly_length' %(FLAGS.db_dir), 'a')
    assemblyPathWriter = open('%s/assembly_path' %(FLAGS.db_dir), 'a')
    assemblyTaxidWriter = open('%s/assembly_tax_id' %(FLAGS.db_dir), 'a')
    sequenceSummaryWriter = open('%s/sequence_summary' %(FLAGS.db_dir), 'a')
    with open(FLAGS.decoy_fasta, 'rt') as fi:
        for record in SeqIO.parse(fi, 'fasta'):
            totalLength += len(record)
            sequenceSummaryWriter.write("%s\t%d\t%s\n"%(record.id, len(record), decoy_name))
    assemblyLengthWriter.write("%s\t%d\n"%(decoy_name, totalLength))
    assemblyPathWriter.write("%s\t%s\n"%(decoy_name, path))
    assemblyTaxidWriter.write("%s\t1000000099\t1000000099\t1000000001\t%s\n"%(decoy_name, '35'))  #arbitrary taxid


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add fasta to decoy db')
    parser.add_argument('--decoy_fasta', required=True)

    cwd=os.path.dirname(os.path.realpath(__file__))
    nano_dir=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    parser.add_argument('--db_dir', default='%s/genomes' %(nano_dir))
    parser.add_argument('--config_folder', help='Config file folder', default='%s/config' %(nano_dir))

    FLAGS, UNPARSED = parser.parse_known_args()

    main()
