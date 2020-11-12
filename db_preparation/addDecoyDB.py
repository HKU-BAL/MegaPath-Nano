import argparse
from shutil import copyfile
import os
from Bio import SeqIO

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add fasta to decoy db')
    parser.add_argument('--decoy_fasta', required=True)

    cwd=os.path.dirname(os.path.realpath(__file__))
    nano_dir=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    parser.add_argument('--assembly_dir', default='%s/genomes' %(nano_dir))
    parser.add_argument('--config_folder', help='Config file folder', default='%s/config' %(nano_dir))

    FLAGS = parser.parse_args()

    decoy_name=os.path.basename(os.path.splitext(FLAGS.decoy_fasta))
    path= '%s/refseq/plasmid/%s' %(assembly_dir,decoy_name)
    copyfile(FLAGS.decoy_fasta,path)
    with open('%s/plasmid.genome_set' %(FLAGS.config_folder),'a') as config_plasmid:
        config_plasmid.write(decoy_name+'\n')
    assemblyLengthWriter = open('%s/assembly_length' %(FLAGS.assembly_dir), 'a')
    assemblyPathWriter = open('%s/assembly_path' %(FLAGS.assembly_dir), 'a')
    assemblyTaxidWriter = open('%s/assembly_tax_id' %(FLAGS.assembly_dir), 'a')
    sequenceSummaryWriter = open('%s/sequence_summary' %(FLAGS.assembly_dir), 'a')
    with open(path, 'rt') as fi:
        for record in SeqIO.parse(fi, 'fasta'):
            totalLength += len(record)
            sequenceSummaryWriter.write("%s\t%d\t%s\n"%(record.id, len(record), decoy_name))
    assemblyLengthWriter.write("%s\t%d\n"%(decoy_name, totalLength))
    assemblyPathWriter.write("%s\t%s\n"%(decoy_name, path))
    assemblyTaxidWriter.write("%s\t1000000099\t1000000099\t1000000001\t%s\n"%(decoy_name, '35'))  #arbitrary taxid

