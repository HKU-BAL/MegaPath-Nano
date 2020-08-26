import argparse
from os import listdir, path, system
import sys
import pandas
import gzip
from Bio import SeqIO

FLAGS = None

rankLookUpTable = {
    "no rank": 35,
    "domain": 34,
    "superkingdom": 33,
    "kingdom": 32,
    "subkingdom": 31,
    "superphylum": 30,
    "phylum": 29,
    "subphylum": 28,
    "superclass": 27,
    "class": 26,
    "subclass": 25,
    "infraclass": 24,
    "cohort": 23,
    "subcohort": 22,
    "superorder": 21,
    "order": 20,
    "parvorder": 19,
    "suborder": 18,
    "infraorder": 17,
    "superfamily": 16,
    "family": 15,
    "subfamily": 14,
    "tribe": 13,
    "subtribe": 12,
    "genus": 11,
    "subgenus": 10,
    "section": 9,
    "subsection": 8,
    "series": 7,
    "species group": 6,
    "species subgroup": 5,
    "species": 4,
    "subspecies": 3,
    "varietas": 2,
    "forma": 1
}

def parseNodesDmp(nodesDmp):
    taxid2parent = {}
    taxid2rank = {}
    with open(nodesDmp, 'r') as f:
        for line in f:
            eles = line.strip().split('|')
            taxid = eles[0].strip()
            parentTaxid = eles[1].strip()
            rank = eles[2].strip()
            taxid2parent[taxid] = parentTaxid
            taxid2rank[taxid] = rank

    return taxid2parent, taxid2rank


def summary(genome, nodesDmp, assemblyLength, assemblyPath, assemblyTaxid, sequenceSummary, db_dir):
    taxid2parent, taxid2rank = parseNodesDmp(nodesDmp)

    assembly_summary_path = os.path.join(db_dir, genome+'/assembly_summary.txt')
    assembly_summary = pandas.read_csv(assembly_summary_path, dtype=str, sep='\t',
                                       header=1)

    assemblyLengthWriter = open(assemblyLength, 'w') if assemblyLength is not None else sys.stdout
    assemblyPathWriter = open(assemblyPath, 'w') if assemblyPath is not None else sys.stdout
    assemblyTaxidWriter = open(assemblyTaxid, 'w') if assemblyTaxid is not None else sys.stdout
    sequenceSummaryWriter = open(sequenceSummary, 'w') if sequenceSummary is not None else sys.stdout

    for i in range(assembly_summary.shape[0]):
        assemblyAccession = assembly_summary['# assembly_accession'][i]
        taxid = assembly_summary['taxid'][i]
        speciesTaxid = assembly_summary['species_taxid'][i]
        prefix = assembly_summary['ftp_path'][i][assembly_summary['ftp_path'][i].find(assemblyAccession):]
        path = os.path.join(db_dir, "%s/%s/%s_genomic.fna.gz"%(genome, assemblyAccession, prefix))
        assemblyReportPath = os.path.join(db_dir, '%s/%s/%s_assembly_report.txt'%(genome, assemblyAccession, prefix))
        totalLength = 0
        with gzip.open(path, 'rt') as f:
            for record in SeqIO.parse(f, 'fasta'):
                totalLength += len(record)
                sequenceSummaryWriter.write("%s\t%d\t%s\n"%(record.id, len(record), assemblyAccession))
        
        assemblyLengthWriter.write("%s\t%d\n"%(assemblyAccession, totalLength))
        assemblyPathWriter.write("%s\t%s\n"%(assemblyAccession, path))
        if speciesTaxid in taxid2parent:
            if rankLookUpTable[taxid2rank[speciesTaxid]] > rankLookUpTable[taxid2rank[taxid2parent[speciesTaxid]]]:
                sys.stderr.write("Error rank: %s\t%s\t%s\t%s\n"%(speciesTaxid, taxid2parent[speciesTaxid], taxid2rank[speciesTaxid], taxid2rank[taxid2parent[speciesTaxid]]))
            assemblyTaxidWriter.write("%s\t%s\t%s\t%s\t%s\n"%(assemblyAccession, taxid, speciesTaxid, taxid2parent[speciesTaxid], rankLookUpTable[taxid2rank[taxid2parent[speciesTaxid]]]))
        else:
            sys.stderr.write("%s\t%s\t%s\n"%(assemblyAccession, taxid, speciesTaxid))

    if assemblyLength is not None:
        assemblyLengthWriter.close()
    if assemblyPath is not None:
        assemblyPathWriter.close()
    if assemblyTaxid is not None:
        assemblyTaxidWriter.close()
    if sequenceSummary is not None:
        sequenceSummaryWriter.close()

def summaryPlasmid(assemblyLength, assemblyPath, assemblyTaxid, sequenceSummary, num, db_dir):
    assemblyLengthWriter = open(assemblyLength, 'w') if assemblyLength is not None else sys.stdout
    assemblyPathWriter = open(assemblyPath, 'w') if assemblyPath is not None else sys.stdout
    assemblyTaxidWriter = open(assemblyTaxid, 'w') if assemblyTaxid is not None else sys.stdout
    sequenceSummaryWriter = open(sequenceSummary, 'w') if sequenceSummary is not None else sys.stdout
    for i in range(1, num+1):
        assemblyAccession = 'PLA_00000000%d.1'%(i)
        path = os.path.join(db_dir, 'plasmid/PLA_00000000%d.1/plasmid.%d.1.genomic.fna.gz'%(i, i))

        totalLength = 0
        with gzip.open(path, 'rt') as f:
            for record in SeqIO.parse(f, 'fasta'):
                totalLength += len(record)
                sequenceSummaryWriter.write("%s\t%d\t%s\n"%(record.id, len(record), assemblyAccession))

        assemblyLengthWriter.write("%s\t%d\n"%(assemblyAccession, totalLength))
        assemblyPathWriter.write("%s\t%s\n"%(assemblyAccession, path))
        assemblyTaxidWriter.write("%s\t100000000%d\t100000000%d\t1000000000\t%s\n"%(assemblyAccession, i, i, rankLookUpTable["no rank"]))

    if assemblyLength is not None:
        assemblyLengthWriter.close()
    if assemblyPath is not None:
        assemblyPathWriter.close()
    if assemblyTaxid is not None:
        assemblyTaxidWriter.close()
    if sequenceSummary is not None:
        sequenceSummaryWriter.close()

def main():
    if FLAGS.viral:
        print("Update viral")
        summary("viral", FLAGS.nodesDmp, FLAGS.assemblyLength, FLAGS.assemblyPath, FLAGS.assemblyTaxid, FLAGS.sequenceSummary, FLAGS.db_dir)
    if FLAGS.bacteria:
        print("Update bacteria")
        summary("bacteria", FLAGS.nodesDmp, FLAGS.assemblyLength, FLAGS.assemblyPath, FLAGS.assemblyTaxid, FLAGS.sequenceSummary, FLAGS.db_dir)
    if FLAGS.fungi:
        print("Update fungi")
        summary("fungi", FLAGS.nodesDmp, FLAGS.assemblyLength, FLAGS.assemblyPath, FLAGS.assemblyTaxid, FLAGS.sequenceSummary, FLAGS.db_dir)
    if FLAGS.protozoa:
        print("Update protozoa")
        summary("protozoa", FLAGS.nodesDmp, FLAGS.assemblyLength, FLAGS.assemblyPath, FLAGS.assemblyTaxid, FLAGS.sequenceSummary, FLAGS.db_dir)
    if FLAGS.archaea:
        print("Update archaea")
        summary("archaea", FLAGS.nodesDmp, FLAGS.assemblyLength, FLAGS.assemblyPath, FLAGS.assemblyTaxid, FLAGS.sequenceSummary, FLAGS.db_dir)
    if FLAGS.plasmid:
        print("Update plasmid")
        summaryPlasmid(FLAGS.assemblyLength, FLAGS.assemblyPath, FLAGS.assemblyTaxid, FLAGS.sequenceSummary, FLAGS.num, FLAGS.db_dir)
    if FLAGS.vertebrate_mammalian:
        print("Update vertebrate_mammalian")
        summary("vertebrate_mammalian", FLAGS.nodesDmp, FLAGS.assemblyLength, FLAGS.assemblyPath, FLAGS.assemblyTaxid, FLAGS.sequenceSummary, FLAGS.db_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Refseq assembly info update')

    parser.add_argument('--viral', dest='viral', default=False, action='store_true')
    parser.add_argument('--bacteria', dest='bacteria', default=False, action='store_true')
    parser.add_argument('--fungi', dest='fungi', default=False, action='store_true')
    parser.add_argument('--protozoa', dest='protozoa', default=False, action='store_true')
    parser.add_argument('--archaea', dest='archaea', default=False, action='store_true')
    parser.add_argument('--plasmid', dest='plasmid', default=False, action='store_true')
    parser.add_argument('--vertebrate_mammalian', dest='vertebrate_mammalian', default=False, action='store_true')
    parser.add_argument('--db_dir', default='./')

    parser.add_argument('--nodesDmp', default=None, required=False)
    parser.add_argument('--assemblyLength', default=None, required=False)
    parser.add_argument('--assemblyPath', default=None, required=False)
    parser.add_argument('--assemblyTaxid', default=None, required=False)
    parser.add_argument('--sequenceSummary', default=None, required=False)
    parser.add_argument('--num', type=int, default=8, required=False)

    FLAGS, UNPARSED = parser.parse_known_args()

    main()
