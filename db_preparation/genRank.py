import argparse
import sys
import os


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

def main(ranks):
    fileWriter = open(ranks, 'w')
    for key in rankLookUpTable:
        fileWriter.write("%s\t%s\n"%(key, rankLookUpTable[key]))
    fileWriter.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='rank table')
    parser.add_argument('--ranks', default=None, required=False)
    FLAGS, UNPARSED = parser.parse_known_args()
    main(ranks = FLAGS.ranks)
