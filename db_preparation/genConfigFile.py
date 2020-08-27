import sys
import os
import argparse

class Assembly:
    #assemblyAccession, refseqCategory, taxid, speciesTaxid, assemblyLevel, organismName
    def __init__(self, assemblyAccession, refseqCategory, taxid, speciesTaxid, assemblyLevel, organismName):
        self.assemblyAccession = assemblyAccession
        self.refseqCategory = refseqCategory
        self.taxid = taxid
        self.speciesTaxid = speciesTaxid
        self.assemblyLevel = assemblyLevel
        self.organismName = organismName
    
    def __str__(self):
        return "%s\t%s\t%s\t%s"%(self.assemblyAccession, self.speciesTaxid, self.taxid, self.organismName)

def main(assemblySummary, outputFile, outputFile2):
    # assemblySummary are already sorted by species taxid
    fileWriter = open(outputFile, 'w') if outputFile is not None else sys.stdout
    fileWriter2 = open(outputFile2, 'w') if outputFile2 is not None else sys.stdout
    prevSpeciesTaxid = None
    assemblyDict = {}
    referenceFound = False
    representativeFound = False
    assemblySummaryList = list(filter(None, assemblySummary.split(',')))
    for assemblySummary in assemblySummaryList:
        with open(assemblySummary, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                eles = line.strip().split('\t')
                speciesTaxid = eles[6]
                if prevSpeciesTaxid is None:
                    prevSpeciesTaxid = speciesTaxid
                if prevSpeciesTaxid != speciesTaxid:
                    for assemblyKey in assemblyDict:
                        assemblyInfo = assemblyDict[assemblyKey]
                        if assemblyInfo.refseqCategory == 'reference genome':
                            fileWriter.write('%s\n'%(assemblyInfo))
                        if assemblyInfo.refseqCategory == 'representative genome':
                            if referenceFound:
                                if assemblyInfo.assemblyLevel == 'Scaffold' or assemblyInfo.assemblyLevel == 'Contig': 
                                    continue
                                fileWriter2.write('%s\n'%(assemblyInfo))
                            else:
                                fileWriter.write('%s\n'%(assemblyInfo))
                        if assemblyInfo.refseqCategory == 'na':
                            if referenceFound or representativeFound:
                                if assemblyInfo.assemblyLevel == 'Scaffold' or assemblyInfo.assemblyLevel == 'Contig':
                                    continue
                                fileWriter2.write('%s\n'%(assemblyInfo))
                            else:
                                fileWriter.write('%s\n'%(assemblyInfo))                
                    assemblyDict.clear()
                    referenceFound = False
                    representativeFound = False
                    prevSpeciesTaxid = speciesTaxid
            
                assemblyInfo = Assembly(eles[0], eles[4], eles[5], eles[6], eles[11], eles[7])
                if eles[4] == 'reference genome':
                    referenceFound = True
                if eles[4] == 'representative genome':
                    representativeFound = True
                
                assemblyDict[eles[0]] = assemblyInfo

    for assemblyKey in assemblyDict:
        assemblyInfo = assemblyDict[assemblyKey]
        #sys.stderr.write('%s\n'%(assemblyInfo))
        if assemblyInfo.refseqCategory == 'reference genome':
            fileWriter.write('%s\n'%(assemblyInfo))
        if assemblyInfo.refseqCategory == 'representative genome':
            if referenceFound:
                if assemblyInfo.assemblyLevel == 'Scaffold' or assemblyInfo.assemblyLevel == 'Contig':
                    continue
                fileWriter2.write('%s\n'%(assemblyInfo))
            else:
                fileWriter.write('%s\n'%(assemblyInfo))
        if assemblyInfo.refseqCategory == 'na':
            if referenceFound or representativeFound:
                if assemblyInfo.assemblyLevel == 'Scaffold' or assemblyInfo.assemblyLevel == 'Contig':
                    continue
                fileWriter2.write('%s\n'%(assemblyInfo))
            else:
                fileWriter.write('%s\n'%(assemblyInfo))


    if outputFile != None:
        fileWriter.close()
    if outputFile2 != None:
        fileWriter2.close()


def humanGenomeSet(assemblySummary, outputFile):
    fileWriter = open(outputFile, 'w') if outputFile is not None else sys.stdout
    with open(assemblySummary, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            eles = line.strip().split('\t')
            assemblyInfo = Assembly(eles[0], eles[4], eles[5], eles[6], eles[11], eles[7])
            fileWriter.write("%s\n"%(assemblyInfo))
    if outputFile is not None:
        fileWriter.close()

def plasmidGenomeSet(num, outputFile):
    fileWriter = open(outputFile, 'w') if outputFile is not None else sys.stdout
    for i in range(1, num+1):
        fileWriter.write("PLA_00000000%d.1"%(i))
    if outputFile is not None:
        fileWriter.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='generate assembly_id.genome_set and assembly_id.genome_set')
    parser.add_argument('--assemblySummary', help='input assembly summary file', required=False, default=None)
    parser.add_argument('--outputFile', help='output file for species_id.genome_set', required=False)
    parser.add_argument('--outputFile2', help='output file for assembly_id.genome_set', required=False)
    parser.add_argument('--num', type=int)
    parser.add_argument('--function', help='1-main, 2-generate human genome set', type=int, default=None, required=True)
    FLAGS = parser.parse_args()
    if FLAGS.function == 1:
        main(assemblySummary = FLAGS.assemblySummary, outputFile = FLAGS.outputFile, outputFile2 = FLAGS.outputFile2)
    elif FLAGS.function == 2:
        humanGenomeSet(assemblySummary = FLAGS.assemblySummary, outputFile = FLAGS.outputFile)
    elif FLAGS.function == 3:
        plasmidGenomeSet(num = FLAGS.num, outputFile = FLAGS.outputFile)
