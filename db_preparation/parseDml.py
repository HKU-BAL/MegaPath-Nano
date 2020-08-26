import sys
import os
import argparse

def nodesMain(nodesDmp, outputFile):
    fileWriter = open(outputFile, 'w') if outputFile is not None else sys.stdout
    with open(nodesDmp, 'r') as f:
        for line in f:
            eles = line.split('|')
            for i in range(len(eles)):
                eles[i] = eles[i].strip()
                eles[i] = eles[i].replace('"', '')
            fileWriter.write("%s\t1\t1\n"%('\t'.join(eles[:5])))
    
    if outputFile is not None:
        fileWriter.close()

def namesMain(namesDmp, outputFile):
    fileWriter = open(outputFile, 'w') if outputFile is not None else sys.stdout
    prevNames = ''
    with open(namesDmp, 'r') as f:
        for line in f:
            eles = line.split('|')
            for i in range(len(eles)):
                eles[i] = eles[i].strip()
                eles[i] = eles[i].replace('"', '')
            if eles[0] == prevNames:
                fileWriter.write("%s\t1\t0\t1\n"%('\t'.join(eles[:4])))
            else:
                fileWriter.write("%s\t1\t1\t1\n"%('\t'.join(eles[:4])))
                prevNames = eles[0]

    if outputFile is not None:
        fileWriter.close()

def plasmidNames(outputFile, num):
    fileWriter = open(outputFile, 'w') if outputFile is not None else sys.stdout
    fileWriter.write("1000000000\tPLA_000000000\tPLA_000000000.1\t\t1\t1\t1\n")
    for i in range(1, num+1):
        fileWriter.write("100000000%d\tPLA_00000000%d.1\tPLA_00000000%d.1\t\t1\t1\t1\n"%(i, i, i))
    if outputFile is not None:
        fileWriter.close()

def plasmidNodes(outputFile, num):
    fileWriter = open(outputFile, 'w') if outputFile is not None else sys.stdout
    fileWriter.write("1000000000\t0\t35\t\t\t1\t1\n")
    for i in range(1, num+1):
        fileWriter.write("100000000%d\t1000000000\t35\t\t\t1\t1\n"%(i))
    if outputFile is not None:
        fileWriter.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='parse dmp')
    parser.add_argument('--dmp', default=None, required=False)
    parser.add_argument('--outputFile', default=None, required=False)
    parser.add_argument('--num', type=int)
    parser.add_argument('--function', type=int, help='1-names dmp, 2-nodes dmp, 3-plasmid names, 4-plasmid nodes', required=True)
    FLAGS, UNPARSED = parser.parse_known_args()
    if FLAGS.function == 1:
        namesMain(namesDmp = FLAGS.dmp, outputFile = FLAGS.outputFile)
    elif FLAGS.function == 2:
        nodesMain(nodesDmp = FLAGS.dmp, outputFile = FLAGS.outputFile)
    elif FLAGS.function == 3:
        plasmidNames(outputFile = FLAGS.outputFile, num = FLAGS.num)
    elif FLAGS.function == 4:
        plasmidNodes(outputFile = FLAGS.outputFile, num = FLAGS.num)
