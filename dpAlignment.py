import argparse
#from os import makedirs
#from os.path import join, isdir
from Bio import SeqIO

#from openpyxl import load_workbook
#import pyard



def parseArgs():
    print('Parsing commandline arguments..')
    parser = argparse.ArgumentParser()

    parser.add_argument("-v", "--verbose", help="verbose operation", action="store_true")
    parser.add_argument("-f", "--fasta", required=True, help="FASTA of DP protein sequences", type=str)


    args = parser.parse_args()

    return args


def filterAlleles(inputFastaFileName=None, minimumLength = 258, newline='\r\n'):
    print('Filtering allele file ' + str(inputFastaFileName))

    alleles = {}

    for record in SeqIO.parse(inputFastaFileName, "fasta"):
        fullName=record.description

        fullNameTokens = fullName.split(' ')
        alleleName = fullNameTokens[1]
        seqLength = int(fullNameTokens[2])
        currentSequence = record.seq
        print('Allele ' + str(alleleName) + ' is of length ' + str(seqLength))
        #print('Sequence:' + str(currentSequence))

        if(seqLength >= minimumLength):
            nomenclatureTokens = alleleName.split(':')
            twoFieldName = str(nomenclatureTokens[0]) + ':' + str(nomenclatureTokens[1])
            print('2field:' + str(twoFieldName))

            if(twoFieldName in alleles.keys()):
                if(seqLength > len(alleles[twoFieldName])):
                    alleles[twoFieldName] = currentSequence
                else:
                    print('This allele is somehow shorter, i did not expect that to happen.')
            else:
                alleles[twoFieldName] = currentSequence

    outputFileName = inputFastaFileName.replace('.fasta','_filtered.fasta')
    with open(outputFileName, 'w') as outputFile:
        for alleleName in sorted(list(alleles.keys())):
            print('Writing allele' + str(alleleName))
            outputFile.write('>' + alleleName + newline)
            outputFile.write(str(alleles[alleleName]) + newline)


if __name__ == '__main__':
    args=parseArgs()

    verbose = args.verbose
    if verbose:
        print("verbose mode active")


    filterAlleles(inputFastaFileName=args.fasta)

    print('All done.')
