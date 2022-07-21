import argparse
from time import sleep
from urllib import request
import ast

def parseArgs():
    print('Parsing commandline arguments..')
    parser = argparse.ArgumentParser()

    parser.add_argument("-v", "--verbose", help="verbose operation", action="store_true")
    parser.add_argument("-i", "--inputfile", required=True, help="csv input file name. Row 1=Headers ; Separator=,", type=str)

    args = parser.parse_args()

    return args


def cleanup(rawAlleleText=None):
    cleanAlleleText=rawAlleleText.strip().replace('"','')
    return cleanAlleleText


def readInputFile(inputFileName=None, header=True, delimiter=','):
    print('Reading input file:' + str(inputFileName))
    patientTypes={}
    donorTypes={}

    with open(inputFileName,'r') as inputFile:
        for lineIndex, line in enumerate(inputFile.readlines()):
            #if (not header or lineIndex > 0) and lineIndex<100:
            if (not header or lineIndex > 0):

                if verbose:
                    print('line#' + str(lineIndex) + '=' + str(line.strip()))

                alleles = line.strip().split(delimiter)
                #if not len(alleles) == 4:
                #    raise Exception('There are not exactly 4 alleles on line# ' + str(lineIndex))

                patientTypes[lineIndex] = [cleanup(alleles[13]), cleanup(alleles[14])]
                donorTypes[lineIndex] = [cleanup(alleles[26]), cleanup(alleles[27])]

                if verbose:
                    print('Patient DP Types: ' + str(patientTypes[lineIndex]))
                    print('Donor DP Types: ' + str(donorTypes[lineIndex]))
    return patientTypes, donorTypes

def isInteger(queryText=None):
    try:
        integerText=int(queryText)
        return True
    except Exception as e:
        return False

def isUnambiguousTyping(queryTyping=None):
    #print('Checking typing:' + str(queryTyping))
    looksValid=True
    if(len(queryTyping.strip())==0):
        print('0-length string.')
        return True
    typingTokens=queryTyping.strip().split(':')
    for typingToken in typingTokens:
        if not isInteger(typingToken):
            looksValid=False
    '''
    if looksValid:
        print('Excellent:' + str(queryTyping))
    else:
        print('No good:' + str(queryTyping))
    '''
    return looksValid



def queryPairs(patientTypes=None, donorTypes=None, sleepSeconds=.01):
    print('Querying TCE compatibility for ' + str(len(patientTypes.keys())) + ' patient/donor pairs')

    patientTceGroups = {}
    donorTceGroups = {}
    matchGrades = {}

    for lineIndex in patientTypes.keys():
        try:
            # TODO: Do both typings exist?
            # TODO: Are either a Mac code?

            currentPatientTyping = patientTypes[lineIndex]
            currentDonorTyping = donorTypes[lineIndex]

            if(len(str(currentPatientTyping[1]).strip()))==0:
                currentPatientTyping[1]=currentPatientTyping[0]
            if(len(str(currentDonorTyping[1]).strip()))==0:
                currentDonorTyping[1]=currentDonorTyping[0]

            # Is it STILL empty? Then both types are missing:
            if((len(str(currentPatientTyping[1]).strip()))==0 or
                (len(str(currentDonorTyping[1]).strip()))==0):
                patientDpb1Allele = '?'
                patientDpb2Allele = '?'
                donorDpb1Allele = '?'
                donorDpb2Allele = '?'

                patientDpb1TceGroup = '?'
                patientDpb2TceGroup = '?'
                donorDpb1TceGroup = '?'
                donorDpb2TceGroup = '?'

                tcePredictedMatchGrade = '?'

            elif( not isUnambiguousTyping(currentPatientTyping[0])
                or not isUnambiguousTyping(currentPatientTyping[1])
                or not isUnambiguousTyping(currentDonorTyping[0])
                or not isUnambiguousTyping(currentDonorTyping[1])
                ):
                print('Skipping line')

                patientDpb1Allele = '?'
                patientDpb2Allele = '?'
                donorDpb1Allele   = '?'
                donorDpb2Allele   = '?'

                patientDpb1TceGroup = '?'
                patientDpb2TceGroup = '?'
                donorDpb1TceGroup   = '?'
                donorDpb2TceGroup   = '?'

                tcePredictedMatchGrade = '?'

            else:

                # Main URL
                requestUrl = 'https://www.ebi.ac.uk/cgi-bin/ipd/pl/hla/dpb_v2.cgi'

                # Patient info
                # ?pid=1P&patdpb1=01:01&patdpb2=02:01
                requestUrl = (requestUrl + '?pid=' + str(lineIndex) + 'P&patdpb1='
                    + str(patientTypes[lineIndex][0]) )

                requestUrl = (requestUrl + '&patdpb2='
                    + str(patientTypes[lineIndex][1]))

                # Donor Info
                # &did=2D&dondpb1=03:01&dondpb2=04:01
                requestUrl = (requestUrl + '&did=' + str(lineIndex) + 'D&dondpb1='
                    + str(donorTypes[lineIndex][0]) )

                requestUrl = (requestUrl + '&dondpb2='
                    + str(donorTypes[lineIndex][1]))

                if(verbose):
                    print('On line ' + str(lineIndex) + ' I found this URL:\n' + requestUrl)

                # Send Request
                tokenRequest = request.Request(url=requestUrl)
                requestResponse = request.urlopen(tokenRequest)
                responseData = ast.literal_eval(requestResponse.read().decode("UTF-8"))

                if(verbose):
                    print('Response:' + str(responseData))

                # Interpret JSON
                patientDpb1Allele = responseData['HLA-DPB1_TCE_report_V2.0']['patient']['gene']['HLA-DPB1'][0]['allele'].split('*')[1].strip()
                patientDpb2Allele = responseData['HLA-DPB1_TCE_report_V2.0']['patient']['gene']['HLA-DPB1'][1]['allele'].split('*')[1].strip()
                donorDpb1Allele   = responseData['HLA-DPB1_TCE_report_V2.0']['patient']['donor'][0]['gene']['HLA-DPB1'][0]['allele'].split('*')[1].strip()
                donorDpb2Allele   = responseData['HLA-DPB1_TCE_report_V2.0']['patient']['donor'][0]['gene']['HLA-DPB1'][1]['allele'].split('*')[1].strip()

                patientDpb1TceGroup = responseData['HLA-DPB1_TCE_report_V2.0']['patient']['gene']['HLA-DPB1'][0]['tce_group'].strip()
                patientDpb2TceGroup = responseData['HLA-DPB1_TCE_report_V2.0']['patient']['gene']['HLA-DPB1'][1]['tce_group'].strip()
                donorDpb1TceGroup   = responseData['HLA-DPB1_TCE_report_V2.0']['patient']['donor'][0]['gene']['HLA-DPB1'][0]['tce_group'].strip()
                donorDpb2TceGroup   = responseData['HLA-DPB1_TCE_report_V2.0']['patient']['donor'][0]['gene']['HLA-DPB1'][1]['tce_group'].strip()

                tcePredictedMatchGrade = responseData['HLA-DPB1_TCE_report_V2.0']['patient']['donor'][0]['results']['tce_prediction']

            '''
            # Make sure my data lines up, for error detection...
            if(patientDpb1Allele != str(patientTypes[lineIndex][0])
                or patientDpb2Allele != str(patientTypes[lineIndex][1])
                or donorDpb1Allele != str(donorTypes[lineIndex][0])
                or donorDpb2Allele != str(donorTypes[lineIndex][1])):
                # I hope this never happens.
                print('Patient Alleles Input:' + str(patientTypes[lineIndex]))
                print('Patient Alleles Output:' + str([patientDpb1Allele,patientDpb2Allele]))
                print('Donor Alleles Input:' + str(donorTypes[lineIndex]))
                print('Donor Alleles Output:' + str([donorDpb1Allele,donorDpb2Allele]))
                raise Exception('Apparently the data does not add up, check the inputs and outputs on line# ' + str(lineIndex))
            '''


            '''
            print('patientDpb1Allele=' + str(patientDpb1Allele))
            print('patientDpb2Allele=' + str(patientDpb2Allele))
            print('donorDpb1Allele=' + str(donorDpb1Allele))
            print('donorDpb2Allele=' + str(donorDpb2Allele))

            print('patientDpb1TceGroup=' + str(patientDpb1TceGroup))
            print('patientDpb2TceGroup=' + str(patientDpb2TceGroup))
            print('donorDpb1TceGroup=' + str(donorDpb1TceGroup))
            print('donorDpb2TceGroup=' + str(donorDpb2TceGroup))
            
            print('Predicted Match Grade: ' + str(tcePredictedMatchGrade))            
            '''

            # Store immunogenicities & grade.
            patientTceGroups[lineIndex] = [patientDpb1TceGroup,patientDpb2TceGroup]
            donorTceGroups[lineIndex] = [donorDpb1TceGroup,donorDpb2TceGroup]
            matchGrades[lineIndex] = tcePredictedMatchGrade

            if(verbose):
                print('Sleeping for ' + str(sleepSeconds) + ' seconds')
            sleep(sleepSeconds)
        except Exception as e:
            print('!!!Warning! there is a problem interpreting line # ' + str(lineIndex) + ':' + str(e))
            raise(e)

    return patientTceGroups, donorTceGroups, matchGrades


def writeOutputFile(inputFileName=None, patientTypes=None, donorTypes=None, patientTceGroups=None, donorTceGroups=None, matchGrades=None, delimiter = ',' ,newline='\r\n'):
    outputFileName = str(inputFileName) + '.Results.csv'
    print('Writing output file:' + str(outputFileName))

    with open(outputFileName, 'w') as outputFile:
        # header
        outputFile.write('P1' + delimiter + 'P2' + delimiter + 'D1' + delimiter + 'D2' + delimiter + 'tce_prediction' + newline)
        for lineIndex in patientTypes.keys():
            try:
                currentDataLine = (str(patientTypes[lineIndex][0]) + ' - Group ' + str(patientTceGroups[lineIndex][0])
                    + delimiter + str(patientTypes[lineIndex][1]) + ' - Group ' + str(patientTceGroups[lineIndex][1])
                    + delimiter + str(donorTypes[lineIndex][0]) + ' - Group ' + str(donorTceGroups[lineIndex][0])
                    + delimiter + str(donorTypes[lineIndex][1]) + ' - Group ' + str(donorTceGroups[lineIndex][1])
                    + delimiter + str(matchGrades[lineIndex]) + newline)
                outputFile.write(currentDataLine)
            except Exception as e:
                print('!!!Warning! there is a problem writing outputfile on line # ' + str(lineIndex) + ':' + str(e))
                raise(e)

if __name__ == '__main__':
    args=parseArgs()

    verbose = args.verbose
    if verbose:
        print("verbose mode active")

    patientTypes,donorTypes = readInputFile(inputFileName=args.inputfile)

    patientTceGroups, donorTceGroups, matchGrades = queryPairs(patientTypes=patientTypes, donorTypes=donorTypes)

    writeOutputFile(inputFileName=args.inputfile, patientTypes=patientTypes,donorTypes=donorTypes, patientTceGroups=patientTceGroups, donorTceGroups=donorTceGroups, matchGrades=matchGrades)

    print('All done.')
