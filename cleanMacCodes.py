import argparse
#from os import makedirs
#from os.path import join, isdir

from openpyxl import load_workbook
import pyard



def parseArgs():
    print('Parsing commandline arguments..')
    parser = argparse.ArgumentParser()

    parser.add_argument("-v", "--verbose", help="verbose operation", action="store_true")
    parser.add_argument("-x", "--excel", required=True, help="excel input file name. Row 1=Headers. Row2:X=HLA Typings", type=str)
    parser.add_argument("-c", "--columns", required=True, help="Columns that contain HLA typings with MAC Codes, separate by commas ex. G,H", type=str)

    args = parser.parse_args()

    return args

def isInteger(text=None):
    try:
        int(text)
        return True
    except ValueError:
        return False

def convertAlleleString(cellData=None, locus=None, ardObject=None):
    try:
        if(cellData is None or len(cellData) < 1):
            return ''
        #print('Converting cellData:' + str(cellData))

        nomenclatureTokens = cellData.split(':')
        alleleString = locus + '*' + cellData
        #print('allele string:' + str(alleleString))


        if (len(nomenclatureTokens)==2 and isInteger(nomenclatureTokens[0]) and nomenclatureTokens[1] == 'XX'):
            #print('This is an XX allele.')
            return cellData
        elif (len(nomenclatureTokens) == 2 and isInteger(nomenclatureTokens[0]) and not isInteger(nomenclatureTokens[1]) and nomenclatureTokens[1] != 'XX'):
            #print('This is a MAC Code.')
            expandMac = ardObject.expand_mac(alleleString)

            # Get rid of the locus again. Not sure if this is correct.
            for alleleIndex, allele in enumerate(expandMac):
                expandMac[alleleIndex] = allele.split('*')[1]

            return '|'.join(expandMac)
        else:
            # This should be a normal allele. Re-return it.
            return cellData

    except Exception:
        raise Exception('What should i do with this cellData:' + str(cellData))

def loadExcelFile(excelFileName=None):
    if (verbose):
        print('Loading Donor File:' + str(excelFileName))
    workBook = load_workbook(excelFileName)
    return workBook


def exportFile(outputFileName=None, cleanedWorkbook=None):
    #outputFileName=join(outputDirectory,'ResultWithMatchScore.xlsx')
    if (verbose):
        print('Exporting Donor File:' + str(outputFileName))
    cleanedWorkbook.save(outputFileName)


def getArdObject():
    if (verbose):
        print('Creating ARD object for latest IMGT/HLA Release..')
    ard = pyard.ARD()
    return ard


def parseLocusName(rawLocusName=None, delimiter='_'):
    locusTokens = str(rawLocusName).strip().upper().split(delimiter)
    locus = locusTokens[0] + '-' + locusTokens[1]
    return locus


def cleanMacCodes(excelFileName=None, columns=None, delimiter=',', headers=True):
    columnList = columns.split(delimiter)
    if (verbose):
        print('Cleaning Mac Codes in file ' + excelFileName + ' for ' + str(len(columnList)) + ' columns: ' + str(columnList))
    excelData = loadExcelFile(excelFileName=excelFileName)
    ardObject = getArdObject()

    firstSheet=excelData[excelData.sheetnames[0]]
    if(headers):
        startRowIndex=2
    else:
        startRowIndex=1
    #endRowIndex=18 # For testing.
    endRowIndex=None

    for columnName in columnList:
        rawLocusName=str(firstSheet[str(columnName) + str(1)].value)
        print('Column ' + str(columnName) + ' has locus name ' + str(rawLocusName))
        locusName = parseLocusName(rawLocusName=rawLocusName)
        #print('Locus:' + str(locusName))

        for rowIndexRaw, row in enumerate(firstSheet.iter_rows(min_row=startRowIndex, max_row=endRowIndex, values_only=True)):
            excelRowNumber = rowIndexRaw + startRowIndex
            #print('Raw Index ' + str(rowIndexRaw) + ', Excel Index ' + str(excelRowNumber) + ':' + str(row))

            columnName=columnName.strip()
            #print('Checking Cell ' + str(columnName) + str(excelRowNumber))
            cellData=firstSheet[str(columnName) + str(excelRowNumber)].value
            #print('Found Data:' + str(cellData))
            cleanedCellData = convertAlleleString(cellData=cellData, locus=locusName, ardObject=ardObject)
            firstSheet[str(columnName) + str(excelRowNumber)] =cleanedCellData

    return excelData







if __name__ == '__main__':
    args=parseArgs()

    verbose = args.verbose
    if verbose:
        print("verbose mode active")


    cleanedFileName = args.excel.replace('.xlsx','.cleanedMacCodes.xlsx')

    cleanedWorkbook = cleanMacCodes(excelFileName=args.excel, columns=args.columns)
    exportFile(outputFileName=cleanedFileName, cleanedWorkbook=cleanedWorkbook)

    print('All done.')
