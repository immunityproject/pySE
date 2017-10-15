import os
import sys
import csv


if len(sys.argv) != 2:
    print 'usage: %s {top level foldx output directory}'
    sys.exit(1)

rootDirName = sys.argv[1]
outputFilename = rootDirName + '-summary.csv'
out = csv.writer(open(os.path.join(rootDirName, outputFilename), 'wb'))
out.writerow(('epitope', 'site', 'remnant', 'energy', 'dwt'))


def listdir(dirName):
    return filter(lambda subdirName: subdirName != outputFilename and not subdirName.startswith('.'),
                  os.listdir(dirName))


wtEnergy = None

for proteinDirShortName in listdir(rootDirName):
    proteinDirName = os.path.join(rootDirName, proteinDirShortName)
    for epitopeName in listdir(proteinDirName):
        epitopeDirName = os.path.join(proteinDirName, epitopeName)
        for fileName in listdir(epitopeDirName):
            baseFileName, site, proteinName, proteinCode, variant = fileName.replace(
                '.txt', '').split('_')
            filePath = os.path.join(epitopeDirName, fileName)
            for line in open(filePath):
                remnantFileName, energy = line.split()[0:2]
                remnant = remnantFileName.split('_')[0]
                if remnant == 'WTref':
                    wtEnergy = float(energy)
                else:
                    out.writerow((epitopeName, site, remnant,
                                  energy, float(energy) - wtEnergy))
