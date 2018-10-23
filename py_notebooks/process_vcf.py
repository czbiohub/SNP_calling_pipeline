#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
# script: process_vcf.py
# author: Lincoln 
# date: 10.11.18
#
# Want to turn my jupyter notebook into a python script
#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
import vcf
import numpy as np
import VCF # comes from Kamil Slowikowski
import os
import csv
import pandas as pd
import sys

#////////////////////////////////////////////////////////////////////
# getFileNames()
#	Get file names based on the specified path
#
#////////////////////////////////////////////////////////////////////
def getFileNames():
	files = []
	for file in os.listdir("../vcf/"):
		if file.endswith(".vcf"):
			fullPath = (os.path.join("../vcf/", file))
			files.append(fullPath)
    
	return files

#////////////////////////////////////////////////////////////////////
# getRawCounts()
#	Creates dictionary obj with raw counts for GATK hits w/in a given set of vcf files
#
#////////////////////////////////////////////////////////////////////
def getRawCounts(fileNames):
	print('getting raw counts...')
	cells_dict = {}

	for f in fileNames:
		cell = f.replace("../vcf/", "")
		cell = cell.replace(".vcf", "")
    
		df = VCF.dataframe(f)
		unique = len(np.unique(df.POS))
    
		cells_dict.update({cell : unique})
	print('finished!')
	return cells_dict

#////////////////////////////////////////////////////////////////////
# getGenomePos()
#	Returns a genome position sting that will match against the ones w/in COSMIC db
#
#////////////////////////////////////////////////////////////////////
def getGenomePos(sample):
	chr = sample[0]
	chr = chr.replace("chr", "")
	pos = sample[1]
	genomePos = chr + ':' + str(pos) + '-' + str(pos)

	return(genomePos)

#////////////////////////////////////////////////////////////////////
# getFilterCountsBasic()
#	Creates dictionry obj with COSMIC filtered GATK hits w/in a given set of vcfs 
#
#////////////////////////////////////////////////////////////////////
def getFilterCountsBasic(fileNames):
	print('getting filter counts basic...')
	cells_dict_filter = {}
	genomePos_db = pd.Series(database['Mutation genome position'])

	for f in fileNames:
		cell = f.replace("../vcf/", "")
		cell = cell.replace(".vcf", "")
		print(cell)
		df = VCF.dataframe(f)

		genomePos_query = df.apply(getGenomePos, axis=1)
    
		shared = list(set(genomePos_query) & set(genomePos_db))
		cells_dict_filter.update({cell : len(shared)})
    
		print(cells_dict_filter)
	print('finished!')
	return cells_dict_filter

#////////////////////////////////////////////////////////////////////
# getLAUD_db()
#	Return the cosmic database after lung adeno filter
#
#////////////////////////////////////////////////////////////////////
def getLAUD_db():
	print('setting up LAUD filtered database...')
	pHistList = database.index[database['Primary histology'] == 'carcinoma'].tolist()
	pSiteList = database.index[database['Primary site'] == 'lung'].tolist()
	shared = list(set(pHistList) & set(pSiteList))
	database_filter = database.iloc[shared]
	return database_filter

#////////////////////////////////////////////////////////////////////
# getFilterCountsLAUD()
#	Creates dictionry obj with COSMIC filtered GATK hits w/in a given set of vcfs 
#
#////////////////////////////////////////////////////////////////////
def getFilterCountsLAUD(fileNames):
	print('getting filter counts LAUD...')
	cells_dict_laud = {}
	genomePos_laud_db = pd.Series(database_laud['Mutation genome position'])

	for f in fileNames:
		cell = f.replace("../vcf/", "")
		cell = cell.replace(".vcf", "")
		#print(cell)
		df = VCF.dataframe(f)
    
		genomePos_query = df.apply(getGenomePos, axis=1) # apply function for every row in df
    
		shared = list(set(genomePos_query) & set(genomePos_laud_db))
		cells_dict_laud.update({cell : len(shared)})
        
		#print(cells_dict_laud)
	print('finished!')
	return cells_dict_laud

#////////////////////////////////////////////////////////////////////
# writeCSV()
#	Writes the contents of a dictionary object to a csv
#
#////////////////////////////////////////////////////////////////////
def writeCSV(dictObj, outFile):
	with open(outFile, 'w') as csv_file:
		writer = csv.writer(csv_file)
		for key, value in dictObj.items():
			writer.writerow([key, value])

#////////////////////////////////////////////////////////////////////
# main()
#	Main logic here. Comment out code blocks depending on which output file you want, nonImmune_GATK_hits_raw.csv, 
#	nonImmune_GATK_hits_COSMIC_filter.csv, or nonImmune_GATK_hits_COSMIC_filter_adv.csv
#////////////////////////////////////////////////////////////////////

global database
global database_laud

if len(sys.argv) != 2:
	print('usage: python3 process_vcf.py [-h] [1] [2] [3]')
	print('  ')
	print('		1 - getRawCounts')
	print('		2 - getFilterCountsBasic')
	print('		3 - getFilterCountsLAUD')
	print('  ')
	sys.exit()

# raw counts
if sys.argv[1] == '1':
	fNames = getFileNames()
	rawDict = getRawCounts(fNames)
	print("raw counts done!")
	print('writing csv')
	writeCSV(rawDict, "nonImmune_GATK_hits_raw.csv")

# filter counts (basic)
if sys.argv[1] == '2':
	print('setting up COSMIC database...')
	database = pd.read_csv("../CosmicGenomeScreensMutantExport.tsv", delimiter = '\t')
	filterDict = getFilterCountsBasic(fNames)
	print("filter counts (basic) done!")
	print('writing csv')
	writeCSV(filterDict, "nonImmune_GATK_hits_COSMIC_filter.csv")

# filter counts LAUD
if sys.argv[1] == '3':
	print('setting up COSMIC database...')
	database = pd.read_csv("../CosmicGenomeScreensMutantExport.tsv", delimiter = '\t')
	database_laud = getLAUD_db()
	filterDict1 = getFilterCountsLAUD(fNames) 
	print("filter counts (LAUD) done!")
	print('writing csv')
	writeCSV(filterDict1, "nonImmune_GATK_hits_LAUD_filter.csv")

#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
