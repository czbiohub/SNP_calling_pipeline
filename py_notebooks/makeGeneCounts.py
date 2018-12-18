#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
# script: makeGeneCounts.py
# author: Lincoln 
# date: 12.18.18
#
# TODO: ADD DESCRIPTION 
#
# usage: 
#			python3 makeGeneCounts.py
#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
import vcf
import numpy as np
import VCF # comes from Kamil Slowikowski
import os
import csv
import pandas as pd
import sys
import itertools

#////////////////////////////////////////////////////////////////////
# getFileNames()
#	Get file names based on the specified path
#
#////////////////////////////////////////////////////////////////////
def getFileNames():
	files = []
	for file in os.listdir("../vcf_test/"):
		if file.endswith(".vcf"):
			fullPath = (os.path.join("../vcf_test/", file))
			files.append(fullPath)
    
	return files

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
# getGenomePos()
#	Returns a genome position sting that will match against the ones w/in COSMIC db
#
#////////////////////////////////////////////////////////////////////
def getGenomePos(sample):
	chr = sample[0]
	chr = chr.replace("chr", "")
	pos = sample[1]
	ref = sample[3]
	alt = sample[4]
	
	if (len(ref) == 1) & (len(alt) == 1): # most basic case
		secondPos = pos
		genomePos = chr + ':' + str(pos) + '-' + str(secondPos)
	elif (len(ref) > 1) & (len(alt) == 1):
		secondPos = pos + len(ref)
		genomePos = chr + ':' + str(pos) + '-' + str(secondPos)
	elif (len(alt) > 1) & (len(ref) == 1):
		secondPos = pos + len(alt)
		genomePos = chr + ':' + str(pos) + '-' + str(secondPos)
	else: # BOTH > 1 .... not sure what to do here. does this actually happen? 
		secondPos = 'dummy'
		genomePos = chr + ':' + str(pos) + '-' + str(secondPos)

	return(genomePos)

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
		cell = f.replace("../vcf_test/", "")
		cell = cell.replace(".vcf", "")

		df = VCF.dataframe(f)
		genomePos_query = df.apply(getGenomePos, axis=1) # apply function for every row in df
    
		shared = list(set(genomePos_query) & set(genomePos_laud_db))
		cells_dict_laud.update({cell : len(shared)})

	print('finished!')
	return cells_dict_laud

#////////////////////////////////////////////////////////////////////
# writeCSV()
#	Writes the contents of a dictionary object to a csv
#
#////////////////////////////////////////////////////////////////////
def writeCSV(dictObj, outFile):
	print('writing csv')
	with open(outFile, 'w') as csv_file:
		writer = csv.writer(csv_file)
		for key, value in dictObj.items():
			writer.writerow([key, value])

#////////////////////////////////////////////////////////////////////
# main()
#	TODO: ADD DESCRIPTION
#////////////////////////////////////////////////////////////////////

global database
global database_laud

# filter counts LAUD
print('setting up COSMIC database...')
database = pd.read_csv("../CosmicGenomeScreensMutantExport.tsv", delimiter = '\t')
database_laud = getLAUD_db()
fNames = getFileNames()
filterDict1 = getFilterCountsLAUD(fNames) 
print("filter counts (LAUD) done!")
print('writing csv')
writeCSV(filterDict1, "foo.csv")

#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
