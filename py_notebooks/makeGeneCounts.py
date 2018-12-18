#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
# script: makeGeneCounts.py
# author: Lincoln 
# date: 12.18.18
#
# TODO: ADD DESCRIPTION 
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
	for file in os.listdir("../vcf/"):
		if file.endswith(".vcf"):
			fullPath = (os.path.join("../vcf/", file))
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


#print('writing csv')
#writeCSV(dummy, "foo.csv")

#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
