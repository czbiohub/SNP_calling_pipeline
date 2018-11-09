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

		df = VCF.dataframe(f)
		genomePos_query = df.apply(getGenomePos, axis=1) # apply function for every row in df
    
		shared = list(set(genomePos_query) & set(genomePos_laud_db))
		cells_dict_laud.update({cell : len(shared)})

	print('finished!')
	return cells_dict_laud

#////////////////////////////////////////////////////////////////////
# hitSearchFunc()
#	Performs the actual search
#
#////////////////////////////////////////////////////////////////////
def hitSearchFunc(sample):
	match = 0
	currChrom = sample.split(':')[0]
	if currChrom == queryChrom:
		sub0 = sample.split('-')[0] # split on -
		sub1 = sample.split('-')[1] # this guy is good
		sub00 = sub0.split(':')[1] # split on :, need to get rid of chrom

		try:
			lPosCurr = sub00
			rPosCurr = sub1

			if (lPosCurr >= lPosQuery) & (lPosCurr <= rPosQuery): # left position good
				if (rPosCurr >= lPosQuery) & (rPosCurr <= rPosQuery): # right position good
					match = 1
		except IndexError:
			print('index error')

	return match

#////////////////////////////////////////////////////////////////////
# hitSearchFunc_coords()
#	Performs the actual search, and returns coords
#
#////////////////////////////////////////////////////////////////////
def hitSearchFunc_coords(sample):
	match = ""
	currChrom = sample.split(':')[0]
	if currChrom == queryChrom:
		sub0 = sample.split('-')[0] # split on -
		sub1 = sample.split('-')[1] # this guy is good
		sub00 = sub0.split(':')[1] # split on :, need to get rid of chrom

		try:
			lPosCurr = sub00
			rPosCurr = sub1

			if (lPosCurr >= lPosQuery) & (lPosCurr <= rPosQuery): # left position good
				if (rPosCurr >= lPosQuery) & (rPosCurr <= rPosQuery): # right position good
					match = lPosCurr
					#print(lPosCurr) # print out the actual SNP genome coord
		except IndexError:
			print('index error')

	return match

#////////////////////////////////////////////////////////////////////
# getGOIHits()
#	Creates dictionry obj with hits to a specific Gene of Interest
#
#////////////////////////////////////////////////////////////////////
def getGOIHits(fileNames, chrom, pos1, pos2):
	print('getting hits to GOI')

	global queryChrom, lPosQuery, rPosQuery # dont like this
	genomePos_laud_db = pd.Series(database_laud['Mutation genome position'])

	cells_dict_GOI = {}
	queryChrom = chrom
	lPosQuery = pos1
	rPosQuery = pos2

	for f in fileNames:
		numMatches = 0
		cell = f.replace("../vcf/", "")
		cell = cell.replace(".vcf", "")	

		df = VCF.dataframe(f)
		genomePos_query = df.apply(getGenomePos, axis=1) # apply function for every row in df
		
		shared = list(set(genomePos_query) & set(genomePos_laud_db)) # get the LAUD filter set
		shared1 = pd.Series(shared) # what if i convert this guy to a pandas object? 

		numMatches = shared1.apply(hitSearchFunc) # another apply call 

		cells_dict_GOI.update({cell : sum(numMatches)})
	
	return cells_dict_GOI

#////////////////////////////////////////////////////////////////////
# getGOIHit_coords()
#	Creates dictionry obj with genome coords for hits to specific GOI
#
#////////////////////////////////////////////////////////////////////
def getGOIHit_coords(fileNames, chrom, pos1, pos2):
	print('getting coords to GOI hits')

	global queryChrom, lPosQuery, rPosQuery # dont like this
	genomePos_laud_db = pd.Series(database_laud['Mutation genome position'])

	cells_dict_GOI_coords = {}
	queryChrom = chrom
	lPosQuery = pos1
	rPosQuery = pos2

	for f in fileNames:
		numMatches = 0
		cell = f.replace("../vcf/", "")
		cell = cell.replace(".vcf", "")	

		df = VCF.dataframe(f)
		genomePos_query = df.apply(getGenomePos, axis=1) # apply function for every row in df
		shared = list(set(genomePos_query) & set(genomePos_laud_db)) # get the LAUD filter set

		shared1 = pd.Series(shared) # what if i convert this guy to a pandas object? 
		matches = shared1.apply(hitSearchFunc_coords) # another apply call 

		# delete empty dict keys
		for k in matches.keys():
			try:
				if len(matches[k])<1:
					del matches[k]
			except: pass

		cells_dict_GOI_coords.update({cell : list(matches.values)})
		#print(list(matches.values))
	return cells_dict_GOI_coords

#////////////////////////////////////////////////////////////////////
# getMutationAA()
#	Pass in a dict of {cell, list(genomePos)} items and it returns a
#	dict of {cell, list(Mutation.AA)}
# 
#////////////////////////////////////////////////////////////////////
def getMutationAA(d, chr):
	print('AA searching')
	newDict = {}

	for k in d:
		valuesList = d.get(k) # can now handle values with multiple entries
		newValues = []

		for entry in valuesList:	
			chrStr = chr + ':' + entry + '-' + entry
			filter = database_laud["Mutation genome position"]==chrStr
			sub = database_laud.where(filter).dropna(axis=0, how='all')
			currMut = sub['Mutation AA']

			for item in currMut:		# really shouldnt have a for loop here
				item = item.replace("p.", "")

			newValues.append(item)
		
		newDict.update({k : newValues})

	return newDict

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
#	Main logic here. Prints usage, parses cmd line args, calls 
#		appropriate driver function 
#		
#		EGFR test: 
#			python3 process_vcf.py 4 7 55019101 55211628 egfr_out.csv
#////////////////////////////////////////////////////////////////////

global database
global database_laud

if len(sys.argv) == 1:
	print('usage: python3 process_vcf.py [-h] [1] [2] [3] [4]')
	print('  ')
	print('		1 - getRawCounts')
	print('		2 - getFilterCountsBasic')
	print('		3 - getFilterCountsLAUD')
	print('		4 - getGeneOfInterest')
	print('			needs [chrom] [pos1] [pos2] [outFile]')
	print('			ie. python3 process_vcf.py 4 7 500050 50010 myOut.csv')
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
	fNames = getFileNames()
	filterDict = getFilterCountsBasic(fNames)
	print("filter counts (basic) done!")
	print('writing csv')
	writeCSV(filterDict, "nonImmune_GATK_hits_COSMIC_filter.csv")

# filter counts LAUD
if sys.argv[1] == '3':
	print('setting up COSMIC database...')
	database = pd.read_csv("../CosmicGenomeScreensMutantExport.tsv", delimiter = '\t')
	database_laud = getLAUD_db()
	fNames = getFileNames()
	filterDict1 = getFilterCountsLAUD(fNames) 
	print("filter counts (LAUD) done!")
	print('writing csv')
	writeCSV(filterDict1, "nonImmune_GATK_hits_LAUD_filter.csv")

# Gene of interest
if sys.argv[1] == '4':

	if len(sys.argv) != 6:
		print('  ')
		print('USER ERROR')
		print('This function requires [chrom] [pos1] [pos2] [outFile]')
		print('	ie. python3 process_vcf.py 4 7 500000 50010 myFile.csv')
		print(' ')
		sys.exit()
	
	print('setting up COSMIC database...')
	database = pd.read_csv("../CosmicGenomeScreensMutantExport.tsv", delimiter = '\t')
	database_laud = getLAUD_db()
	fNames = getFileNames()
	
	chromo = sys.argv[2]
	position1 = sys.argv[3]
	position2 = sys.argv[4]	

	goiDict = getGOIHits(fNames, chromo, position1, position2) # standard call - get raw counts
	#goiDict = getGOIHit_coords(fNames, chromo, position1, position2) # get genome coords
	print("GOI search done!")
	
	outFilePref = sys.argv[5]
	writeCSV(goiDict, './out/' + outFilePref + '.csv')

	#goiDict_AA = getMutationAA(goiDict, chromo)
	#print('AA search done')
	#writeCSV(goiDict_AA, './out/' + outFilePref + '_AA.csv')

#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
