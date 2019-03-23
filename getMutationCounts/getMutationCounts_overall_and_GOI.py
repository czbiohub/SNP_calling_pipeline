#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
# script: getMutationCounts_overall_and_GOI.py
# author: Lincoln 
# date: 10.11.18
#
# This script performs some basic analysis on vcf files, as output by
# my SNP_detection_pipeline. It has 4 separate run modes:
#	1. get raw mutation counts, for every cell
#	2. get mutation counts, after filtering through COSMIC database
#	3. get mutation counts, '                                      ', 
#			with the specific LAUD annotation
#	4. for a given GOI, which cells have mutations, and what are those
#		mutations, on the amino acid level? This creates the necessary
#		input for all of the lolliplot stuff. As well as for 
#			makeSummaryTable.ipynb 
# 
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
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

#////////////////////////////////////////////////////////////////////
# getFileNames()
#	Get file names based on the specified path
#
#////////////////////////////////////////////////////////////////////
def getFileNames():
	files = []
	for file in os.listdir("vcf_old_test/"):
		if file.endswith(".vcf"):
			fullPath = (os.path.join("vcf_old_test/", file))
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
		cell = f.replace("vcf_old_test/", "")
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
# getFilterCountsBasic()
#	Creates dictionry obj with COSMIC filtered GATK hits w/in a given set of vcfs 
#
#////////////////////////////////////////////////////////////////////
def getFilterCountsBasic(fileNames):
	print('getting filter counts basic...')
	cells_dict_filter = {}
	genomePos_db = pd.Series(database['Mutation genome position'])

	for f in fileNames:
		cell = f.replace("vcf_old_test/", "")
		cell = cell.replace(".vcf", "")
		print(cell)
		df = VCF.dataframe(f)

		genomePos_query = df.apply(getGenomePos, axis=1)
    
		shared = list(set(genomePos_query) & set(genomePos_db))
		cells_dict_filter.update({cell : len(shared)})
    
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
		cell = f.replace("vcf_old_test/", "")
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
#    	REMEMBER: `match` is just a boolean
#////////////////////////////////////////////////////////////////////
def hitSearchFunc(sample):
	match = 0
	currChrom = sample.split(':')[0]
	if currChrom == queryChrom:
		sub0 = sample.split('-')[0] # split on `-`
		sub1 = sample.split('-')[1] # this guy is good
		sub00 = sub0.split(':')[1] # split on :, need to get rid of chrom

		try:
			lPosCurr = sub00
			rPosCurr = sub1
			# rPosQuery and lPosQuery are GLOBALs
			if (lPosCurr >= lPosQuery) & (lPosCurr <= rPosQuery): # left position good
				if (rPosCurr >= lPosQuery) & (rPosCurr <= rPosQuery): # right position good
					match = 1
		except IndexError:
			print('index error')

	return match

#////////////////////////////////////////////////////////////////////
# hitSearchFunc_coords()
#	given a list of shared entries between an individual cell's VCF
# 	and the COSMIC LAUD db, searches for hits to the GOI, as specified
#	with cmd line option '4'
#
#		REMEMBER: `match` is NOT a bool here
#
# 	passing in *args so that i can tell what cell im finding indels in!!!
#////////////////////////////////////////////////////////////////////
def hitSearchFunc_coords(sample, *args):
	cell_ = args[0]
	match = ""

	currChrom = sample.split(':')[0]
	if currChrom == queryChrom:
		sub0 = sample.split('-')[0] # split on `-`
		sub1 = sample.split('-')[1] # this guy is good
		sub00 = sub0.split(':')[1] # split on :, need to get rid of chrom

		try:
			lPosCurr = sub00
			rPosCurr = sub1
			# keep in mind rPosQuery and lPosQuery are GLOBALs
			if (lPosCurr >= lPosQuery) & (lPosCurr <= rPosQuery): # left pos GOI match
				if (rPosCurr >= lPosQuery) & (rPosCurr <= rPosQuery): # right pos GOI match
					if lPosCurr == rPosCurr: # SNP
						match = lPosCurr
					else: 		# found an indel!!
						match = lPosCurr + '-' + rPosCurr
						#print(cell_)

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
		cell = f.replace("vcf_old_test/", "")
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
		cell = f.replace("vcf_old_test/", "")
		cell = cell.replace(".vcf", "")	

		df = VCF.dataframe(f)
		genomePos_query = df.apply(getGenomePos, axis=1) # apply function for every row in df
		# get the entries shared between curr cells VCF and the LAUD filter set
		#	remember, these are general, and NOT gene specific
		genomePos_query_expand = expandSet(set(genomePos_query))

		shared = list(set(genomePos_query_expand) & set(genomePos_laud_db))
		shared1 = pd.Series(shared) # convert to pandas obj
		matches = shared1.apply(hitSearchFunc_coords, args=(cell,)) # another apply call 

		# delete empty dict keys
		for k in matches.keys():
			try:
				if len(matches[k])<1:
					del matches[k]
			except: pass

		cells_dict_GOI_coords.update({cell : list(matches.values)})

	return cells_dict_GOI_coords

#////////////////////////////////////////////////////////////////////
# getMutationAA()
#	Pass in a dict of {cell, list(genomePos)} items and it returns a
#	dict of {cell, list(Mutation.AA)}. going BACK to database_laud here
# 
#////////////////////////////////////////////////////////////////////
def getMutationAA(d, chr):
	print('AA searching...')
	newDict = {}

	for k in d:
		valuesList = d.get(k) # can now handle values with multiple entries
		newValues = []

		for entry in valuesList:
			testSplit = entry.split('-') # if its a SNP it wont have '-' at all	

			### CASE 1 -- SNP
			if len(testSplit) == 1:
				chrStr = chr + ':' + entry + '-' + entry
				filter_df = database_laud[database_laud["Mutation genome position"].str.contains(chrStr)==True]
				currMuts = filter_df['Mutation AA']

				for item in currMuts:		# really shouldnt have a for loop here
					item = item.replace("p.", "")

				newValues.append(item) 		# effectively just taking the last item in the list
			
			### CASE 2 -- INDEL 
			else:
				chrStr = chr + ':' + entry
				filter_df = database_laud[database_laud["Mutation genome position"].str.contains(chrStr)==True]
				currMuts = filter_df['Mutation AA']

				for item in currMuts:		# really shouldnt have a for loop here
					item = item.replace("p.", "")

				newValues.append(item)		# effectively just taking the last item in the list

		newDict.update({k : newValues})

	return newDict

#////////////////////////////////////////////////////////////////////
# expandSet()
#	Pass in a set of genome coords, and it will 'expand' the indels
# 	within that set by adding +/- 3 bp copies for each one
#
#////////////////////////////////////////////////////////////////////
def expandSet(mySet):
	returnSet = []

	for entry in mySet:
		l0 = []
		l1 = []
		try:
			sub0 = entry.split('-')[0] # split on `-`
			sub1 = entry.split('-')[1] # this guy is good
			sub00 = sub0.split(':')[1] # split on :, need to get rid of chrom
			chrom = sub0.split(':')[0]
		
			if sub00 != sub1: # got an indel 
				sub00_1 = int(sub00) + 1
				sub00_2 = int(sub00) + 2
				sub00_3 = int(sub00) + 3
				sub00_4 = int(sub00) - 1
				sub00_5 = int(sub00) - 2
				sub00_6 = int(sub00) - 3

				l0.extend((sub00_1, sub00_2, sub00_3, sub00_4, sub00_5, sub00_6))
				
				try:
					sub1_1 = int(sub1) + 1
					sub1_2 = int(sub1) + 2
					sub1_3 = int(sub1) + 3
					sub1_4 = int(sub1) - 1
					sub1_5 = int(sub1) - 2
					sub1_6 = int(sub1) - 3

					l1.extend((sub1_1, sub1_2, sub1_3, sub1_4, sub1_5, sub1_6))
				
				except ValueError:
					continue

				coord_combos = list(itertools.product(l0, l1))
				for pair in coord_combos:
					toAdd = chrom + ':' + str(pair[0]) + '-' + str(pair[1])
					returnSet.append(toAdd)

			else:
				returnSet.append(entry)
		
		except IndexError:
			continue

	return returnSet

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
#			python3 getMutationCounts_overall_and_GOI.py 4 7 55019101 55211628 egfr_out.csv
#////////////////////////////////////////////////////////////////////

global database
global database_laud

if len(sys.argv) == 1:
	print('usage: python3 getMutationCounts_overall_and_GOI.py [-h] [1] [2] [3] [4]')
	print('  ')
	print('		1 - getRawCounts')
	print('		2 - getFilterCountsBasic')
	print('		3 - getFilterCountsLAUD')
	print('		4 - getGeneOfInterest')
	print('			needs [chrom] [pos1] [pos2] [outFile]')
	print('			ie. python3 getMutationCounts_overall_and_GOI.py 4 7 500050 50010 myOut.csv')
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
	database = pd.read_csv("CosmicGenomeScreensMutantExport.tsv", delimiter = '\t')
	fNames = getFileNames()
	filterDict = getFilterCountsBasic(fNames)
	print("filter counts (basic) done!")
	print('writing csv')
	writeCSV(filterDict, "nonImmune_GATK_hits_COSMIC_filter.csv")

# filter counts LAUD
if sys.argv[1] == '3':
	print('setting up COSMIC database...')
	database = pd.read_csv("CosmicGenomeScreensMutantExport.tsv", delimiter = '\t')
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
		print('	ie. python3 getMutationCounts_overall_and_GOI.py 4 7 500000 50010 myFile.csv')
		print(' ')
		sys.exit()
	
	print('setting up COSMIC database...')
	database = pd.read_csv("CosmicGenomeScreensMutantExport.tsv", delimiter = '\t')
	database_laud = getLAUD_db()
	#database_laud.to_csv('database_laud.csv') # test write
	fNames = getFileNames()
	
	chromo = sys.argv[2]
	position1 = sys.argv[3]
	position2 = sys.argv[4]	

	#goiDict = getGOIHits(fNames, chromo, position1, position2) # standard call - get raw counts
	goiDict = getGOIHit_coords(fNames, chromo, position1, position2) # get genome coords
	print("GOI search done!")
	
	outFilePref = sys.argv[5]
	writeCSV(goiDict, outFilePref + '.csv')

	goiDict_AA = getMutationAA(goiDict, chromo)
	print('AA search done')
	writeCSV(goiDict_AA, outFilePref + '_AA.csv')

#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
