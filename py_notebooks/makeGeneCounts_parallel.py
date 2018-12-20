#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
# script: makeGeneCounts.py
# author: Lincoln 
# date: 12.18.18
#
# Creates a gene/cell table for the mutations found in a given 
# population of cells. Trying to implement parallelization here
#
# usage: 
#			python3 makeGeneCounts_parallel.py
#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
import vcf
import numpy as np
import VCF # comes from Kamil Slowikowski
import os
import csv
import pandas as pd
import sys
import multiprocessing as mp

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
# getGeneName()
#	want to return the gene name from a given genome position string
#   (ie. '1:21890111-21890111'), by querying the hg38-plus.gtf
#
#////////////////////////////////////////////////////////////////////
def getGeneName(posString):
	# work on posString
	chrom = posString.split(':')[0]
	posString_remove = posString.split(':')[1]
	lPosition = posString_remove.split('-')[0] 
	rPosition = posString_remove.split('-')[1] 

	# work on hg38_gtf
	chromStr = 'chr' + str(chrom)
	hg38_gtf_filt = hg38_gtf.where(hg38_gtf[0] == chromStr).dropna()
	hg38_gtf_filt = hg38_gtf_filt.where(hg38_gtf_filt[3] <= int(lPosition)).dropna() # lPos good
	hg38_gtf_filt = hg38_gtf_filt.where(hg38_gtf_filt[4] >= int(rPosition)).dropna() # rPos good
	
	try:
		returnStr = str(hg38_gtf_filt.iloc[0][8])	# keep just the gene name / meta data col
		returnStr = returnStr.split(';')[1]
		returnStr = returnStr.strip(' gene_name')
		returnStr = returnStr.strip(' ')
		returnStr = returnStr.strip('"')
	except IndexError:
		returnStr = ''
	
	return returnStr

#////////////////////////////////////////////////////////////////////
# getGeneCellCounts()
#	Creates dictionry obj where every key is a cell and every value is
#	a list of the genes we found mutations in for that cell. 
#////////////////////////////////////////////////////////////////////
def getGeneCellMutCounts(f):
	tup = [] # not really a tuple, just a list, i guess

	cell = f.replace("../vcf_test/", "")
	cell = cell.replace(".vcf", "")
	print(cell) # to see where we are
	
	df = VCF.dataframe(f)
	genomePos_query = df.apply(getGenomePos, axis=1) # apply function for every row in df

	# this solution retains duplicates
	items = set(genomePos_query) # genomePos_query (potentially) has dups
	shared = [i for i in genomePos_laud_db if i in items]

	shared_series = pd.Series(shared)
	sharedGeneNames = shared_series.apply(getGeneName)

	tup = [cell, sharedGeneNames]

	return(tup)

#////////////////////////////////////////////////////////////////////
# formatDataFrame()
#	logic for creating the cell/mutation counts table from the raw 
#	output that getGeneCellMutCounts provides
#////////////////////////////////////////////////////////////////////
def formatDataFrame(raw_df):
	cellNames = list(raw_df.index)

	genesList = []
	for i in range(0, raw_df.shape[0]):
		currList = list(raw_df.iloc[i].unique()) # unique genes for curr_cell 

		for elm in currList:	
			if elm not in genesList:
				genesList.append(elm)

	genesList1 = pd.Series(genesList)

	df = pd.DataFrame(columns=genesList1, index=cellNames) # initialize blank dataframe
	for col in df.columns: # set everybody to zero
		df[col] = 0

	for i in range(0,raw_df.shape[0]): # loop over raw dataframe
		currCell = raw_df.index[i]
		currRow = raw_df.iloc[i]

		for currGene in currRow:	# find corresponding entry in formatted dataframe
			df[currGene][currCell] += 1		# add

	return(df)

#////////////////////////////////////////////////////////////////////
# main()
#	
#////////////////////////////////////////////////////////////////////
global database
global database_laud
global hg38_gtf
global genomePos_laud_db

database = pd.read_csv("../CosmicGenomeScreensMutantExport.tsv", delimiter = '\t')
database_laud = getLAUD_db()
genomePos_laud_db = pd.Series(database_laud['Mutation genome position'])
hg38_gtf = pd.read_csv('../hg38-plus.gtf', delimiter = '\t', header = None)
fNames = getFileNames()

print('creating pool')

p = mp.Pool(processes=12)

try:
	cells_list = p.map(getGeneCellMutCounts, fNames, chunksize=1) # default chunksize=1
finally:
	p.close()
	p.join()

# convert to dictionary
cells_dict = {}

for item in cells_list:
    cells_dict.update({item[0]:item[1]})

print('writing file')

filterDict_pd = pd.DataFrame.from_dict(cells_dict, orient="index") # orient refers to row/col orientation 
filterDict_format = formatDataFrame(filterDict_pd)
filterDict_format.to_csv("foo.csv")

#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////