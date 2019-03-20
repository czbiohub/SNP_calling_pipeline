#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
# script: makeGeneCounts.py
# author: Lincoln 
# date: 12.18.18
#
# Creates a gene/cell table for the mutations found in a given 
# population of cells. Trying to implement parallelization here
# Run this on a big-ass machine
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
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

#////////////////////////////////////////////////////////////////////
# getGermlineFilteredCellsList()
#	just returns a list of the cells we already have germline filter
#    for
#////////////////////////////////////////////////////////////////////
def getGermlineFilteredCellsList():
	filterDir = '/home/ubuntu/code/SNP_calling_pipeline/bulkAnalysis/filteredOut/'
	filterDir_list = os.listdir(filterDir)

	filteredCells = []
	for f in filterDir_list:
		cell = f.strip('_unique.vcf')
		filteredCells.append(cell)

	return filteredCells

#////////////////////////////////////////////////////////////////////
# getFileNames()
#	Get file names based on the specified path
#
#////////////////////////////////////////////////////////////////////
def getFileNames():
	files = []
	for file in os.listdir("vcf_test/"):
		if file.endswith(".vcf"):
			fullPath = (os.path.join("vcf_test/", file))
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
	try:
		chr = str(sample[0])
		chr = chr.replace("chr", "")
		pos = int(sample[1])
		ref = str(sample[3])
		alt = str(sample[4])
	
		if (len(ref) == 1) & (len(alt) == 1): # most basic case
			secondPos = pos
			genomePos = chr + ':' + str(pos) + '-' + str(secondPos)
		elif (len(ref) > 1) & (len(alt) == 1):
			secondPos = pos + len(ref)
			genomePos = chr + ':' + str(pos) + '-' + str(secondPos)
		elif (len(alt) > 1) & (len(ref) == 1):
			secondPos = pos + len(alt)
			genomePos = chr + ':' + str(pos) + '-' + str(secondPos)
		else: # multibase-for-multibase substitution
			secondPos = '1'
			genomePos = chr + ':' + str(pos) + '-' + str(secondPos)
	except:
		genomePos = 'ERROR'

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
#
#////////////////////////////////////////////////////////////////////
def getGeneCellMutCounts(f):
	tup = [] # not really a tuple, just a list, i guess

	cell = f.replace("vcf_test/", "")
	cell = cell.replace(".vcf", "")
	print(cell) # to see where we are
	
	df = VCF.dataframe(f)
	genomePos_query = df.apply(getGenomePos, axis=1) # apply function for every row in df

	items = set(genomePos_query) # genomePos_query (potentially) has dups

	# COSMIC filter
	shared = [i for i in genomePos_laud_db if i in items] # retains dups

	shared_series = pd.Series(shared)
	sharedGeneNames = shared_series.apply(getGeneName)
	tup = [cell, sharedGeneNames]

	return(tup)

#////////////////////////////////////////////////////////////////////
# formatDataFrame()
#	logic for creating the cell/mutation counts table from the raw 
#	output that getGeneCellMutCounts provides
#
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

	for i in range(0,raw_df.shape[0]):
		currCell = raw_df.index[i]
		currRow = raw_df.iloc[i]

		for currGene in currRow:
			df[currGene][currCell] += 1

	return(df)

#////////////////////////////////////////////////////////////////////
# main()
#	
#////////////////////////////////////////////////////////////////////
global database
global database_laud
global hg38_gtf
global genomePos_laud_db
global germlineFilterCells

database = pd.read_csv("CosmicGenomeScreensMutantExport.tsv", delimiter = '\t')
database_laud = getLAUD_db()
genomePos_laud_db = pd.Series(database_laud['Mutation genome position'])
hg38_gtf = pd.read_csv('hg38-plus.gtf', delimiter = '\t', header = None)
fNames = getFileNames()

germlineFilterCells = getGermlineFilteredCellsList()

print('creating pool')

p = mp.Pool(processes=14)

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
filterDict_format.to_csv("geneCellMutationCounts_test_toggle.csv")

#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////

#if cell in germlineFilterCells: # DONT do cosmic filter
#	shared = list(items)
#	print('GERMLINE FILTER. length hits list: %d' % len(shared))
#else:  # DO cosmic filter
#	shared = [i for i in genomePos_laud_db if i in items] # retains dups
#	print('COSMIC FILTER. length hits list: %d' % len(shared))
