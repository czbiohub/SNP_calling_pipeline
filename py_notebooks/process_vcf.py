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
import VCF
import os
import csv
import pandas as pd

# getFileNames()
#	Get file names based on the specified path
#

def getFileNames():
	files = []
	for file in os.listdir("/home/ubuntu/expansionVol2/04-GATK_out/all"):
		if file.endswith(".vcf"):
			fullPath = (os.path.join("/home/ubuntu/expansionVol2/04-GATK_out/all", file))
			files.append(fullPath)
    
	return files

# getRawCounts()
#	Creates dictionary obj with raw counts for GATK hits w/in a given set of vcf files
#

def getRawCounts(fileNames):
	cells_dict = {}

	for f in fileNames:
		cell = f.replace("/home/ubuntu/expansionVol2/04-GATK_out/all/", "")
		cell = cell.replace(".vcf", "")
    
		df = VCF.dataframe(f)
		unique = len(np.unique(df.POS))
    
		cells_dict.update({cell : unique})

	return cells_dict

# getGenomePos()
#	Returns a genome position sting that will match against the ones w/in COSMIC db
#

def getGenomePos(sample):
	chr = sample[0]
	chr = chr.replace("chr", "")
	pos = sample[1]
	genomePos = chr + ':' + str(pos) + '-' + str(pos)

	return(genomePos)

# getFilterCountsBasic()
#	Creates dictionry obj with COSMIC filtered GATK hits w/in a given set of vcfs 
#

def getFilterCountsBasic(fileNames):
	cells_dict_filter = {}
	genomePos_db = pd.Series(database['Mutation genome position'])

	for f in fileNames:
		cell = f.replace("/home/ubuntu/expansionVol2/04-GATK_out/all/", "")
		cell = cell.replace(".vcf", "")
		print(cell)
		df = VCF.dataframe(f)
		nrow = df.shape[0]
		genomePos_query = []
    
		for i in range(1, nrow):
			found = False
			currRow = df.iloc[i,]
			toAdd = getGenomePos(currRow)
			genomePos_query.append(toAdd)
        
		shared = list(set(genomePos_query) & set(genomePos_db))
		cells_dict_filter.update({cell : len(shared)})
    
		#print(cells_dict_filter)
	return cells_dict_filter

# writeCSV()
#	Writes the contents of a dictionary object to a csv
#

def writeCSV(dictObj, outFile):
	with open(outFile, 'w') as csv_file:
		writer = csv.writer(csv_file)
		for key, value in dictObj.items():
			writer.writerow([key, value])


# main()
#	Main logic here. Comment out code blocks depending on which output file you want, nonImmune_GATK_hits_raw.csv, 
#	nonImmune_GATK_hits_COSMIC_filter.csv, or nonImmune_GATK_hits_COSMIC_filter_adv.csv

database = pd.read_csv("/home/rstudio/17-variantCallingR/05-cosmicDB/CosmicGenomeScreensMutantExport.tsv", delimiter = '\t')
fNames = getFileNames()

#rawDict = getRawCounts(fNames)
#print("raw counts done!")
#writeCSV(rawDict, "nonImmune_GATK_hits_raw.csv")

filterDict = getFilterCountsBasic(fNames)
print("filter counts (basic) done!")
writeCSV(filterDict, "nonImmune_GATK_hits_COSMIC_filter.csv")

#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
