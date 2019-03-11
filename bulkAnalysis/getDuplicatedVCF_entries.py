#/////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////
# author: Lincoln
# date: 3/11/19
# script: getDuplicatedVCF_entries.py
# 
# can we do the germline filter, for all the 'germline' variants for
# a given patient? i think we're pretty close here
# 		turning our python notebook into a script here
#
# lets just write this for a single patient and see how things go
#
# usage:
#    python3 getDuplicatedVCF_entries.py [patientID]
#
#/////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////
import pandas as pd
import VCF
import os
import shutil
import sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

#////////////////////////////////////////////////////////////////////
# getCellsList
#	get the list of cell names from a given patient
#////////////////////////////////////////////////////////////////////
def getPatientCellsList(scVCF_list_, patientID):
	currPatient_cells_ = []

	for item in scVCF_list:
		currCell = item.strip('.vcf')
		currPlate = currCell.split('_')[1]
    
		rowToKeep = patientMetadata['plate'] == currPlate
    
		try:
			currPatient = patientMetadata['patient_id'][rowToKeep]
			currPatientVal = currPatient.item()

			if currPatientVal == patientID:
				currPatient_cells_.append(currCell)
		except:
			continue

	return currPatient_cells_

#////////////////////////////////////////////////////////////////////
# getUniqueVCF_entries()
#     do the germline filter, and return a dataframe with only the
#     UNIQUE entries for a given cell 
#////////////////////////////////////////////////////////////////////
def getUniqueVCF_entries(patient, cell):
	basePATH = '/Users/lincoln.harris/code/SNP_calling_pipeline/bulkAnalysis/'
	patientPATH = basePATH + 'bulkVCF/' + patient + '.vcf'
	cellPATH = basePATH + 'scVCF/' + cell + '.vcf'
    
	patient_df = VCF.dataframe(patientPATH)
	cell_df = VCF.dataframe(cellPATH)
    
	patient_df_trimmed = patient_df[['CHROM', 'POS', 'ID', 'REF', 'ALT']]
	cell_df_trimmed = cell_df[['CHROM', 'POS', 'ID', 'REF', 'ALT']]
    
	# get whats SHARED between patient and cell 
	#    FIND GERMLINE MUTATIONS
	patient_cell_concat = pd.concat([patient_df_trimmed, cell_df_trimmed])
	rowsToKeep = patient_cell_concat.duplicated()
	patient_cell_shared = patient_cell_concat[rowsToKeep]
	patient_cell_shared = patient_cell_shared.reset_index(drop=True)
    
	# now go back to the original cell df, pull out whats UNIQUE 
	#     THIS IS THE GERMLINE FILTER!!
	cell_cell_concat = pd.concat([cell_df_trimmed, patient_cell_shared])
	cell_cell_concat_noDups = cell_cell_concat.drop_duplicates(keep=False)
	cell_cell_concat_noDups = cell_cell_concat_noDups.reset_index(drop=True)
    
	return(cell_cell_concat_noDups)

#////////////////////////////////////////////////////////////////////
# main()
#	get the patient name from the cmdline, set up patientMetadata, 
#   call subroutines to get list of per-patient cells, then call the 
#   filtering func and write output to new csv
#////////////////////////////////////////////////////////////////////

global patientMetadata

# parse cmdline arg
currPatient = sys.argv[1]

# read in patient metadata
patientMetadata = pd.read_csv('../cDNA_plate_metadata.csv')

# get a list of all the single-cell VCF files
wrkDir = '/Users/lincoln.harris/code/SNP_calling_pipeline/bulkAnalysis/scVCF/'
scVCF_list = os.listdir(wrkDir)

currPatient_cells = getPatientCellsList(scVCF_list, currPatient)

for currCell in currPatient_cells:
	currCell_unique = getUniqueVCF_entries(currPatient, currCell)
	outStr = './filteredOut/' + currCell + '_unique.csv'
	currCell_unique.to_csv(outStr, index=False)

#/////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////