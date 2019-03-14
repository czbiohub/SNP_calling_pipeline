#/////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////
# author: Lincoln
# date: 3/11/19
# script: getDuplicatedVCF_entries.py
# 
# can we do the germline filter, for all the 'germline' variants for
# a given patient? i think we're pretty close here
# 		turning our jupyter notebook into a script
#
# trying to write this for ALL patients now
#
# usage:
#    python3 getDuplicatedVCF_entries.py
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

	for item in scVCF_list_:
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
	print('numCells: %d' % len(currPatient_cells_))
	return currPatient_cells_

#////////////////////////////////////////////////////////////////////
# getUniqueVCF_entries()
#     do the germline filter, and return a dataframe with only the
#     UNIQUE entries for a given cell 
#////////////////////////////////////////////////////////////////////
def getUniqueVCF_entries(patient, cell):
	basePATH = os.getcwd()
	#patientPATH = basePATH + '/bulkVCF/' + patient + '.vcf'
	patientPATH = basePATH + '/bulkVCF/' + patient
	cellPATH = basePATH + '/scVCF/' + cell + '.vcf'
	try:
		patient_df = VCF.dataframe(patientPATH)
		cell_df = VCF.dataframe(cellPATH)
	except FileNotFoundError:
		print('FILE NOT FOUND: %s' % cellPATH)
		return
    
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

# read in patient metadata
patientMetadata = pd.read_csv('../cDNA_plate_metadata.csv')

# get a list of all the single-cell VCF files
cwd = os.getcwd()
vcfDir = cwd + '/scVCF/'
scVCF_list = os.listdir(vcfDir)

# get list of bulk VCF files
bulkVCF_dir = cwd + '/bulkVCF/'
bulkVCF_list = os.listdir(bulkVCF_dir)

patientsRun = [] # need to keep track of which patients have been run

# outer loop -- by PATIENT
for item in bulkVCF_list:
	currSample = item.strip('.vcf')
	currPatient = currSample.split('_')[0]
	suffix1 = currSample.split('_')[1]
	try:
		suffix2 = currSample.split('_')[2]
	except IndexError:
		suffix2 = ''
	
	if suffix2 != '' and currPatient not in patientsRun:
		print('WHOLE BLOOD FOUND, for %s' % currPatient)
		currPatient_cells = getPatientCellsList(scVCF_list, currPatient)

		# inner loop -- by CELL 
		for currCell in currPatient_cells:
			currCell_unique = getUniqueVCF_entries(item, currCell)
			outStr = './filteredOut/' + currCell + '_unique.csv'
			currCell_unique.to_csv(outStr, index=False)
			#continue
			
		patientsRun.append(currPatient)

#/////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////
