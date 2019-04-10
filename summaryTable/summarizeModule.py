# summarizeLib.py
# Lincoln Harris
# 3.28.19
# 
# module that allow you to create per-cell / per-sample summary tables

import pandas as pd
import numpy as np
import math

# mutationsDF__fillIn()
#    goal is to construct a cell-wise dataframe with mutations to each
#    of EGFR, KRAS and BRAF. the challange is getting the cells to line
#    up, hence the for loop 
#
#    GOI needs to be lowercase
#
def mutationsDF_fillIn(GOI, GOI_df, mutationsDF_):
	mutName = GOI + 'Mut'
	for i in range(0,len(mutationsDF_.index)):
		currCell = mutationsDF_['cell'][i]

		rightIndex = GOI_df['cell'] == currCell
		rightRow = GOI_df[rightIndex]
    
		rightCell = rightRow['cell']
		rightCell = str(rightCell).split()[1]
    
		rightMut = rightRow['mutations']
		rightMut = str(rightMut).split()[1]
    
		mutationsDF_[mutName][i] = rightMut


# removeExtraCharacters_mutationsDF_()
#    essentially converting mutationsDF_ mutation cols from lists to 
#    strings. makes downstream analysis easier
#
#    GOI needs to be lowercase
#
def removeExtraCharacters_mutationsDF(GOI, mutationsDF_):
	mutName = GOI + 'Mut'

	mutationsDF_[mutName] = mutationsDF_[mutName].str.replace("'", "") # remove quotes
	mutationsDF_[mutName] = mutationsDF_[mutName].str.replace("[", "") # remove brackets
	mutationsDF_[mutName] = mutationsDF_[mutName].str.replace("]", "") # remove brackets
	mutationsDF_[mutName] = mutationsDF_[mutName].str.replace(" ", "") # remove whitespace?


# genericSummaryTableFillIn()
#    fills in a given (metadata) field in summaryTable_. pulls from 
#    patientMetadata_ and goes cell-by-cell through 
#    summaryTable_, filling in fields like patientID/driver_gene
#
def genericSummaryTableFillIn(metaField, summaryField, summaryTable_, patientMetadata_):
	for i in range(0,len(summaryTable_.index)):
		currCell = summaryTable_['cell'].iloc[i]
		currPlate = currCell.split('_')[1]
    
		index_to_keep = patientMetadata_['plate'] == currPlate
		keepRow = patientMetadata_[index_to_keep]
		try:
			currField = list(keepRow[metaField])[0]
			summaryTable_[summaryField][i] = currField
		except IndexError:
			continue
			#print('ERROR: plate not found') # these are just the plates were NOT 
        	                                 # including in the analysis


# fusionsFillIn()
#    Takes the existing fusionsDF (which is just a list of the five fusions
#    we looked for, and what cells they're found in) and populates 
#    summaryTable_ with this shit
#
#    this works, but holllllyyyy shitttt we can do better
#
def fusionsFillIn(fusionsDF_, summaryTable_):
	for i in range(0, len(summaryTable_.index)):
		currCell = summaryTable_['cell'][i]
		fusionsListCurr = []
    
		colList0 = list(fusionsDF_['ALK--EML4'])
		colList1 = list(fusionsDF_['ALK_any'])
		colList2 = list(fusionsDF_['EML4_any'])
		colList3 = list(fusionsDF_['NTRK_any'])
		colList4 = list(fusionsDF_['RET_any'])
		colList5 = list(fusionsDF_['ROS1_any'])

		if currCell in colList0:
			fusionsListCurr.append('ALK-EML4')
		elif currCell in colList1:
			fusionsListCurr.append('ALK_any')
		elif currCell in colList2:
			fusionsListCurr.append('EML4_any')
		elif currCell in colList3:
			fusionsListCurr.append('NTRK_any')
		elif currCell in colList4:
			fusionsListCurr.append('RET_any')
		elif currCell in colList5:
			fusionsListCurr.append('ROS1_any')
		else:
			fusionsListCurr = ""
        
		fusionsListCurr = str(fusionsListCurr)
		fusionsListCurr = fusionsListCurr.strip(']')
		fusionsListCurr = fusionsListCurr.strip('[')
		fusionsListCurr = fusionsListCurr.strip("'")
		fusionsListCurr = fusionsListCurr.strip(" ")
 
		summaryTable_['fusions_found'][i] = fusionsListCurr


# translatedMutsFillIn_EGFR()
#    need to make a 'mutations_found_translated' field that converts our
#    'raw' mutation calls to something that more resembles those reported
#    in our clinical cols. Need a seperate func for EGFR, bc there are 
#    so many potential variants to account for
#
def translatedMutsFillIn_EGFR(summaryTable_):
	for i in range(0,len(summaryTable_.index)):
		translatedList = []
		currCell = summaryTable_['cell'].iloc[i]
		currMuts_egfr = summaryTable_['mutations_found_EGFR'].iloc[i]
		currMuts_egfr_split = currMuts_egfr.split(',')
		for item in currMuts_egfr_split:
			if 'delELR' in item:
				translatedList.append('EGFR del19')
			elif '745_' in item:
				translatedList.append('EGFR del19')
			elif '746_' in item:
				translatedList.append('EGFR del19')
			elif 'ins' in item:
				translatedList.append('EGFR ins20')
			elif item != '':
				translatedList.append('EGFR ' + item)
        
		summaryTable_['mutations_found_translated'][i] = translatedList


# translatedMutsFillIn_nonEGFR()
#    need to make a 'mutations_found_translated' field that converts our
#    'raw' mutation calls to something that more resembles those reported
#    in our clinical cols. This func handles BRAF and KRAS, bc there are
#    only like 2 possible clinically reported muts for them, so we'd might
#    as well keep everything
#
#    want GOI to be capitilized here
def translatedMutsFillIn_nonEGFR(GOI, summaryTable_):
	colName = 'mutations_found_' + GOI
	for i in range(0,len(summaryTable_.index)):
		translatedList = []
		currCell = summaryTable_['cell'].iloc[i]
		currMuts = summaryTable_[colName].iloc[i]
		currMuts_split = currMuts.split(',')
		for item in currMuts_split:
			if item != '' and '?' not in item:
				translatedList.append(GOI + ' ' + item)

		summaryTable_['mutations_found_translated'][i] = summaryTable_['mutations_found_translated'][i] + translatedList


# translatedMutsFillIn_fusions()
# 	 need to make a 'mutations_found_translated' field that converts our
#    'raw' mutation calls to something that more resembles those reported
#    in our clinical cols. for fusions this time
#
def translatedMutsFillIn_fusions(summaryTable_):
	for i in range(0,len(summaryTable_.index)):
		translatedList = []
		currCell = summaryTable_['cell'].iloc[i]
		currFus = summaryTable_['fusions_found'].iloc[i]
		currFus_split = currFus.split(',')
		for item in currFus_split:
			if item == 'ALK-EML4':
				translatedList.append('ALK fusion')
				translatedList.append('EML4 fusion')
				translatedList.append('ALK-EML4 fusion')
			elif item != '' and '?' not in item:
				item = item.split('_')[0]
				translatedList.append(item + ' fusion')

		summaryTable_['mutations_found_translated'][i] = summaryTable_['mutations_found_translated'][i] + translatedList


# convertToString()
#    really just taking this mutations_found_translated col and converting
#    it from a list to a string. makes taking set() easier, but since
#    this is a script now, maybe i dont even need this 
#
def convertToString(summaryTable_):
	for i in range(0,len(summaryTable_.index)):
		currStr = str(summaryTable_['mutations_found_translated'][i])
		currStr = currStr.replace("'", "")
		currStr = currStr.replace("]", "")
		currStr = currStr.replace("[", "")
		summaryTable_['mutations_found_translated'][i] = currStr


# clinMutFound_fillIn()
#    want to fill in this clin_mut_found_bool col with 1 if the clinically
#    reported mutation is found, 0 if else
#
def clinMutFound_fillIn(summaryTable_):
	for i in range(0,len(summaryTable_.index)):
		currCell = summaryTable_['cell'][i]
		currMuts = summaryTable_['mutations_found_translated'][i]
		currClinGene = summaryTable_['clinical_driver_gene'][i]
		currClinMut = summaryTable_['clinical_mutation'][i]
		currClinMut_str = str(currClinGene) + ' ' + str(currClinMut)
    
		if currClinMut_str in currMuts:
			summaryTable_['clin_mut_found_bool'][i] = 1
		else:
			summaryTable_['clin_mut_found_bool'][i] = 0


# tumorCellBoolFillIn()
#    want to fill in this tumorCell_bool with 1 if we're calling that
#    cell a tumor cell in our seurat obj, 0 if else
#
def tumorCellBoolFillIn(summaryTable_):
	# NEED TO READ IN SEURAT METADATA, SO WE CAN SET tumorCell_bool
	metaPATH = '/Users/lincoln.harris/Desktop/LAUD_important_shit/metadataSeurat.csv'
	metadataSeurat = pd.read_csv(metaPATH)

	myCols = list(metadataSeurat.columns)
	myCols[0] = 'cell'
	metadataSeurat.columns = myCols
	
	indicies = metadataSeurat['inferCNV_annotation'] == 'perturbed'
	metadataSeurat_pert = metadataSeurat[indicies]
	
	tumorCellsList = list(metadataSeurat_pert['cell'])

	# now fill in 'tumorCell_bool' for summaryTable_
	for i in range(0, len(summaryTable_.index)):
		currCell = summaryTable_['cell'][i]
		if currCell in tumorCellsList:
			summaryTable_['tumorCell_bool'][i] = 1
		else:
			summaryTable_['tumorCell_bool'][i] = 0


# getNonZeroCovROI()
#    takes a given coverageByCell dataframe and filters for the non-zero 
#    vals. coverage dfs come from checkCoverage_parallel.py
#
def getNonZeroCovROI(gene, mut):
	fPATH = '/Users/lincoln.harris/code/SNP_calling_pipeline/coverage/out/' + gene + '_' + mut + '_' + 'coverageByCell.csv'
	cov = pd.read_csv(fPATH)
	indices = cov['depth_gvcf'] != 0
	cov_nonZero = cov[indices]

	return(cov_nonZero)


# ROI_coverage_fillIn()
#    fills in coverage for a given ROI, for summaryTable_
#
def ROI_coverage_fillIn(coverage_df, queryGene, queryMutation, summaryTable_):
    for i in range(0, len(summaryTable_.index)):
        currCell = summaryTable_['cell'][i]
        currDriver = summaryTable_['clinical_driver_gene'][i]
        currMut = summaryTable_['clinical_mutation'][i]
    
        if currDriver == queryGene and currMut == queryMutation:
            if currCell in list(coverage_df['cellName']):
                index_cov_nonZero = coverage_df['cellName'] == currCell
                currRow_cov_nonZero = coverage_df[index_cov_nonZero]
                currDepth_gvcf = int(currRow_cov_nonZero['depth_gvcf'])
        
                summaryTable_['coverage_to_ROI'][i] = currDepth_gvcf
            else:
                summaryTable_['coverage_to_ROI'][i] = 0


# validationTable_metadata_fillIn()
#    fills in metadata field for the validationTable
#              
def validationTable_metadata_fillIn(metaField, validationField, validationTable_, patientMetadata_):
	for i in range(0, len(validationTable_.index)):
		currSample = validationTable_['sample'][i]
		try:
			rowToKeep = patientMetadata_['sample_name'] == currSample
			patientRows = patientMetadata_[rowToKeep] # will return MULTIPLE rows
			patientRows = patientRows.reset_index(drop=True)

			fillField = patientRows[metaField][0]
       
			validationTable_[validationField][i] = fillField
		except:
			continue
			#print('ERROR')


# validationTable_dict_muts()
#    returns a dictionary that holds values for all of the
#    mutations to a given cell. 
# 
def validationTable_dict_muts(validationTable_, summaryTable_):
	d = {}
	samplesList = validationTable_['sample']

	for item in samplesList:
		d.update({item:''})

	for i in range(0, len(summaryTable_.index)):
		currSample = summaryTable_['sample_name'][i]
		currMuts = summaryTable_['mutations_found'][i]
		currMuts = str(currMuts)
		currMutsSplit = currMuts.split(',')

		currDictVal = d[currSample]
    
		for item in currMutsSplit:
			if item not in currDictVal and item != 'nan':
				updateVal = currDictVal + item + ', '
				d.update({currSample:updateVal})

	return(d)


# validationTable_dict_generic()
#    returns a dict with values for num cells that are tumor/
#    have coverage to a given ROI, depending on what field 
#    value is passed in.
#         
def validationTable_dict_generic(validationTable_, summaryTable_, field):
	d = {}
	samplesList = validationTable_['sample']
	for item in samplesList:
		d.update({item:0})

	for i in range(0, len(summaryTable_.index)):
		currSample = summaryTable_['sample_name'][i]
		currBool = summaryTable_[field][i]

		currDictVal = d[currSample]  

		if not math.isnan(currBool) and currBool != 0:
			updateVal = currDictVal + 1
			d.update({currSample:updateVal})

	return(d)
    