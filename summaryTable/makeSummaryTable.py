#/////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////
# script: makeSummaryTable.py
# author: Lincoln
# date: 3.28.19
#
# Takes as input GOI_out_AA.csv files (from getMutationCounts_overall_and_GOI.py),
# patient metadata, seurat metadata, fusionsDF, and creates a BY CELL 
# summaryTable. The goal with this table is to provide an answer to questions like
# 'which patients have which mutations?', and 'how many cells have clinically relevant
# mutations?'. Currently i've got an ipynb that accomplishes this, but its like
# super long and unwieldy, so i though converting it to a script would be a good
# idea. Lets try and make this more modular and flowy. 
#/////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////
import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None  # want to disable this SettingWithCopyWarning

#////////////////////////////////////////////////////////////////////
# mutationsDF_fillIn()
#    goal is to construct a cell-wise dataframe with mutations to each
#    of EGFR, KRAS and BRAF. the challange is getting the cells to line
#    up, hence the for loop 
#
#    GOI needs to be lowercase
#////////////////////////////////////////////////////////////////////
def mutationsDF_fillIn(GOI, GOI_df):
	mutName = GOI + 'Mut'
	for i in range(0,len(mutationsDF.index)):
		currCell = mutationsDF['cell'][i]

		rightIndex = GOI_df['cell'] == currCell
		rightRow = GOI_df[rightIndex]
    
		rightCell = rightRow['cell']
		rightCell = str(rightCell).split()[1]
    
		rightMut = rightRow['mutations']
		rightMut = str(rightMut).split()[1]
    
		mutationsDF[mutName][i] = rightMut

#////////////////////////////////////////////////////////////////////
# removeExtraCharacters_mutationsDF()
#    essentially converting mutationsDF mutation cols from lists to 
#    strings. makes downstream analysis easier
#
#    GOI needs to be lowercase
#////////////////////////////////////////////////////////////////////
def removeExtraCharacters_mutationsDF(GOI):
	mutName = GOI + 'Mut'

	mutationsDF[mutName] = mutationsDF[mutName].str.replace("'", "") # remove quotes
	mutationsDF[mutName] = mutationsDF[mutName].str.replace("[", "") # remove brackets
	mutationsDF[mutName] = mutationsDF[mutName].str.replace("]", "") # remove brackets
	mutationsDF[mutName] = mutationsDF[mutName].str.replace(" ", "") # remove whitespace?

#////////////////////////////////////////////////////////////////////
# genericSummaryTableFillIn()
#    fills in a given (metadata) field in summaryTable. pulls from 
#    patientMetadata (which is a global) and goes cell-by-cell through 
#    summaryTable, filling in fields like patientID/driver_gene
#
#////////////////////////////////////////////////////////////////////
def genericSummaryTableFillIn(metaField, summaryField):
	for i in range(0,len(summaryTable.index)):
		currCell = summaryTable['cell'].iloc[i]
		currPlate = currCell.split('_')[1]
    
		index_to_keep = patientMetadata['plate'] == currPlate
		keepRow = patientMetadata[index_to_keep]
		try:
			currField = list(keepRow[metaField])[0]
			summaryTable[summaryField][i] = currField
		except IndexError:
			continue
			#print('ERROR: plate not found') # these are just the plates were NOT 
        	                                 # including in the analysis

#////////////////////////////////////////////////////////////////////
# fusionsFillIn()
#    Takes the existing fusionsDF (which is just a list of the five fusions
#    we looked for, and what cells they're found in) and populates 
#    summaryTable with this shit
#
#    this works, but holllllyyyy shitttt we can do better
#////////////////////////////////////////////////////////////////////
def fusionsFillIn(fusionsDF_):
	for i in range(0, len(summaryTable.index)):
		currCell = summaryTable['cell'][i]
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
 
		summaryTable['fusions_found'][i] = fusionsListCurr

#////////////////////////////////////////////////////////////////////
# translatedMutsFillIn_EGFR()
#    need to make a 'mutations_found_translated' field that converts our
#    'raw' mutation calls to something that more resembles those reported
#    in our clinical cols. Need a seperate func for EGFR, bc there are 
#    so many potential variants to account for
#
#////////////////////////////////////////////////////////////////////
def translatedMutsFillIn_EGFR():
	for i in range(0,len(summaryTable.index)):
		translatedList = []
		currCell = summaryTable['cell'].iloc[i]
		currMuts_egfr = summaryTable['mutations_found_EGFR'].iloc[i]
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
        
		summaryTable['mutations_found_translated'][i] = translatedList

#////////////////////////////////////////////////////////////////////
# translatedMutsFillIn_nonEGFR()
#    need to make a 'mutations_found_translated' field that converts our
#    'raw' mutation calls to something that more resembles those reported
#    in our clinical cols. This func handles BRAF and KRAS, bc there are
#    only like 2 possible clinically reported muts for them, so we'd might
#    as well keep everything
#
#    want GOI to be capitilized here
#////////////////////////////////////////////////////////////////////
def translatedMutsFillIn_nonEGFR(GOI):
	colName = 'mutations_found_' + GOI
	for i in range(0,len(summaryTable.index)):
		translatedList = []
		currCell = summaryTable['cell'].iloc[i]
		currMuts = summaryTable[colName].iloc[i]
		currMuts_split = currMuts.split(',')
		for item in currMuts_split:
			if item != '' and '?' not in item:
				translatedList.append(GOI + ' ' + item)

		summaryTable['mutations_found_translated'][i] = summaryTable['mutations_found_translated'][i] + translatedList

#////////////////////////////////////////////////////////////////////
# translatedMutsFillIn_fusions()
# 	 need to make a 'mutations_found_translated' field that converts our
#    'raw' mutation calls to something that more resembles those reported
#    in our clinical cols. for fusions this time
#
#////////////////////////////////////////////////////////////////////
def translatedMutsFillIn_fusions():
	for i in range(0,len(summaryTable.index)):
		translatedList = []
		currCell = summaryTable['cell'].iloc[i]
		currFus = summaryTable['fusions_found'].iloc[i]
		currFus_split = currFus.split(',')
		for item in currFus_split:
			if item == 'ALK-EML4':
				translatedList.append('ALK fusion')
				translatedList.append('EML4 fusion')
				translatedList.append('ALK-EML4 fusion')
			elif item != '' and '?' not in item:
				item = item.split('_')[0]
				translatedList.append(item + ' fusion')

		summaryTable['mutations_found_translated'][i] = summaryTable['mutations_found_translated'][i] + translatedList

#////////////////////////////////////////////////////////////////////
# convertToString()
#    really just taking this mutations_found_translated col and converting
#    it from a list to a string. makes taking set() easier, but since
#    this is a script now, maybe i dont even need this 
#
#////////////////////////////////////////////////////////////////////
def convertToString():
	for i in range(0,len(summaryTable.index)):
		currStr = str(summaryTable['mutations_found_translated'][i])
		currStr = currStr.replace("'", "")
		currStr = currStr.replace("]", "")
		currStr = currStr.replace("[", "")
		summaryTable['mutations_found_translated'][i] = currStr

#////////////////////////////////////////////////////////////////////
# clinMutFound_fillIn()
#    want to fill in this clin_mut_found_bool col with 1 if the clinically
#    reported mutation is found, 0 if else
#
#////////////////////////////////////////////////////////////////////
def clinMutFound_fillIn():
	for i in range(0,len(summaryTable.index)):
		currCell = summaryTable['cell'][i]
		currMuts = summaryTable['mutations_found_translated'][i]
		currClinGene = summaryTable['clinical_driver_gene'][i]
		currClinMut = summaryTable['clinical_mutation'][i]
		currClinMut_str = str(currClinGene) + ' ' + str(currClinMut)
    
		if currClinMut_str in currMuts:
			summaryTable['clin_mut_found_bool'][i] = 1
		else:
			summaryTable['clin_mut_found_bool'][i] = 0

#////////////////////////////////////////////////////////////////////
# tumorCellBoolFillIn()
#    want to fill in this tumorCell_bool with 1 if we're calling that
#    cell a tumor cell in our seurat obj, 0 if else
#
#////////////////////////////////////////////////////////////////////
def tumorCellBoolFillIn():
	# NEED TO READ IN SEURAT METADATA, SO WE CAN SET tumorCell_bool
	metaPATH = '/Users/lincoln.harris/Desktop/LAUD_important_shit/metadataSeurat.csv'
	metadataSeurat = pd.read_csv(metaPATH)

	myCols = list(metadataSeurat.columns)
	myCols[0] = 'cell'
	metadataSeurat.columns = myCols
	
	indicies = metadataSeurat['inferCNV_annotation'] == 'perturbed'
	metadataSeurat_pert = metadataSeurat[indicies]
	
	tumorCellsList = list(metadataSeurat_pert['cell'])

	# now fill in 'tumorCell_bool' for summaryTable
	for i in range(0, len(summaryTable.index)):
		currCell = summaryTable['cell'][i]
		if currCell in tumorCellsList:
			summaryTable['tumorCell_bool'][i] = 1
		else:
			summaryTable['tumorCell_bool'][i] = 0

#////////////////////////////////////////////////////////////////////
# getNonZeroCovROI()
#    takes a given coverageByCell dataframe and filters for the non-zero 
#    vals. coverage dfs come from checkCoverage_parallel.py
#
#////////////////////////////////////////////////////////////////////
def getNonZeroCovROI(gene, mut):
	fPATH = '../coverage/out/' + gene + '_' + mut + '_' + 'coverageByCell.csv'
	cov = pd.read_csv(fPATH)
	indices = cov['depth_gvcf'] != 0
	cov_nonZero = cov[indices]

	return(cov_nonZero)

#////////////////////////////////////////////////////////////////////
# ROI_coverage_fillIn()
#    fills in coverage for a given ROI, for summaryTable
#
#////////////////////////////////////////////////////////////////////
def ROI_coverage_fillIn(coverage_df, queryGene, queryMutation):
    for i in range(0, len(summaryTable.index)):
        currCell = summaryTable['cell'][i]
        currDriver = summaryTable['clinical_driver_gene'][i]
        currMut = summaryTable['clinical_mutation'][i]
    
        if currDriver == queryGene and currMut == queryMutation:
            if currCell in list(coverage_df['cellName']):
                index_cov_nonZero = coverage_df['cellName'] == currCell
                currRow_cov_nonZero = coverage_df[index_cov_nonZero]
                currDepth_gvcf = int(currRow_cov_nonZero['depth_gvcf'])
        
                summaryTable['coverage_to_ROI'][i] = currDepth_gvcf
            else:
                summaryTable['coverage_to_ROI'][i] = 0

#////////////////////////////////////////////////////////////////////
# main()
#   reads in our inputs, calls routines, does filtering and finally
#   writes to a file 
#////////////////////////////////////////////////////////////////////
global mutationsDF
global summaryTable
global patientMetadata

print('running...')

# READ IN ALL OF THESE BY-GENE AMINO-ACID LEVEL MUTATION COUNTS OBJECTS
mutsPATH = '/Users/lincoln.harris/code/SNP_calling_pipeline/getMutationCounts/'
egfrPATH = mutsPATH + 'egfr_germline_out_AA.csv'
brafPATH = mutsPATH + 'braf_germline_out_AA.csv'
krasPATH = mutsPATH + 'kras_germline_out_AA.csv'

egfr_df = pd.read_csv(egfrPATH, header=None, names=['cell', 'mutations'])
braf_df = pd.read_csv(brafPATH, header=None, names=['cell', 'mutations'])
kras_df = pd.read_csv(krasPATH, header=None, names=['cell', 'mutations'])

# FIRST STEP IS TO GENERATE THE mutationsDF
print('creating mutationsDF')
mutationsDF = pd.DataFrame(columns=['cell', 'brafMut', 'egfrMut', 'krasMut'])
mutationsDF['cell'] = egfr_df['cell']
mutationsDF['egfrMut'] = egfr_df['mutations'] # fill in EGFR first -- this is ok bc
                                              #  the cell order is based on egfr_df

mutationsDF_fillIn('braf', braf_df) 
mutationsDF_fillIn('kras', kras_df)

removeExtraCharacters_mutationsDF('egfr')
removeExtraCharacters_mutationsDF('braf')
removeExtraCharacters_mutationsDF('kras')

# READ IN patientMetadata
patientMetadata = pd.read_csv('../cDNA_plate_metadata.csv')
patientMetadata = patientMetadata.drop([0,1]) # first two rows are wierd

# INIT THE SUMMARY TABLE
cols = ['cell', 'patient', 'clinical_driver_gene', 'clinical_mutation', 'coverage_to_ROI', 'clin_mut_found_bool', 'mutations_found_EGFR', 'mutations_found_BRAF', 'mutations_found_KRAS', 'fusions_found', 'tumorCell_bool']
summaryTable = pd.DataFrame(columns=cols)
summaryTable['cell'] = mutationsDF['cell']

# FILL IN VARIOUS METADATA COLS
genericSummaryTableFillIn('patient_id', 'patient')
genericSummaryTableFillIn('driver_gene', 'clinical_driver_gene')
genericSummaryTableFillIn('driver_mutation', 'clinical_mutation')

# FILL IN MUTATIONS FOUND COL 
summaryTable['mutations_found_EGFR'] = mutationsDF['egfrMut']
summaryTable['mutations_found_KRAS'] = mutationsDF['krasMut']
summaryTable['mutations_found_BRAF'] = mutationsDF['brafMut']

# READ IN FUSIONS DATAFRAME, THEN FILL IN summaryTable
fusionsDF = pd.read_csv('./fusion_dataframe.csv')
fusionsFillIn(fusionsDF)

# SET UP A COL TO TRANSLATE 'RAW' MUTATION CALLS TO 'CLINICAL'
print('translating mutations')
summaryTable['mutations_found_translated'] = ""
translatedMutsFillIn_EGFR()
translatedMutsFillIn_nonEGFR('KRAS')
translatedMutsFillIn_nonEGFR('BRAF')
translatedMutsFillIn_fusions()

# CONVERT LISTS TO STRING, SO I CAN GET SET -- probably not necessary, actually 
convertToString()

# FILL IN clin_mut_found_bool
clinMutFound_fillIn()

# FILL IN  tumorCellBool
tumorCellBoolFillIn()

# GET PER-CELL ROI COVERAGE DFs
print('getting coverage to ROIs...')
braf_V600E_cov_nonZero = getNonZeroCovROI('braf', 'V600E')
egfr_L858R_cov_nonZero = getNonZeroCovROI('egfr', 'L858R')
egfr_exon19del_cov_nonZero = getNonZeroCovROI('egfr', 'exon19del')
egfr_exon20ins_cov_nonZero = getNonZeroCovROI('egfr', 'exon20ins') # this guy is totally empty...
egfr_G719X_cov_nonZero = getNonZeroCovROI('egfr', 'G719X')
egfr_L861Q_cov_nonZero = getNonZeroCovROI('egfr', 'L861Q')
egfr_S768I_cov_nonZero = getNonZeroCovROI('egfr', 'S768I')
egfr_T790M_cov_nonZero = getNonZeroCovROI('egfr', 'T790M')
kras_G12C_cov_nonZero = getNonZeroCovROI('kras', 'G12C')
kras_G13X_cov_nonZero = getNonZeroCovROI('kras', 'G13X')
kras_Q61X_cov_nonZero = getNonZeroCovROI('kras', 'Q61X')

# FIX UP SOME OF THE WEIRD ONES
kras_G13X_cov_nonZero['depth_gvcf'][4202] = 34	
kras_Q61X_cov_nonZero['depth_gvcf'][6431] = 92
egfr_exon19del_cov_nonZero['depth_gvcf'] = egfr_exon19del_cov_nonZero['depth_gvcf'].str.strip('[')
egfr_exon19del_cov_nonZero['depth_gvcf'] = egfr_exon19del_cov_nonZero['depth_gvcf'].str.strip(']')
egfr_exon19del_cov_nonZero['depth_gvcf'] = egfr_exon19del_cov_nonZero['depth_gvcf'].str.strip("'")

# FILL IN ROI COVERAGE TO SUMMARY TABLE
ROI_coverage_fillIn(braf_V600E_cov_nonZero, 'BRAF', 'V600E')
ROI_coverage_fillIn(egfr_G719X_cov_nonZero, 'EGFR', 'G719X')
ROI_coverage_fillIn(egfr_L858R_cov_nonZero, 'EGFR', 'L858R')
ROI_coverage_fillIn(egfr_L861Q_cov_nonZero, 'EGFR', 'L861Q')
ROI_coverage_fillIn(egfr_S768I_cov_nonZero, 'EGFR', 'S768I')
ROI_coverage_fillIn(egfr_T790M_cov_nonZero, 'EGFR', 'T790M')
ROI_coverage_fillIn(kras_G12C_cov_nonZero, 'KRAS', 'G12C')
ROI_coverage_fillIn(kras_G13X_cov_nonZero, 'KRAS', 'G13X')
ROI_coverage_fillIn(kras_Q61X_cov_nonZero, 'KRAS', 'Q61X')
ROI_coverage_fillIn(egfr_exon19del_cov_nonZero, 'EGFR', 'del19')
ROI_coverage_fillIn(egfr_exon20ins_cov_nonZero, 'EGFR', 'ins20')

# TRIM IT DOWN
print('trimming down...')
summaryTable_trimmed = summaryTable[['cell', 'patient', 'clinical_driver_gene', 'clinical_mutation', 'coverage_to_ROI', 'clin_mut_found_bool', 'tumorCell_bool', 'mutations_found_translated']]
summaryTable_trimmed.columns = ['cell', 'patient', 'clinical_driver_gene', 'clinical_mutation', 'coverage_to_ROI', 'clinical_mutation_found_bool', 'tumorCell_bool', 'mutations_found']
summaryTable_trimmed = summaryTable_trimmed[['cell', 'patient', 'clinical_driver_gene', 'clinical_mutation', 'mutations_found', 'coverage_to_ROI', 'clinical_mutation_found_bool', 'tumorCell_bool']]

# write this fucker
print('writing...')
summaryTable_trimmed.to_csv('/Users/lincoln.harris/Desktop/validationTable_TEST.csv', index=False)

print('done!')

#/////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////