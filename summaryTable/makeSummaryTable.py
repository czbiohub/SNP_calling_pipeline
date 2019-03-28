#/////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////
# script: makeSummaryTable.py
# author: Lincoln
# date: 3.28.19
#
# TODO: add top level description here!!
# 
# the idea is to convert my existing ipynb to a script, bc its gotten super
# 																	unwieldy
# lets try and make this more modular and flowy
#/////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////
import pandas as pd
import numpy as np

#////////////////////////////////////////////////////////////////////
# mutationsDF_fillIn()
#    TODO: what does this function do? 
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
    # think i can do this without returning anything

#////////////////////////////////////////////////////////////////////
# removeExtraCharacters_mutationsDF()
#    TODO: what does this function do? 
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
#    TODO: what does this function do? 
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
#    TODO: what does this function do? 
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
#    TODO: what does this function do? 
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
#    TODO: what does this function do? 
#
#////////////////////////////////////////////////////////////////////
def translatedMutsFillIn_nonEGFR(GOI):
	for i in range(0,len(summaryTable.index)):
    	translatedList = []
		currCell = summaryTable['cell'].iloc[i]
		currMuts = summaryTable['mutations_found_BRAF'].iloc[i]
		currMuts_split = currMuts.split(',')
		for item in currMuts_split:
			if item != '' and '?' not in item:
				translatedList.append('BRAF ' + item)

		summaryTable['mutations_found_translated'][i] = summaryTable['mutations_found_translated'][i] + translatedList

#////////////////////////////////////////////////////////////////////
# main()
#   not sure how we're gonna structure this, yet
#   but pretty sure its gonna be super long
#////////////////////////////////////////////////////////////////////
global mutationsDF
global summaryTable

# READ IN ALL OF THESE BY-GENE AMINO-ACID LEVEL MUTATION COUNTS OBJECTS
mutsPATH = '/Users/lincoln.harris/code/SNP_calling_pipeline/getMutationCounts/'
egfrPATH = mutsPATH + 'egfr_germline_out_AA.csv'
brafPATH = mutsPATH + 'braf_germline_out_AA.csv'
krasPATH = mutsPATH + 'kras_germline_out_AA.csv'

egfr_df = pd.read_csv(egfrPATH, header=None, names=['cell', 'mutations'])
braf_df = pd.read_csv(brafPATH, header=None, names=['cell', 'mutations'])
kras_df = pd.read_csv(krasPATH, header=None, names=['cell', 'mutations'])

# FIRST STEP IS TO GENERATE THE mutationsDF
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
summaryTable['mutations_found_translated'] = ""
translatedMutsFillIn_EGFR()
translatedMutsFillIn_nonEGFR('kras')
translatedMutsFillIn_nonEGFR('braf')

#/////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////