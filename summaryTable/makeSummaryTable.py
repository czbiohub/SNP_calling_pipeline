#/////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////
# script: makeSummaryTable.py
# author: Lincoln
# date: 3.28.19
#
# TODO: add top level description here!!
# 
# the idea is to convert my existing ipynb to a script, bc its gotten super
# unwieldy 
#/////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////
import pandas as pd
import numpy as np
   
#////////////////////////////////////////////////////////////////////
# fill_in_mutationsDF()
#    TODO: what does this function do? 
#    GOI needs to be lowercase
#////////////////////////////////////////////////////////////////////
def fill_in_mutationsDF(GOI, GOI_df):
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
# main()
#   not sure how we're gonna structure this, yet
#////////////////////////////////////////////////////////////////////
global mutationsDF

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

fill_in_mutationsDF('braf', braf_df) 
fill_in_mutationsDF('kras', kras_df)

#/////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////