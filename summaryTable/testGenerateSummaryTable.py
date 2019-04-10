# This notebook takes as input GOI_out_AA.csv files (from getMutationCounts_overall_and_GOI.py),
# patient metadata, seurat metadata, fusionsDF, and creates a BY CELL 
# summaryTable. The goal with this table is to provide an answer to questions like
# 'which patients have which mutations?', and 'how many cells have clinically relevant
# mutations?' 
import summarizeLib
import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None  # want to disable this SettingWithCopyWarning

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
mutationsDF['egfrMut'] = egfr_df['mutations'] # fill in EGFR first -- this is ok bc the cell order is based on egfr_df

summarizeLib.mutationsDF_fillIn('braf', braf_df, mutationsDF) 
summarizeLib.mutationsDF_fillIn('kras', kras_df, mutationsDF)

# CONVERTING LISTS INTO STRS. MAKES DOWNSTEAM ANALYSIS EASIER
summarizeLib.removeExtraCharacters_mutationsDF('egfr', mutationsDF)
summarizeLib.removeExtraCharacters_mutationsDF('braf', mutationsDF)
summarizeLib.removeExtraCharacters_mutationsDF('kras', mutationsDF)

# READ IN patientMetadata
patientMetadata = pd.read_csv('/Users/lincoln.harris/code/SNP_calling_pipeline/cDNA_plate_metadata.csv')
patientMetadata = patientMetadata.drop([0,1]) # first two rows are wierd

# INIT THE SUMMARY TABLE
cols = ['cell', 'patient', 'clinical_driver_gene', 'clinical_mutation', 'coverage_to_ROI', 'clin_mut_found_bool', 'mutations_found_EGFR', 'mutations_found_BRAF', 'mutations_found_KRAS', 'fusions_found', 'tumorCell_bool']
summaryTable = pd.DataFrame(columns=cols)
summaryTable['cell'] = mutationsDF['cell']

# FILL IN VARIOUS METADATA COLS
summarizeLib.genericSummaryTableFillIn('patient_id', 'patient', summaryTable, patientMetadata)
summarizeLib.genericSummaryTableFillIn('driver_gene', 'clinical_driver_gene', summaryTable, patientMetadata)
summarizeLib.genericSummaryTableFillIn('driver_mutation', 'clinical_mutation', summaryTable, patientMetadata)

# FILL IN MUTATIONS FOUND COL 
summaryTable['mutations_found_EGFR'] = mutationsDF['egfrMut']
summaryTable['mutations_found_KRAS'] = mutationsDF['krasMut']
summaryTable['mutations_found_BRAF'] = mutationsDF['brafMut']

# READ IN FUSIONS DATAFRAME, THEN FILL IN summaryTable
fusionsDF = pd.read_csv('/Users/lincoln.harris/code/SNP_calling_pipeline/summaryTable/fusion_dataframe.csv')

summarizeLib.fusionsFillIn(fusionsDF, summaryTable)

# SET UP A COL TO TRANSLATE 'RAW' MUTATION CALLS TO 'CLINICAL'
summaryTable['mutations_found_translated'] = ""
summarizeLib.translatedMutsFillIn_EGFR(summaryTable)
summarizeLib.translatedMutsFillIn_nonEGFR('KRAS', summaryTable)
summarizeLib.translatedMutsFillIn_nonEGFR('BRAF', summaryTable)
summarizeLib.translatedMutsFillIn_fusions(summaryTable)

# CONVERT LISTS TO STRING, SO I CAN GET SET -- probably not necessary, actually 
summarizeLib.convertToString(summaryTable)

# FILL IN clin_mut_found_bool
summarizeLib.clinMutFound_fillIn(summaryTable)

# FILL IN  tumorCellBool
summarizeLib.tumorCellBoolFillIn(summaryTable)

# GET PER-CELL ROI COVERAGE DFs
braf_V600E_cov_nonZero = summarizeLib.getNonZeroCovROI('braf', 'V600E')
egfr_L858R_cov_nonZero = summarizeLib.getNonZeroCovROI('egfr', 'L858R')
egfr_exon19del_cov_nonZero = summarizeLib.getNonZeroCovROI('egfr', 'exon19del')
egfr_exon20ins_cov_nonZero = summarizeLib.getNonZeroCovROI('egfr', 'exon20ins') # this guy is totally empty...
egfr_G719X_cov_nonZero = summarizeLib.getNonZeroCovROI('egfr', 'G719X')
egfr_L861Q_cov_nonZero = summarizeLib.getNonZeroCovROI('egfr', 'L861Q')
egfr_S768I_cov_nonZero = summarizeLib.getNonZeroCovROI('egfr', 'S768I')
egfr_T790M_cov_nonZero = summarizeLib.getNonZeroCovROI('egfr', 'T790M')
kras_G12C_cov_nonZero = summarizeLib.getNonZeroCovROI('kras', 'G12C')
kras_G13X_cov_nonZero = summarizeLib.getNonZeroCovROI('kras', 'G13X')
kras_Q61X_cov_nonZero = summarizeLib.getNonZeroCovROI('kras', 'Q61X')

# FIX UP SOME OF THE WEIRD ONES
kras_G13X_cov_nonZero['depth_gvcf'][4202] = 34  
kras_Q61X_cov_nonZero['depth_gvcf'][6431] = 92
egfr_exon19del_cov_nonZero['depth_gvcf'] = egfr_exon19del_cov_nonZero['depth_gvcf'].str.strip('[')
egfr_exon19del_cov_nonZero['depth_gvcf'] = egfr_exon19del_cov_nonZero['depth_gvcf'].str.strip(']')
egfr_exon19del_cov_nonZero['depth_gvcf'] = egfr_exon19del_cov_nonZero['depth_gvcf'].str.strip("'")

# FILL IN ROI COVERAGE TO SUMMARY TABLE
summarizeLib.ROI_coverage_fillIn(braf_V600E_cov_nonZero, 'BRAF', 'V600E', summaryTable)
summarizeLib.ROI_coverage_fillIn(egfr_G719X_cov_nonZero, 'EGFR', 'G719X', summaryTable)
summarizeLib.ROI_coverage_fillIn(egfr_L858R_cov_nonZero, 'EGFR', 'L858R', summaryTable)
summarizeLib.ROI_coverage_fillIn(egfr_L861Q_cov_nonZero, 'EGFR', 'L861Q', summaryTable)
summarizeLib.ROI_coverage_fillIn(egfr_S768I_cov_nonZero, 'EGFR', 'S768I', summaryTable)
summarizeLib.ROI_coverage_fillIn(egfr_T790M_cov_nonZero, 'EGFR', 'T790M', summaryTable)
summarizeLib.ROI_coverage_fillIn(kras_G12C_cov_nonZero, 'KRAS', 'G12C', summaryTable)
summarizeLib.ROI_coverage_fillIn(kras_G13X_cov_nonZero, 'KRAS', 'G13X', summaryTable)
summarizeLib.ROI_coverage_fillIn(kras_Q61X_cov_nonZero, 'KRAS', 'Q61X', summaryTable)
summarizeLib.ROI_coverage_fillIn(egfr_exon19del_cov_nonZero, 'EGFR', 'del19', summaryTable)
summarizeLib.ROI_coverage_fillIn(egfr_exon20ins_cov_nonZero, 'EGFR', 'ins20', summaryTable)

# TRIM IT DOWN
summaryTable_trimmed = summaryTable[['cell', 'patient', 'clinical_driver_gene', 'clinical_mutation', 'coverage_to_ROI', 'clin_mut_found_bool', 'tumorCell_bool', 'mutations_found_translated']]
summaryTable_trimmed.columns = ['cell', 'patient', 'clinical_driver_gene', 'clinical_mutation', 'coverage_to_ROI', 'clinical_mutation_found_bool', 'tumorCell_bool', 'mutations_found']
summaryTable_trimmed = summaryTable_trimmed[['cell', 'patient', 'clinical_driver_gene', 'clinical_mutation', 'mutations_found', 'coverage_to_ROI', 'clinical_mutation_found_bool', 'tumorCell_bool']]

# WRITE TO FILE
summaryTable_trimmed.to_csv('/Users/lincoln.harris/Desktop/validationTable_TEST.csv', index=False)

# ADD SAMPLE_NAME COL TO SUMMARYTABLE
summaryTable_trimmed['sample_name'] = ''
summarizeLib.genericSummaryTableFillIn('sample_name', 'sample_name', summaryTable_trimmed, patientMetadata)

# GET MIN SET OF SAMPLE NAMES
relevantSamplesSet = set(summaryTable_trimmed['sample_name'])
relevantSamplesList = list(relevantSamplesSet)
relevantSamplesSeries = pd.Series(relevantSamplesList)

# INIT VALIDATIONTABLE_SAMPLES
cols = ['sample', 'patient', 'driver_gene', 'driver_mutation', 'mutations_found', 'numCells', 'numTumorCells', 'numTumorCells_w_coverage_to_ROI', 'numTumorCells_clinMut_found']
validationTable_samples = pd.DataFrame(columns=cols)
validationTable_samples['sample'] = relevantSamplesSeries

# FILL IN METADATA FIELDS
summarizeLib.validationTable_metadata_fillIn('patient_id', 'patient', validationTable_samples, patientMetadata)
summarizeLib.validationTable_metadata_fillIn('driver_gene', 'driver_gene', validationTable_samples, patientMetadata)
summarizeLib.validationTable_metadata_fillIn('driver_mutation', 'driver_mutation', validationTable_samples, patientMetadata)

# FILL IN MUTATIONS FOUND
muts_dict = summarizeLib.validationTable_dict_muts(validationTable_samples, summaryTable_trimmed)
validationTable_samples['mutations_found'] = muts_dict.values()

# FILL IN NUMTUMORCELLS
tc_dict = summarizeLib.validationTable_dict_generic(validationTable_samples, summaryTable_trimmed, 'tumorCell_bool')
tc_cov_dict = summarizeLib.validationTable_dict_generic(validationTable_samples, summaryTable_trimmed, 'coverage_to_ROI')
clinMut_dict = summarizeLib.validationTable_dict_generic(validationTable_samples, summaryTable_trimmed, 'clinical_mutation_found_bool')

validationTable_samples['numTumorCells'] = tc_dict.values()
validationTable_samples['numTumorCells_w_coverage_to_ROI'] = tc_cov_dict.values()
validationTable_samples['numTumorCells_clinMut_found'] = clinMut_dict.values()

# CLEAN UP A BIT
validationTable_samples = validationTable_samples.drop([26])
cols = ['sample', 'patient', 'driver_gene', 'driver_mutation', 'mutations_found', 'numTumorCells', 'numTumorCells_w_coverage_to_ROI', 'numTumorCells_clinMut_found']
validationTable_samples = validationTable_samples[cols]

# write this bitch
validationTable_samples.to_csv('./validationTable_samples_TEST.csv', index=False)

