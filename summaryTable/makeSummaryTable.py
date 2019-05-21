# This notebook takes as input GOI_out_AA.csv files (from getMutationCounts_overall_and_GOI.py),
# patient metadata, seurat metadata, fusionsDF, and creates both by CELL and 
# by SAMPLE summaryTables. The goal with this table is to provide an answer to questions like
# 'which patients have which mutations?', and 'how many cells have clinically relevant
# mutations?' 

import summarizeModule
import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None  # want to disable this SettingWithCopyWarning

# FIRST STEP IS TO GENERATE THE mutationsDF
mutationsDF = pd.DataFrame(columns=['cell', 'AKT1_mut', 'ALK_mut', 'BAP1_mut', 
                                    'BRAF_mut', 'DDR2_mut', 'DROSHA_mut', 'EGFR_mut',
                                    'ERBB2_mut', 'ERBB4_mut', 'FGFR2_mut', 'GRIN2A_mut',
                                    'HIF1A_mut', 'KDR_mut', 'KEAP1_mut', 'KRAS_mut', 'MAP2K1_mut', 
                                    'MAP2K2_mut', 'MYCL_mut', 'NFE2L2_mut', 'NKX21_mut', 'NOTCH1_mut',
                                    'PIK3CB_mut', 'PTPN13_mut', 'PTPRT_mut', 'RAD21_mut', 'RB1_mut',
                                    'RBM10_mut', 'RET_mut', 'SMARCA4_mut', 'SOX2_mut', 'STK11_mut',
                                    'TP53_mut', 'TP63_mut',])

muts_path = '/Users/lincoln.harris/code/SNP_calling_pipeline/getMutationCounts/out/'
genesList = ['AKT1', 'ALK', 'BAP1', 'BRAF', 'DDR2', 'DROSHA', 'EGFR', 'ERBB2', 'ERBB4', 'FGFR2', 'GRIN2A',
             'HIF1A', 'KDR', 'KEAP1', 'KRAS', 'MAP2K1', 'MAP2K2', 'MYCL', 'NFE2L2', 'NKX21', 'NOTCH1',
             'PIK3CB', 'PTPN13', 'PTPRT', 'RAD21', 'RB1', 'RBM10', 'RET', 'SMARCA4', 'SOX2', 'STK11', 'TP53',
             'TP63',]

EGFR_path = muts_path + 'EGFR_out_AA.csv'
EGFR_df = pd.read_csv(EGFR_path, header=None, names=['cell', 'mutations'])
mutationsDF['cell'] = EGFR_df['cell']
mutationsDF['EGFR_mut'] = EGFR_df['mutations'] # fill in EGFR first -- this is ok bc the cell order is based on egfr_df

for gene in genesList:
    gene_path = muts_path + gene + '_out_AA.csv'
    gene_df = pd.read_csv(gene_path, header=None, names=['cell', 'mutations'])
    summarizeModule.mutationsDF_fillIn(gene, gene_df, mutationsDF)
    summarizeModule.removeExtraCharacters_mutationsDF(gene, mutationsDF)

# READ IN patientMetadata
patientMetadata = pd.read_csv('/Users/lincoln.harris/code/SNP_calling_pipeline/metadata_all_cells_4.10.19.csv')
patientMetadata = patientMetadata.drop([0,1]) # first two rows are wierd
patientMetadata

# INIT THE SUMMARY TABLE
cols = ['cell', 'patient', 'clinical_driver_gene', 'clinical_mutation', 'coverage_to_ROI', 
             'clin_mut_found_bool', 'mutations_found_AKT1', 'mutations_found_ALK',
             'mutations_found_BAP1', 'mutations_found_BRAF', 'mutations_found_DDR2', 'mutations_found_DROSHA', 
             'mutations_found_EGFR', 'mutations_found_ERBB2','mutations_found_ERBB4', 'mutations_found_FGFR2',
             'mutations_found_GRIN2A','mutations_found_HIF1A', 
             'mutations_found_KDR', 'mutations_found_KEAP1', 'mutations_found_KRAS', 'mutations_found_MAP2K1', 
             'mutations_found_MAP2K2', 'mutations_found_MYCL', 'mutations_found_NFE2L2', 'mutations_found_NKX21',
             'mutations_found_NOTCH1','mutations_found_PIK3CB', 'mutations_found_PTPN13', 'mutations_found_PTPRT',
             'mutations_found_RAD21', 'mutations_found_RB1', 'mutations_found_RBM10', 'mutations_found_RET',
             'mutations_found_SMARCA4', 'mutations_found_SOX2', 'mutations_found_STK11', 'mutations_found_TP53',
             'mutations_found_TP63', 'fusions_found', 'tumorCell_bool']

summaryTable = pd.DataFrame(columns=cols)
summaryTable['cell'] = mutationsDF['cell']

# FILL IN VARIOUS METADATA COLS
summarizeModule.genericSummaryTableFillIn('patient_id', 'patient', summaryTable, patientMetadata)
summarizeModule.genericSummaryTableFillIn('driver_gene', 'clinical_driver_gene', summaryTable, patientMetadata)
summarizeModule.genericSummaryTableFillIn('driver_mutation', 'clinical_mutation', summaryTable, patientMetadata)

# FILL IN MUTATIONS FOUND COL 
for gene in genesList:
    summary_col = 'mutations_found_' + gene
    mutsdf_col = gene + '_mut'
    summaryTable[summary_col] = mutationsDF[mutsdf_col]

# READ IN FUSIONS DATAFRAME, THEN FILL IN summaryTable
fusionsDF = pd.read_csv('/Users/lincoln.harris/code/SNP_calling_pipeline/summaryTable/fusion_dataframe.csv')

summarizeModule.fusionsFillIn(fusionsDF, summaryTable)