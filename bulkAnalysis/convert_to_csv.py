#/////////////////////////////////////////////////////////////////////////
# script: convert_to_csv.py
# author: Lincoln
# date: 3.18.19
#
# want to convert any remaining vcfs to csv
#/////////////////////////////////////////////////////////////////////////
import pandas as pd
import VCF
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

filterDir = '/home/ubuntu/code/SNP_calling_pipeline/bulkAnalysis/scVCF_filtered_all/'
filterDir_list = os.listdir(filterDir)

for f in filterDir_list:
	if '.vcf' in f:
		currPATH = filterDir + f
		df = VCF.dataframe(currPATH)
		df_trimmed = df[['CHROM', 'POS', 'ID', 'REF', 'ALT']]

		cellName = f.strip('.vcf')
		outStr = filterDir + cellName + '.csv'
		df_trimmed.to_csv(outStr, index=False)

#/////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////
