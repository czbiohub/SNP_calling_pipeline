""" creates a cell-wise dataframe for variant reads normalized by total 
	reads to a given loci. starting point is raw vcfs and a genomic 
	coordinate batch file that delineates every ROI. 

	parallel whooo, whoo! """ 
from collections import OrderedDict
import gzip
import pandas as pd
import os
import multiprocessing as mp
import warnings
import numpy as np
warnings.simplefilter(action='ignore', category=FutureWarning)


def count_comments(filename):
	""" Count comment lines (those that start with "#") 
		cribbed from slowko """
	comments = 0
	fn_open = gzip.open if filename.endswith('.gz') else open
	with fn_open(filename) as fh:
		for line in fh:
			if line.startswith('#'):
				comments += 1
			else:
				break
	return comments



def vcf_to_dataframe(filename):
	""" Open a VCF file and return a pandas.DataFrame with
		each INFO field included as a column in the dataframe 
		cribbed from slowko """
	VCF_HEADER = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', '20']

	# Count how many comment lines should be skipped.
	comments = count_comments(filename)
	tbl = pd.read_table(filename, compression=None, skiprows=comments,
							names=VCF_HEADER, usecols=range(10))
	
	return(tbl)



def ROI_df_subset(vcf, chrom, start, end):
	""" subset a single vcf based on genomic coords"""
	chrStr = 'chr' + str(chrom)
    
	keep0 = vcf['CHROM'] == chrStr
	vcf_sub0 = vcf[keep0]

	keep1 = vcf_sub0['POS'] >= start
	vcf_sub1 = vcf_sub0[keep1]

	keep2 = vcf_sub1['POS'] <= end
	vcf_sub2 = vcf_sub1[keep2]
    
	return(vcf_sub2)



def coverage_search(df):
	""" given single vcf entry, search for AD (DepthPerAlleleBySample) col """ 
	for i in range(0, len(df.index)):
		row = df.iloc[i]
		extra_col = row['20']
		AD = extra_col.split(':')[1]

		wt_count = int(AD.split(',')[0])
		variant_count = int(AD.split(',')[1])
		total_count = wt_count + variant_count

		ret = str(variant_count) + ':' + str(total_count)
		return(ret)



def driver(cellFile):
	""" search for given ROI, across all cells """
	currCell = cellFile.split('/')[7]
	currCell = currCell.strip('.vcf')
	currCell = currCell.strip('.') # why do cell names retain this period? 
	print(currCell)

	vcf_ = vcf_to_dataframe(cellFile)
	ROI_counts = []

	for i in range(0,len(ROI_df.index)):
		row = ROI_df.iloc[i]
		start_ = int(row['start_pos'])
		chrom_ = int(row['chrom'])
		end_ = int(row['end_pos'])
		ROI_ = row['outfile']

		vcf_sub = ROI_df_subset(vcf_, chrom_, start_, end_)
    
		if not vcf_sub.empty:
			counts = coverage_search(vcf_sub)
		else:
			counts = '0:0'

		ROI_counts.append(counts)

	tup = [currCell, ROI_counts]

	return(tup)



""" main. search for each ROI. """
global currPATH
global ROI_df
num_proc = 16

currPATH = os.getcwd()
vcf_list = os.listdir('vcf/')
cell_files_list = [currPATH + '/vcf/' + s for s in vcf_list] # list comprehension

cells_list = []
for item in cell_files_list: # just want the cell names
	currCell = item.split('/')[7]
	currCell = currCell.strip('.vcf')
	cells_list.append(currCell)

ROI_df = pd.read_csv('batch_files/coverageBatch_comb.csv')
ROIs = list(ROI_df['outfile']) # define columns list

print('creating pool')
p = mp.Pool(processes=num_proc)
print('running...')

try:
	outList = p.map(driver, cell_files_list, chunksize=1) # default chunksize=1
finally:
	p.close()
	p.join()

print('writing file')
# convert to dictionary
cells_dict = {}
naSeries = pd.Series([np.nan])

# add each item to dict
for item in outList:
	cell = item[0]
	cov = item[1]
		
	toAdd = {cell:cov}
	cells_dict.update(toAdd)

t = pd.DataFrame.from_dict(cells_dict, orient="index") # orient refers to row/col orientation 

t.columns = ROIs
t.insert(0, 'cells', cells_list)
t.to_csv('big_test.csv', index=False)
