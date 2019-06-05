""" TODO add description""" 
from collections import OrderedDict
import gzip
import pandas as pd
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


def count_comments(filename):
	""" Count comment lines (those that start with "#") in an optionally
	gzipped file """
	comments = 0
	fn_open = gzip.open if filename.endswith('.gz') else open
	with fn_open(filename) as fh:
		for line in fh:
			if line.startswith('#'):
				comments += 1
			else:
				break
	return comments



def vcf_to_dataframe(filename, large=True):
	"""Open a VCF file and return a pandas.DataFrame with
	each INFO field included as a column in the dataframe.

	:param filename:    An optionally gzipped VCF file.
	:param large:       Use this with large VCF files to skip the ## lines and
                        leave the INFO fields unseparated as a single column.
    """
	VCF_HEADER = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', '20']

	if large:
		# Set the proper argument if the file is compressed.
		comp = 'gzip' if filename.endswith('.gz') else None
		# Count how many comment lines should be skipped.
		comments = count_comments(filename)
		# Return a simple DataFrame without splitting the INFO column.
		return pd.read_table(filename, compression=comp, skiprows=comments,
							names=VCF_HEADER, usecols=range(10))

	# Each column is a list stored as a value in this dict. The keys for this
	# dict are the VCF column names and the keys in the INFO column.
	result = OrderedDict()
	# Parse each line in the VCF file into a dict.
	for i, line in enumerate(lines(filename)):
		for key in line.keys():
			# This key has not been seen yet, so set it to None for all
			# previous lines.
			if key not in result:
				result[key] = [None] * i
			# Ensure this row has some value for each column.
		for key in result.keys():
			result[key].append(line.get(key, None))

	return pd.DataFrame(result)



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
		print(AD)



def ROI_search(row):
	""" search for given ROI, across all cells """
	chrom_ = row['chrom']
	start_ = row['start_pos']
	end_ = row['end_pos']
    
	for f in cell_files_list:
		vcf_ = vcf_to_dataframe(f)
		vcf_sub = ROI_df_subset(vcf_, chrom_, start_, end_)
    
		if not vcf_sub.empty:
			coverage_search(vcf_sub)


""" main. search for each ROI. """
global cell_files_list

currPATH = os.getcwd()
cell_files_list = os.listdir('vcf/')
cell_files_list = [currPATH + '/vcf/' + s for s in cell_files_list] # list comprehension

ROI_df = pd.read_csv('/Users/lincoln.harris/code/cerebra/py_notebooks/coverageBatch_v1.csv')
ROI_df.apply(ROI_search, axis=1)


