""" TODO add description"""
import pandas as pd


def driver():
	""" loops through variant_df, matches cells to samples, and fills in 
		samples_x_gene with read count values """

	for i in range(0,len(variant_df.index)):# looping over by-cell df
		currCell = variant_df['cells'].iloc[i]
		keep = meta['cell_id'] == currCell
		meta_row = meta[keep]
	    
		try:
			currSample = list(meta_row['sample_name'])[0]
		except IndexError: # cells not in metadata? 
			currSample = 'NOT_FOUND'
	        
		for j in range(1,len(variant_df.columns)): # starting at 1 bc i dont want the cellnames
			ROI = variant_df.columns[j]
			currGene = ROI.split('_')[0]
	        
			samples_x_gene_sub = samples_x_gene.where(samples_x_gene['gene'] == currGene) # bottleneck
			gene_index = samples_x_gene_sub.index[pd.notna(samples_x_gene_sub['gene'])].tolist()
			gene_index = gene_index[0]
	        
			currVal = variant_df.iloc[i,j]
	        
			samples_x_gene[currSample][gene_index] += currVal



def main():
	""" read in variants-by-cell df as well as metadata, and set up samples-x-gene df """
	global variant_df
	global meta
	global samples_x_gene

	# read in
	variant_df = pd.read_csv('variant_cov_df.csv')
	meta = pd.read_csv('../metadata_all_cells_4.10.19.csv')

	# get set of sample names
	samples_set = set(meta['sample_name'])
	samples_set = list(samples_set)

	# get df containing all ROIs
	ROI_series = pd.Series(variant_df.columns).drop([0])
	ROI_df = pd.DataFrame(ROI_series, columns=['ROI'])
	ROI_df['gene'] = ROI_df['ROI'].str.split('_').str[0]

	# get unique genes, from ROI df
	all_genes = list(ROI_df['gene'])
	unique_genes = set(all_genes)
	unique_genes = list(unique_genes)

	# add 'gene' as the first elem in samples_set
	samples_set.insert(0,'gene')

	# i think this is the df we want to work on
	samples_x_gene = pd.DataFrame(columns=samples_set)
	samples_x_gene['gene'] = unique_genes
	samples_x_gene

	# set everybody to 0
	for col in samples_x_gene.columns:
		if col != 'gene':
			samples_x_gene[col].values[:] = 0

	driver()
	samples_x_gene.to_csv('samples_x_gene_variant.csv', index=False)



if __name__== "__main__": # call main
	main()



