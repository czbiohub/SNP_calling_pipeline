#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
# script: create_runbatch_config.py
# author: Lincoln 
# date: 12.13.18
#
# Seriously how dont i have one of these already? 
# run with 
#		ipython create_runbatch_config.py
#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
import os
import json
import pandas as pd
pd.options.display.max_colwidth = 500
pd.options.mode.chained_assignment = None 

#////////////////////////////////////////////////////////////////////
# writeFunc()
#	Write both samples.csv and config output files
#
#////////////////////////////////////////////////////////////////////
def writeFunc(samples_df):
	# where do you want to write it? 
	out_dir = '../gatk/12.10_run'
	# write samples_df to file
	get_ipython().system(' mkdir -p $out_dir')
	samples_df.to_csv(f'{out_dir}/samples.csv', index=False)
	# write a config file
	config =     {
		"program": "../../reflow/gatk_pipeline.rf",
		"runs_file": "samples.csv"
	}

#////////////////////////////////////////////////////////////////////
# get_bam()
#	what does this fucker do? 
#////////////////////////////////////////////////////////////////////

def get_bam(cell):
    s3_location = f'{prefix}{cell}' #f? 
    lines = get_ipython().getoutput('aws s3 ls $s3_location') 
    
    try:
    	bam_line = [x for x in lines if x.endswith('bam')][0] # get the bam file, specifically
    	bam_basename = bam_line.split()[-1]
    except IndexError:
    	return('dummy') # think this will work so long as i return something

    return f'{s3_location}{bam_basename}'

#////////////////////////////////////////////////////////////////////
# driver()
#	what does this fucker do? 
#////////////////////////////////////////////////////////////////////

def driver(prefix): 
    txt = 'runX_cells.txt'
    get_ipython().getoutput('aws s3 ls $prefix > $txt')
    
    # read into a pandas dataframe
    cells_df = pd.read_table(txt, delim_whitespace=True, header=None, names=['is_prefix', 'cell_name'])
    
    cells_df['input_bam'] = cells_df['cell_name'].map(get_bam) # call get_bam() and add 'input_bam' col
    cells_df['sample_id'] = cells_df.cell_name.str.strip('/') # get rid of forward slashes and add 'sample_id' col
    cells_df['id'] = cells_df['sample_id'] # add an ID col
    
    # add output_prefix col
    cells_df['output_vcf'] = 's3://darmanis-group/singlecell_lungadeno/non_immune/nonImmune_bams_9.27/vcf1/' + cells_df['sample_id'] + '.vcf'               
    
    # subset cells_df by only what we want
    cols_to_keep = ['id', 'input_bam', 'sample_id', 'output_vcf']
    samples_df = cells_df[cols_to_keep]
    
    return(samples_df)

#////////////////////////////////////////////////////////////////////
# main()
#	what does this fucker do? 
#////////////////////////////////////////////////////////////////////

bucketPrefixes = 's3://darmanis-group/singlecell_lungadeno/non_immune/nonImmune_bams_9.27/'
f = 'myCells.txt'
get_ipython().system(' aws s3 ls $bucketPrefixes > $f')
    
# read run prefixes into a pandas df
runs_df = pd.read_table(f, delim_whitespace=True, header=None, names=['is_prefix', 'run_name'])
    
# add a full_path col
runs_df['full_path'] = 's3://darmanis-group/singlecell_lungadeno/non_immune/nonImmune_bams_9.27/' + runs_df['run_name']
    

big_df = pd.DataFrame() # init empty dataframe

for i in range(0, len(runs_df.index)-2): # -2 bc i have two vcf out folders
	global prefix # dont like this
	prefix = runs_df['full_path'][i]
	print(prefix)
	curr_df = driver(prefix)
	toConcat = [big_df, curr_df]
	big_df = pd.concat(toConcat)
	print(big_df.shape)
	writeFunc(big_df) # bc im nervous as FUCK 

writeFunc(big_df)

#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
