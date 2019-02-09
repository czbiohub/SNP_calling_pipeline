#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
# script: checkCoverage.py
# author: Lincoln 
# date: 2.8.18
#
# ooOOOOohhHHHHH baby here we go!!!
# Want to turn my jupyter notebook into a python script
#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
import pandas as pd
import numpy as np
import VCF
import sys
import os

#////////////////////////////////////////////////////////////////////
# getGOI()
#	defome a list of records corresponding to the GOI
#
#////////////////////////////////////////////////////////////////////
def getGOI_records(record, *args):
	chrom = args[0]
	start = args[1]
	end = args[2]

	if record['CHROM'] == chrom:
		if end >= record['POS'] >= start:
			return 1
		else:
			return 0
	else:
		return 0

#////////////////////////////////////////////////////////////////////
# main()
#	Main logic here. Prints usage, parses cmd line args, calls 
#		appropriate driver function 
#		
#		EGFR test: 
#			python3 checkCoverage.py 7 55152337 55207337 egfr_out.csv
#////////////////////////////////////////////////////////////////////

if len(sys.argv) == 1:
	print('usage: python3 checkCoverage [chrom] [start_pos] [end_pos] [vcf] [gvcf]')
	print('			ie. python3 checkCoverage.py 4 7 500050 50010 exampleCell.vcf exampleCell.g.vcf')
	print('  ')
	sys.exit()

chrom_ = sys.argv[2]
start_ = sys.argv[3]
end_ = sys.argv[4]

vcfFilePrefix = sys.argv[5]
gvcfFilePrefix = sys.argv[6]

cwd = os.getcwd()
vcf_path = cwd + '/' + vcfFilePrefix
gvcf_path = cwd + '/' + gvcfFilePrefix

vcf = VCF.dataframe(vcf_path)
gvcf = VCF.dataframe(gvcf_path)

toKeepList_v = vcf.apply(getGOI_record, axis=1, args=(chrom_, start_ ,end_))
toKeepList_g = gvcf.apply(getGOI_record, axis=1, args=(chrom_, start_, end_))

vcf_GOI = df[np.array(toKeepList_v, dtype=bool)]
gvcf_GOI = gvcf[np.array(toKeepList_g, dtype=bool)]

#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
