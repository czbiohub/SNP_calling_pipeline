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
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning) # fuck this message

#////////////////////////////////////////////////////////////////////
# getDepth()
#	what does this fucker do? 
#
#////////////////////////////////////////////////////////////////////
def getDepth(df):

	if len(df.index) != 0: 
		print('record found in %s VCF' % cellName)
		infoStr = df['INFO']
		infoStr = str(infoStr)
		DP = infoStr.split('DP')[1].split(';')[0].strip('=')
		print('sequencing depth: %s' % DP)

	else:
		print('no record found in %s VCF' % cellName)

	print(' ')

#////////////////////////////////////////////////////////////////////
# getDepth_g()
#	what does this fucker do? 
#
#////////////////////////////////////////////////////////////////////
def getDepth_g(df):

	if len(df.index) != 0: 
		print('record found in %s gVCF' % cellName)
		infoStr = df['INFO']
		infoStr = str(infoStr)
		DP = infoStr.split('DP')[1].split(';')[0].strip('=')
		print('sequencing depth: %s' % DP)

	else:
		print('no record found in %s gVCF' % cellName)

	print(' ')

#////////////////////////////////////////////////////////////////////
# getGOI()
#	defome a list of records corresponding to the GOI
#
#////////////////////////////////////////////////////////////////////
def getGOI_record(record, *args):
	chrom = 'chr' + str(args[0])
	start = int(args[1])
	end = int(args[2])

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
#			python3 checkCoverage.py 7 55152337 55207337 D12_B003528.vcf D12_B003528.g.vcf 
#////////////////////////////////////////////////////////////////////

global cellName

if len(sys.argv) != 6:
	print('usage: python3 checkCoverage [chrom] [start_pos] [end_pos] [vcf] [gvcf]')
	print('			ie. python3 checkCoverage.py 7 55152337 55207337 D12_B003528.vcf D12_B003528.g.vcf')
	print('  ')
	sys.exit()

print(' ')
print('this tool should be used for loci specific coverage queries.')
print('it is NOT intended for calculating coverage at the exon/transcript level.')

chrom_ = sys.argv[1]
start_ = sys.argv[2]
end_ = sys.argv[3]

vcfFilePrefix = sys.argv[4]
gvcfFilePrefix = sys.argv[5]

cellName = str(vcfFilePrefix).strip('.vcf')

print('  ')
print('chromosome: %s' % chrom_)
print('start_position: %s' % start_)
print('end_position: %s' % end_)
print('cell name: %s' % cellName)
print(' ')

cwd = os.getcwd()
vcf_path = cwd + '/' + vcfFilePrefix
gvcf_path = cwd + '/' + gvcfFilePrefix

vcf = VCF.dataframe(vcf_path)
gvcf = VCF.dataframe(gvcf_path)

toKeepList_v = vcf.apply(getGOI_record, axis=1, args=(chrom_, start_ ,end_))
toKeepList_g = gvcf.apply(getGOI_record, axis=1, args=(chrom_, start_, end_))

vcf_GOI = vcf[np.array(toKeepList_v, dtype=bool)]
gvcf_GOI = gvcf[np.array(toKeepList_g, dtype=bool)]

getDepth(vcf_GOI)
getDepth_g(gvcf_GOI)

#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
