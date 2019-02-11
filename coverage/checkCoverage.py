#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
# script: checkCoverage.py
# author: Lincoln 
# date: 2.8.18
#
# ooOOOOohhHHHHH baby here we go!!!
#
# this tool compares VCF records to gVCF records, for a given cell. 
# input a genomic loci of interest, and both VCF and gVCF for a given
# cell, and it will spit out whether that loci exists in either file, 
# and the depth of coverage if found. 
#
# im imagining this being used in a low-throughput manner, for loci
# of interest within individual cells. not sure how to scale it up,
# or if thats even possible. 
#
# also keep in mind - this works MUCH BETTER for small loci, ie. 
# individual SNPs / small indels. NOT INTENDED for whole 
# exon or whole transcript queries. 
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
# getDepth_adv()
#	more advanced version of the function(s) below; can give it a 
#	dataframe containing multiple records, and depth will be reported
#	for every record within that df.
#		for VCF file
#////////////////////////////////////////////////////////////////////
def getDepth_adv(df):

	if len(df.index) == 0:
		print('no record in %s VCF' % cellName)
	elif len(df.index) == 1:
		print('record found in %s VCF' % cellName)
		infoStr = df['INFO']
		infoStr = str(infoStr)
		DP = infoStr.split('DP')[1].split(';')[0].strip('=')
		print('sequencing depth: %s' % DP)
	else:
		print('multiple records found in %s VCF' % cellName)
		infoDF = df['INFO']

		for i in range(0, len(infoDF.index)-1):
			line = infoDF.iloc[i]
			line = str(line)
			try:
				DP = line.split('DP')[1].split(';')[0].strip('=')
				print('       sequencing depth (record %d): %s' % (i, DP))
			except IndexError:
				continue

	print(' ')

#////////////////////////////////////////////////////////////////////
# getDepth_adv_g()
#	more advanced version of the function(s) below; can give it a 
#	dataframe containing multiple records, and depth will be reported
#	for every record within that df.
#		for gVCF file
#////////////////////////////////////////////////////////////////////
def getDepth_adv_g(df):

	if len(df.index) == 0:
		print('no record in %s VCF' % cellName)
	elif len(df.index) == 1:
		print('record found in %s VCF' % cellName)
		infoStr = df['INFO']
		infoStr = str(infoStr)
		DP = infoStr.split('DP')[1].split(';')[0].strip('=')
		print('sequencing depth: %s' % DP)
	else:
		print('multiple records found in %s gVCF' % cellName)
		infoDF = df['INFO']

		for i in range(0, len(infoDF.index)-1):
			line = infoDF.iloc[i]
			line = str(line)
			try:
				DP = line.split('DP')[1].split(';')[0].strip('=')
				print('       sequencing depth (record %d): %s' % (i, DP))
			except IndexError:
				continue

	print(' ')

#////////////////////////////////////////////////////////////////////
# getDepth()
#	given a single record, returns depth of coverage to that record
#		for VCF file
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
#	given a single record, returns depth of coverage to that record
#		for gVCF file
#////////////////////////////////////////////////////////////////////
def getDepth_g(df):

	if len(df.index) != 0: 
		print('record found in %s gVCF' % cellName)
		infoStr = df['INFO']
		infoStr = str(infoStr)
		print(infoStr)
		DP = infoStr.split('DP')[1].split(';')[0].strip('=')
		print('sequencing depth: %s' % DP)

	else:
		print('no record found in %s gVCF' % cellName)

	print(' ')

#////////////////////////////////////////////////////////////////////
# getGOI()
#	define a list of records corresponding to the GOI
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
#
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

# get a list of the records we actually care about
toKeepList_v = vcf.apply(getGOI_record, axis=1, args=(chrom_, start_ ,end_))
toKeepList_g = gvcf.apply(getGOI_record, axis=1, args=(chrom_, start_, end_))

# subset by relevant records
vcf_GOI = vcf[np.array(toKeepList_v, dtype=bool)]
gvcf_GOI = gvcf[np.array(toKeepList_g, dtype=bool)]

# get depth of coverage, for relevant records
getDepth_adv(vcf_GOI)
getDepth_adv_g(gvcf_GOI)

#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
