#/////////////////////////////////////////////////////////////////////////
# script: createFinalOutDir.py
# author: Lincoln
# date: 3.18.19
#
# want to make one big dir with all of my filtered/unfiltered cells: 
#         scVCF_filtered_all/
# would be much nicer to do this on a jupyter notebook, buttt thats harder from EC2
# potential problem here is gonna be that we've got vcf AND csv files in this 
#    output dir
#/////////////////////////////////////////////////////////////////////////
import os 
import shutil 

filterDir = '/home/ubuntu/code/SNP_calling_pipeline/bulkAnalysis/filteredOut/'
filterDir_list = os.listdir(filterDir)

filteredCells = []
for f in filterDir_list:
	cell = f.strip('_unique.vcf')
	filteredCells.append(cell)


epiDir = '/home/ubuntu/code/SNP_calling_pipeline/bulkAnalysis/scVCF/'
epiDir_list = os.listdir(epiDir)

epiCells = []
for f in epiDir_list:
	cell = f.strip('.vcf')
	epiCells.append(cell)

# get cells in epiCells but NOT filteredCells
nonFilteredCells = set(epiCells) - set(filteredCells)

nonFilteredCells_list = []
for cell in nonFilteredCells:
	f = cell + '.vcf' 
	nonFilteredCells_list.append(f)


# copy over the non-filtered cells
outPATH = '/home/ubuntu/code/SNP_calling_pipeline/bulkAnalysis/scVCF_filtered_all/'
for file in nonFilteredCells_list:
	src = epiDir + file
	dst = outPATH + file
	shutil.copyfile(src, dst)

# copy over all the filtered cells
for file in filterDir_list:
	f = file.strip('_unique.vcf')
	f = f + '.vcf'
	scr = filterDir + file
	dst = outPATH + f
	shutil.copyfile(src, dst)

#/////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////