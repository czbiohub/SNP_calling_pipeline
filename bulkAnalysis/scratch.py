import os 

filterDir = '/home/ubuntu/code/SNP_calling_pipeline/bulkAnalysis/filteredOut/'
filterDir_list = os.listdir(filterDir)

for f in filterDir_list:
	cell = f.strip('_unique.csv')
	print(cell)