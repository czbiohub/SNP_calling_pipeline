# Makefile
# authors: Lincoln Harris
# date: 2.11.18

test_kras:
	ipython checkCoverage_batch.py 12 25209431 25250803 cellsList_subset.csv testOut_kras.csv

test_actin:
	ipython checkCoverage_batch.py 7 5528150 5530604 cellsList_subset.csv testOut_actin.csv

test_gapdh:
	ipython checkCoverage_batch.py 12 6533927 6538371 cellsList_subset.csv testOut_gapdh.csv

test_egfr:
	ipython checkCoverage_batch.py 7 55152337 55207337 cellsList_subset.csv testOut_egfr.csv

test_egfr_L858R:
	ipython checkCoverage_batch.py 7 55154909 55154911 cellsList.csv testOut_egfr_L858R.csv

run_gapdh_parallel:
	ipython checkCoverage_parallel.py 12 6533927 6538371 14 cellsList.csv gapdh_allCells.csv

run_egfr_L858R_parallel:
	ipython checkCoverage_parallel.py 7 55154909 55154911 14 cellsList.csv egfr_L858R_allCells.csv 

