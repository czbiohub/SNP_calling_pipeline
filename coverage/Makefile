# Makefile
# authors: Lincoln Harris
# date: 2.11.18

test_kras:
	ipython checkCoverage_batch.py 12 25209431 25250803 cellsList_subset.csv testOut_kras.csv

test_actin:
	ipython checkCoverage_batch.py 8 5528150 5530604 cellsList_subset.csv testOut_actin.csv

test_gapdh:
	ipython checkCoverage_batch.py 12 6533927 6538371 cellsList_subset.csv testOut_gapdh.csv

test_egfr:
	ipython checkCoverage_parallel.py 7 55152337 55207337 16 cellsList_subset.csv testOut_egfr.csv

test_egfr_L858R:
	ipython checkCoverage_batch.py 7 55154909 55154911 cellsList.csv testOut_egfr_L858R.csv

test_kras_parallel:
	ipython checkCoverage_parallel.py 12 25209431 25250803 16 cellsList_subset.csv kras_subset.csv

test_gapdh_parallel:
	ipython checkCoverage_parallel.py 12 6533927 6538371 16 cellsList_subset.csv gapdh_subset.csv

test_egfr_L858R_parallel:
	ipython checkCoverage_parallel.py 7 55191821 55191822 16 cellsList_subset.csv egfr_L858R_subset.csv 

run_gapdh_parallel:
	ipython checkCoverage_parallel.py 12 6533927 6538371 16 cellsList.csv gapdh_allCells.csv

run_egfr_L858R_parallel:
	ipython checkCoverage_parallel.py 7 55191821 55191822 16 cellsList.csv egfr_L858R_allCells.csv 

run_egfr_T790M_parallel:
	ipython checkCoverage_parallel.py 7 55181377 55181379 16 cellsList.csv egfr_T790M_allCells.csv 

run_egfr_L861Q_parallel:
	ipython checkCoverage_parallel.py 7 55191830 55191832 16 cellsList.csv egfr_L861Q_allCells.csv 

run_egfr_G719X_parallel:
	ipython checkCoverage_parallel.py 7 55174014 55174016 16 cellsList.csv egfr_G719X_allCells.csv 

run_egfr_S768I_parallel:
	ipython checkCoverage_parallel.py 7 55181311 55181313 16 cellsList.csv egfr_S768I_allCells.csv 

