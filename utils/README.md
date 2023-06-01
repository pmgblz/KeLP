# Files in alphabetical order: 

- ebh_aux: Functions for e-BH and focused e-BH

- emlkf_aux: Functions for e-MLKF. 

- filter_aux: Contains functions for filtering rejections to outer nodes and comparing to levels of resolution and picking the group with highest e-value / fraction. 

- kelp_aux: Function for running kelp, some code snippets inspired by https://github.com/amspector100/blipr/blob/main/R/blip.R
	1. Function to create location constraint matrix 
	2. Function to run kelp
	3. Helper functions for kelp (e.g. for calculating inverse size)


- knockoff_aux: Functions to contruct knockoff variables
	 1. Functions to construct knockoffs with Julia 
	 2. Functions to construct knockoffs for (e-)MLKF
	 3. Construct knockoff e-values for UKB and knockoff helper functions


- load_aux: Script to load pacakges and functions


- multiplier_tuning_aux.R: Function to generate a random matrix with fixed rowsum, which can be used to tune the multiplier "c" (Note: This is not used in our simulations.)


- plot_aux: Function for plotting UKB rejections, using snippets from https://github.com/msesia/knockoffgwas/blob/master/visualization/utils_plotting.R
   	 
- power_fdr_aux: Helper functions for calculating power and fdr.

- sim_data_aux: Functions to simulate data:
	1. Function to group structure
	2. Function to generate data (X, Y and Sigma).
	3. Functions to generate beta-matrix (completely at random, choosing which largest groups should be nonzero, for hierarchical outcome setting)
	4. Function generate covariance matrix with block structure, based on code snippets from  https://github.com/ekatsevi/Focused-BH/blob/c78954d3cf081fbcf22d3d7c2d26c64c17f7c7bb/src/archive/aux_make_design.R
	5. Helper functions for covariance generation



- tune_gamma_sim_aux.R: Functions to tume gamma: 
	1. In genereral simulation setting. 
	2. In hierarchical outcome setting. 
	3. For e-MLKF. 

