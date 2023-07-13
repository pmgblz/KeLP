

**IMPORTANT: Data is not available. We have applied for the data from the UK Biobank. Interested parties can do the same.**

- ## Step 1: Get Phenotypes + Merge with Genotypes

- ## Step 2: Run Lasso for both populations

	- lasso_ukb.R: Run the lasso on knockoffs and original data. 

- ## Step 3: Tune Gamma on White-Non British Population 

	- Run tuning: tune_gamma.R: Script to tune gamma for a particular population, for a particular gamma. 

	- Collect tuning results: find_the_best_gamma.R: Collects simulation results and produces .tex overview table showing the number of rejections for each gamma for each alpha. 

- ## Step 4: Use tuned gamma on British population and run KelP: 

	- eblipr_UKB.R: Run kelp for a particular population, a particular alpha and gamma. 

- ## Step 5: Evaluate Results

	- comparison_yengo_nature_2022_HG19_knockoff.R: Compare height rejections from the knockoff outer nodes with SNPs found in Yengo (2022, Nature)
	- comparison_yengo_nature_2022_HG19.R: Compare height rejections from kelp with SNPs found in Yengo (2022, Nature)



