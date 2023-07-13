# Simulation replication


## Step 1: Run simulation

- block_structure_vary_np.R: File to run kelp for a particular simulation scenario on simulated data. Corresponding batch files in batch/simulation. 
- hierarchical_union.R: File to run simulation with hierarchically structured phenotypes. Corresponding batch files in batch/simulation. 
- emlkf_simulation.R: File to run (e-)MLKF simulations, the corresponding batch files are in batch/simulation. 
- UKB_simulation.R: Code for simulations on the UK Biobank data (the data is not available; we have applied for it from a source (the UKBiobank) and interested parties can also apply to the same source).

## Step 2: Evaluate simulation results

- summarize_np_block.R: Collect simulation results for kelp / fbh and produce power / fdr plots based on output of block_structure_vary_np.R.
- summarize_hierarchical_simulations: Collect simulation results for hierarchically structured output and produce power / fdr plots based on hierarchical_union.R.
- summarize_emlkf_simulation_results: Collect simulation results for (e-)MLKF output and produce power / fdr plots based on emlkf_simulation.R.
