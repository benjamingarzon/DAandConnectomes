# DAandConnectomes
Analysis scripts for article: 

## Processing:

# FMRIPREP
processing/organize_BIDS.sh: Organize data in BIDS format

processing/run_fmriprep.sh: Launch fmriprep pipeline

processing/run_extract_connectivity.sh: Extract connectivity values using atlas 

processing/merge_matrices_shen.m: Merge matrices for all subjects in one file

# DPARSFA (1)
ConnectivityAnalysis/PipelineConnectomeGSR.sh
ConnectivityAnalysis/PipelineConnectomeNOGSR.sh

## Model:

# MCMC/variational
model/run_all_models.sh, model/run_models_parallel.sh: scripts to run models and reorgarnize the large sample files

model/fit_models_connectome.R: wraps the stan model fitting

model/evaluate_models.R: diagnostics on model fitting procedures

model/age_model_H.stan, model/age_model_H_centered.stan: Stan programs for age model, noncentered and centered versions (only noncentered was used)


# MAP estimation

model/run_models.R (3)

model/fit_models_connectome_MAP.R

model/age_model_single.stan: non-hierarchical model for individual connections/ROIs

model/evaluate_models_MAP.R: gather parameters in .csv files


## Main analysis scripts:

stats/analysis.R: Statistical analyses and figure generation

stats/aux_funcs.R: Auxiliary functions

stats/compute_modules.R (2): Average parameter models over anatomical or network modules




Order of execution would be