# DAandConnectomes
Analysis scripts for article: 

Main scripts:

stats/analysis.R: Statistical analyses and figure generation

stats/aux_funcs.R: Auxiliary functions

stats/compute_modules.R: Average parameter models over anatomical or network modules

Model:

model/run_all_models.sh, model/run_models_parallel.sh: scripts to run models and reorgarnize the large sample files

model/fit_models_connectome.R: wraps the stan models

model/age_model_H.stan, model/age_model_H_centered.stan: Stan programs for age model