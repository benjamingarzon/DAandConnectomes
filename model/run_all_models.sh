#!/bin/bash

# Run all the models
HOME="/home/benjamin.garzon/"
INITS=1
USEMCMC=1
FITFC=0
DIR=modellingMCMC

STAN_FILE="age_model_H.stan"

INPUT_FILE="$HOME/Data/DAD/processed/PET/pet_data_ROI.csv"
SAMPLES_DIR="$HOME/Data/DAD/processed/PET/$DIR/"
PRIORS_FILE="priors/PET_priors.csv"
#./run_models_parallel.sh $INPUT_FILE $SAMPLES_DIR $STAN_FILE $PRIORS_FILE $FITFC $USEMCMC $INITS

INPUT_FILE="$HOME/Data/DAD/processed/VBM/vbm_data.csv"
SAMPLES_DIR="$HOME/Data/DAD/processed/VBM/$DIR/"
PRIORS_FILE="priors/VBM_priors.csv"
#./run_models_parallel.sh $INPUT_FILE $SAMPLES_DIR $STAN_FILE $PRIORS_FILE $FITFC $USEMCMC $INITS

#INITS=100
USEMCMC=0
FITFC=1

DIR=modelling_noncentered
STAN_FILE="age_model_H.stan"

INPUT_FILE="$HOME/Data/DAD/processed/fmriprep/connectomes/TAB/zFC_all_150_valid.mat"
SAMPLES_DIR="$HOME/Data/DAD/processed/TAB/$DIR"
./run_models_parallel.sh $INPUT_FILE $SAMPLES_DIR $STAN_FILE $PRIORS_FILE $FITFC $USEMCMC $INITS &

#INPUT_FILE="$HOME/Data/DAD/processed/fmriprep/connectomes/GNG/zFC_all_150_valid.mat"
SAMPLES_DIR="$HOME/Data/DAD/processed/GNG/$DIR/"
#./run_models_parallel.sh $INPUT_FILE $SAMPLES_DIR $STAN_FILE $PRIORS_FILE $FITFC $USEMCMC $INITS &

#INPUT_FILE="$HOME/Data/DAD/processed/fmriprep/connectomes/RS/zFC_all_150_valid.mat"
SAMPLES_DIR="$HOME/Data/DAD/processed/RS/$DIR/"
./run_models_parallel.sh $INPUT_FILE $SAMPLES_DIR $STAN_FILE $PRIORS_FILE $FITFC $USEMCMC $INITS


USEMCMC=1
STAN_FILE="age_model_H.stan"

DIR=modelling_modules_20

# now do it by modules
#old files INPUT_FILE="$HOME/Data/DAD/processed/RS/modules/zFC_all_150_0.4_RS_20.mat"

INPUT_FILE="$HOME/Data/DAD/processed/fmriprep/connectomes/RS/modules/zFC_all_150_RS_20.mat"
SAMPLES_DIR="$HOME/Data/DAD/processed/RS/$DIR/"
./run_models_parallel.sh $INPUT_FILE $SAMPLES_DIR $STAN_FILE $PRIORS_FILE $FITFC $USEMCMC $INITS

INPUT_FILE="$HOME/Data/DAD/processed/fmriprep/connectomes/TAB/modules/zFC_all_150_TAB_20.mat"
SAMPLES_DIR="$HOME/Data/DAD/processed/TAB/$DIR"
./run_models_parallel.sh $INPUT_FILE $SAMPLES_DIR $STAN_FILE $PRIORS_FILE $FITFC $USEMCMC $INITS

INPUT_FILE="$HOME/Data/DAD/processed/fmriprep/connectomes/GNG/modules/zFC_all_150_GNG_20.mat"
SAMPLES_DIR="$HOME/Data/DAD/processed/GNG/$DIR/"
./run_models_parallel.sh $INPUT_FILE $SAMPLES_DIR $STAN_FILE $PRIORS_FILE $FITFC $USEMCMC $INITS

DIR=modelling_modules_70

# now do it by modules
INPUT_FILE="$HOME/Data/DAD/processed/fmriprep/connectomes/RS/modules/zFC_all_150_RS_70.mat"
SAMPLES_DIR="$HOME/Data/DAD/processed/RS/$DIR/"
./run_models_parallel.sh $INPUT_FILE $SAMPLES_DIR $STAN_FILE $PRIORS_FILE $FITFC $USEMCMC $INITS

INPUT_FILE="$HOME/Data/DAD/processed/fmriprep/connectomes/TAB/modules/zFC_all_150_TAB_70.mat"
SAMPLES_DIR="$HOME/Data/DAD/processed/TAB/$DIR"
./run_models_parallel.sh $INPUT_FILE $SAMPLES_DIR $STAN_FILE $PRIORS_FILE $FITFC $USEMCMC $INITS

INPUT_FILE="$HOME/Data/DAD/processed/fmriprep/connectomes/GNG/modules/zFC_all_150_GNG_70.mat"
SAMPLES_DIR="$HOME/Data/DAD/processed/GNG/$DIR/"
./run_models_parallel.sh $INPUT_FILE $SAMPLES_DIR $STAN_FILE $PRIORS_FILE $FITFC $USEMCMC $INITS


# DIR=modelling_centered
# PRIORS_FILE="priors/FC_priors.csv"
# STAN_FILE="age_model_H_centered.stan"
# 
# # old files INPUT_FILE="$HOME/Data/DAD/processed/RS/Connectome0.4/zFC_all_150_0.4_RS_valid.mat"
# INPUT_FILE="$HOME/Data/DAD/processed/connectomes/RS/zFC_all_150_valid.mat"
# SAMPLES_DIR="$HOME/Data/DAD/processed/RS/$DIR/"
# #./run_models_parallel.sh $INPUT_FILE $SAMPLES_DIR $STAN_FILE $PRIORS_FILE $FITFC $USEMCMC $INITS
# 
# INPUT_FILE="$HOME/Data/DAD/processed/connectomes/TAB/zFC_all_150_valid.mat"
# SAMPLES_DIR="$HOME/Data/DAD/processed/TAB/$DIR"
# #./run_models_parallel.sh $INPUT_FILE $SAMPLES_DIR $STAN_FILE $PRIORS_FILE $FITFC $USEMCMC $INITS
# 
# INPUT_FILE="$HOME/Data/DAD/processed/connectomes/GNG/zFC_all_150_valid.mat"
# SAMPLES_DIR="$HOME/Data/DAD/processed/GNG/$DIR/"
# #./run_models_parallel.sh $INPUT_FILE $SAMPLES_DIR $STAN_FILE $PRIORS_FILE $FITFC $USEMCMC $INITS
