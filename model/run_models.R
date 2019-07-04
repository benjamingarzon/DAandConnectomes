setwd('~/Software/DAD/DAandConnectomes/model/')
source('./fit_models_connectome_MAP_covariate.R')

# Run all the models

HOME='/home/benjamin.garzon/'
DIR='modellingMAP'

STAN_FILE='age_model_single_FD.stan'

NCLUSTERS=10
PRIORS_FILE='./priors/priors_MAP.csv'

FITFCS = c(rep(0, 2), rep(1, 18))

INPUT_FILES = c(
  'Data/DAD/processed/PET/pet_data_ROI.csv',
  'Data/DAD/processed/VBM/vbm_data.csv',
  
  'Data/DAD/processedNOGSR/RS/modules/zFC_all_150_RS_20.mat',
  'Data/DAD/processedNOGSR/GNG/modules/zFC_all_150_GNG_20.mat',
  'Data/DAD/processedNOGSR/TAB/modules/zFC_all_150_TAB_20.mat',
  'Data/DAD/processedNOGSR/RS/modules/zFC_all_150_RS_70.mat',
  'Data/DAD/processedNOGSR/GNG/modules/zFC_all_150_GNG_70.mat',
  'Data/DAD/processedNOGSR/TAB/modules/zFC_all_150_TAB_70.mat',
  'Data/DAD/processedNOGSR/RS/zFC_all_150_valid.mat',
  'Data/DAD/processedNOGSR/GNG/zFC_all_150_valid.mat',
  'Data/DAD/processedNOGSR/TAB/zFC_all_150_valid.mat',
  
  'Data/DAD/processedGSR/RS/modules/zFC_all_150_RS_20.mat',
  'Data/DAD/processedGSR/GNG/modules/zFC_all_150_GNG_20.mat',
  'Data/DAD/processedGSR/TAB/modules/zFC_all_150_TAB_20.mat',
  'Data/DAD/processedGSR/RS/modules/zFC_all_150_RS_70.mat',
  'Data/DAD/processedGSR/GNG/modules/zFC_all_150_GNG_70.mat',
  'Data/DAD/processedGSR/TAB/modules/zFC_all_150_TAB_70.mat',
  'Data/DAD/processedGSR/RS/zFC_all_150_valid.mat',
  'Data/DAD/processedGSR/GNG/zFC_all_150_valid.mat',
  'Data/DAD/processedGSR/TAB/zFC_all_150_valid.mat'
)

OUTPUT_DIRS = c(
  'Data/DAD/processed/PET/modellingMAP_FD',
  'Data/DAD/processed/VBM/modellingMAP_FD',
  
  'Data/DAD/processedNOGSR_FD/RS/modelling_modules_20MAP',
  'Data/DAD/processedNOGSR_FD/GNG/modelling_modules_20MAP',
  'Data/DAD/processedNOGSR_FD/TAB/modelling_modules_20MAP',
  'Data/DAD/processedNOGSR_FD/RS/modelling_modules_70MAP',
  'Data/DAD/processedNOGSR_FD/GNG/modelling_modules_70MAP',
  'Data/DAD/processedNOGSR_FD/TAB/modelling_modules_70MAP',
  'Data/DAD/processedNOGSR_FD/RS/modellingMAP',
  'Data/DAD/processedNOGSR_FD/GNG/modellingMAP',
  'Data/DAD/processedNOGSR_FD/TAB/modellingMAP',
  
  'Data/DAD/processedGSR_FD/RS/modelling_modules_20MAP',
  'Data/DAD/processedGSR_FD/GNG/modelling_modules_20MAP',
  'Data/DAD/processedGSR_FD/TAB/modelling_modules_20MAP',
  'Data/DAD/processedGSR_FD/RS/modelling_modules_70MAP',
  'Data/DAD/processedGSR_FD/GNG/modelling_modules_70MAP',
  'Data/DAD/processedGSR_FD/TAB/modelling_modules_70MAP',
  'Data/DAD/processedGSR_FD/RS/modellingMAP',
  'Data/DAD/processedGSR_FD/GNG/modellingMAP',
  'Data/DAD/processedGSR_FD/TAB/modellingMAP'

)

COVARIATE_FILES = c(
  '~/Data/DAD/FD.mean.csv',
  '~/Data/DAD/FD.mean.csv',

  '~/Data/DAD/FD.RS.csv',
  '~/Data/DAD/FD.GNG.csv',
  '~/Data/DAD/FD.TAB.csv',
  '~/Data/DAD/FD.RS.csv',
  '~/Data/DAD/FD.GNG.csv',
  '~/Data/DAD/FD.TAB.csv',
  '~/Data/DAD/FD.RS.csv',
  '~/Data/DAD/FD.GNG.csv',
  '~/Data/DAD/FD.TAB.csv',
  
  '~/Data/DAD/FD.RS.csv',
  '~/Data/DAD/FD.GNG.csv',
  '~/Data/DAD/FD.TAB.csv',
  '~/Data/DAD/FD.RS.csv',
  '~/Data/DAD/FD.GNG.csv',
  '~/Data/DAD/FD.TAB.csv',
  '~/Data/DAD/FD.RS.csv',
  '~/Data/DAD/FD.GNG.csv',
  '~/Data/DAD/FD.TAB.csv'
)


n_models = length(INPUT_FILES)

#for (m in seq(n_models)) {
for (m in seq(2)) {
    dir.create(file.path(HOME, OUTPUT_DIRS[m]))
  fit_models_connectome(
  file.path(HOME, INPUT_FILES[m]),
  file.path(HOME, OUTPUT_DIRS[m], 'params.csv'), 
  STAN_FILE,
  PRIORS_FILE, 
  FITFCS[m], 
  NCLUSTERS,
  COVARIATE_FILES[m])
}