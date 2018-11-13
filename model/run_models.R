setwd('~/Software/DAD/DAandConnectomes/model/')
source('./fit_models_connectome_MAP.R')

# Run all the models

HOME='/home/benjamin.garzon/'
DIR='modellingMAP'

STAN_FILE='age_model_single.stan'

NCLUSTERS=20
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
  'Data/DAD/processed/PET/modellingMAP',
  'Data/DAD/processed/VBM/modellingMAP',
  
  'Data/DAD/processedNOGSR/RS/modelling_modules_20MAP',
  'Data/DAD/processedNOGSR/GNG/modelling_modules_20MAP',
  'Data/DAD/processedNOGSR/TAB/modelling_modules_20MAP',
  'Data/DAD/processedNOGSR/RS/modelling_modules_70MAP',
  'Data/DAD/processedNOGSR/GNG/modelling_modules_70MAP',
  'Data/DAD/processedNOGSR/TAB/modelling_modules_70MAP',
  'Data/DAD/processedNOGSR/RS/modellingMAP',
  'Data/DAD/processedNOGSR/GNG/modellingMAP',
  'Data/DAD/processedNOGSR/TAB/modellingMAP',
  
  'Data/DAD/processedGSR/RS/modelling_modules_20MAP',
  'Data/DAD/processedGSR/GNG/modelling_modules_20MAP',
  'Data/DAD/processedGSR/TAB/modelling_modules_20MAP',
  'Data/DAD/processedGSR/RS/modelling_modules_70MAP',
  'Data/DAD/processedGSR/GNG/modelling_modules_70MAP',
  'Data/DAD/processedGSR/TAB/modelling_modules_70MAP',
  'Data/DAD/processedGSR/RS/modellingMAP',
  'Data/DAD/processedGSR/GNG/modellingMAP',
  'Data/DAD/processedGSR/TAB/modellingMAP'

)

n_models = length(INPUT_FILES)

for (m in seq(n_models)) {
  dir.create(file.path(HOME, OUTPUT_DIRS[m]))
  fit_models_connectome(
  file.path(HOME, INPUT_FILES[m]),
  file.path(HOME, OUTPUT_DIRS[m], 'params.csv'), 
  STAN_FILE,
  PRIORS_FILE, 
  FITFCS[m], 
  NCLUSTERS)
}