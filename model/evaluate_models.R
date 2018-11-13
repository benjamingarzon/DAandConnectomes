#rm(list=ls())
source("./stats/aux_funcs.R")

params = c("intercept_muc_mu", "intercept_muc_sigma", 
           "slope_muc_mu", "slope_muc_sigma", 
           "log_intercept_sigmac_mu", "log_intercept_sigmac_sigma",
           "log_slope_sigmac_mu", "log_slope_sigmac_sigma", "lp__")

#params = c("intercept_muc_mu", "intercept_muc_sigma", 
#           "log_intercept_sigmac_mu", "log_intercept_sigmac_sigma",
#           "lp__")

MODULES_FILE_20="~/Data/DAD/ICNs/Laird/modules.RData"
MODULES_FILE_70="~/Data/DAD/ICNs/Ray/modules.RData"

OUTPUT_FILE="~/Data/DAD/processed/analysis_variability/dataScrubGSR.Rdata"

#SUFFIX=""
#OUTPUT_FILE="~/Data/DAD/processed/analysis_variability/data_GSR.Rdata"


dochecks = function(WD, TASK, USEMCMC, MODEL_DIR, MODULES_FILE=NULL, IS_FC = FALSE, bymodules= FALSE){
  
  FIGS_DIR <<- file.path(WD, TASK, MODEL_DIR, "figs/")
  
  if (!USEMCMC) {
    SAMPLES_FILE="vbsamples"
    
  } else { 
    PARAMS_FILE_LIST=file.path(WD, TASK, MODEL_DIR, paste0("vbsamples", seq(4), "_simplified") )
    check_chains(FIGS_DIR, PARAMS_FILE_LIST, TASK)
    SAMPLES_FILE="vbsamples"
  }
  
  PARAMS_FILE = file.path(WD, TASK, MODEL_DIR, SAMPLES_FILE)
  res = evaluate(FIGS_DIR, PARAMS_FILE, TASK)
  
  if (IS_FC) {
    res$mu.20 = predict_value(0, PARAMS_FILE)
    res$mu.70 = predict_value(50, PARAMS_FILE)
    res$sigma.20 = predict_value(0, PARAMS_FILE, musigma="sigma")
    res$sigma.70 = predict_value(50, PARAMS_FILE, musigma="sigma")
    res$average.cog = evaluate_averages(PARAMS_FILE, MODULES_FILE, ICN=FALSE, bymodules)
    res$average.ICN = evaluate_averages(PARAMS_FILE, MODULES_FILE, ICN=TRUE, bymodules)
  }
  return(res)
  
}

MODEL_DIR = "modellingMCMC"
USEMCMC = T

WD="~/Data/DAD/processed/"
res = dochecks(WD, 'PET', USEMCMC, MODEL_DIR); means.PET = res$means; zeros.PET = res$zeros
res = dochecks(WD, 'VBM', USEMCMC, MODEL_DIR); means.VBM = res$means; zeros.VBM = res$zeros

WD="~/Data/DAD/processedGSR/"
SUFFIX=""

MODEL_DIR=paste0("modelling_modules_20", SUFFIX)
res = dochecks(WD, 'GNG', USEMCMC, MODEL_DIR, MODULES_FILE = MODULES_FILE_20, IS_FC = T, bymodules=T); means.mod.20.GNG = res$means; zeros.mod.20.GNG = res$zeros 
average.cog.GNG.20 = res$average.cog; average.ICN.GNG.20 = res$average.ICN

res = dochecks(WD, 'TAB', USEMCMC, MODEL_DIR, MODULES_FILE = MODULES_FILE_20, IS_FC = T, bymodules=T); means.mod.20.TAB = res$means; zeros.mod.20.TAB = res$zeros
average.cog.TAB.20 = res$average.cog; average.ICN.TAB.20 = res$average.ICN

res = dochecks(WD, 'RS', USEMCMC, MODEL_DIR, MODULES_FILE = MODULES_FILE_20, IS_FC = T, bymodules=T);  means.mod.20.RS = res$means;  zeros.mod.20.RS = res$zeros
average.cog.RS.20 = res$average.cog; average.ICN.RS.20 = res$average.ICN

# MODEL_DIR=paste0("modelling_modules_70", SUFFIX)
# res = dochecks(WD, 'GNG', USEMCMC, MODEL_DIR, MODULES_FILE = MODULES_FILE_70, IS_FC = T, bymodules=T); means.mod.70.GNG = res$means; zeros.mod.70.GNG = res$zeros
# average.cog.GNG.70 = res$average.cog; average.ICN.GNG.70 = res$average.ICN
# 
# res = dochecks(WD, 'TAB', USEMCMC, MODEL_DIR, MODULES_FILE = MODULES_FILE_70, IS_FC = T, bymodules=T); means.mod.70.TAB = res$means; zeros.mod.70.TAB = res$zeros
# average.cog.TAB.70 = res$average.cog; average.ICN.TAB.70 = res$average.ICN
# 
# res = dochecks(WD, 'RS', USEMCMC, MODEL_DIR, MODULES_FILE = MODULES_FILE_70, IS_FC = T, bymodules=T);  means.mod.70.RS = res$means;  zeros.mod.70.RS = res$zeros
# average.cog.RS.70 = res$average.cog; average.ICN.RS.70 = res$average.ICN

USEMCMC = F
MODEL_DIR=paste0("modelling_noncentered", SUFFIX)
res = dochecks(WD, 'GNG', USEMCMC, MODEL_DIR, MODULES_FILE = MODULES_FILE_20, IS_FC=TRUE); means.GNG = res$means; zeros.GNG = res$zeros
FC.GNG.mu.20 = res$mu.20; FC.GNG.mu.70 = res$mu.70; FC.GNG.sigma.20 = res$sigma.20; FC.GNG.sigma.70 = res$sigma.70
average.cog.GNG = res$average.cog; average.ICN.GNG = res$average.ICN

res = dochecks(WD, 'TAB', USEMCMC, MODEL_DIR, MODULES_FILE = MODULES_FILE_20, IS_FC=TRUE); means.TAB = res$means; zeros.TAB = res$zeros
FC.TAB.mu.20 = res$mu.20; FC.TAB.mu.70 = res$mu.70; FC.TAB.sigma.20 = res$sigma.20; FC.TAB.sigma.70 = res$sigma.70
average.cog.TAB = res$average.cog; average.ICN.TAB = res$average.ICN

res = dochecks(WD, 'RS', USEMCMC, MODEL_DIR, MODULES_FILE = MODULES_FILE_20, IS_FC=TRUE);  means.RS = res$means;  zeros.RS = res$zeros
FC.RS.mu.20 = res$mu.20; FC.RS.mu.70 = res$mu.70; FC.RS.sigma.20 = res$sigma.20; FC.RS.sigma.70 = res$sigma.70
average.cog.RS = res$average.cog;  average.ICN.RS = res$average.ICN


save(means.PET, means.VBM, means.TAB, means.GNG, means.RS, 
     zeros.PET, zeros.VBM, zeros.TAB, zeros.GNG, zeros.RS,
     means.mod.20.TAB, means.mod.20.GNG, means.mod.20.RS, 
     zeros.mod.20.TAB, zeros.mod.20.GNG, zeros.mod.20.RS,
     # means.mod.70.TAB, means.mod.70.GNG, means.mod.70.RS, 
     # zeros.mod.70.TAB, zeros.mod.70.GNG, zeros.mod.70.RS,
     average.cog.TAB, average.cog.GNG, average.cog.RS, 
     average.ICN.TAB, average.ICN.GNG, average.ICN.RS,
     average.cog.TAB.20, average.cog.GNG.20, average.cog.RS.20, 
     average.ICN.TAB.20, average.ICN.GNG.20, average.ICN.RS.20,
     # average.cog.TAB.70, average.cog.GNG.70, average.cog.RS.70, 
     # average.ICN.TAB.70, average.ICN.GNG.70, average.ICN.RS.70,
     FC.RS.mu.20, FC.TAB.mu.20, FC.GNG.mu.20,
     # FC.RS.mu.70, FC.TAB.mu.70, FC.GNG.mu.70,
     FC.RS.sigma.20, FC.TAB.sigma.20, FC.GNG.sigma.20,
     # FC.RS.sigma.70, FC.TAB.sigma.70, FC.GNG.sigma.70,
     ## slope_muc.samples.RS, slope_muc.samples.TAB, slope_muc.samples.GNG,
     ## intercept_muc.samples.RS, intercept_muc.samples.TAB, intercept_muc.samples.GNG,
     file = OUTPUT_FILE)

dev.off()

