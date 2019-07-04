#rm(list=ls())
setwd('~/Software/DAD/DAandConnectomes/model/')
source("../stats/aux_funcs.R")

parnames = c('intercept_muc', 
             'slope_muc', 
             'log_intercept_sigmac', 
             'log_slope_sigmac')

MODULES_FILE_20="~/Data/DAD/ICNs/Laird/modules.RData"
MODULES_FILE_70="~/Data/DAD/ICNs/Ray/modules.RData"

#OUTPUT_FILE="~/Data/DAD/analysis_variability/dataScrubNOGSR_FD.Rdata"
#WDFC = "~/Data/DAD/processedNOGSR_FD/"

OUTPUT_FILE="~/Data/DAD/analysis_variability/dataScrubGSR_FD.Rdata"
WDFC = "~/Data/DAD/processedGSR_FD/"

getparams = function(WD, TASK, MODEL_DIR, MODULES_FILE=NULL, IS_FC = FALSE, bymodules= FALSE){
    print('################################')
    print(paste(MODEL_DIR, TASK))
    print('################################')
   
    FIGS_DIR <<- file.path(WD, TASK, MODEL_DIR, "figs/")
    PARAMS_FILE = file.path(WD, TASK, MODEL_DIR, 'params.csv')
    
    params = read.csv(PARAMS_FILE)
    
    means.params = colMeans(params)
    quants.params = apply(params, 2, quantile,  c(.025, .975))
    print('Means')
    print(means.params)
#    print('Quantiles')
#    print(quants.params)
    
    res = list(
      means = params,
      zeros = params == 0
    )
    
    if (IS_FC) {
      res$mu.20 = res$means$intercept_muc
      res$mu.70 = res$means$intercept_muc + res$means$slope_muc * 50
      res$sigma.20 = res$means$log_intercept_sigmac
      res$sigma.70 = res$means$log_intercept_sigmac + res$means$log_slope_sigmac * 50
    }
    
    return(res)    
}


MODEL_DIR = "modellingMAP_FD"
USEMCMC = T

WD="~/Data/DAD/processed/"
res = getparams(WD, 'PET', MODEL_DIR); means.PET = res$means; zeros.PET = res$zeros
res = getparams(WD, 'VBM', MODEL_DIR); means.VBM = res$means; zeros.VBM = res$zeros

WD = WDFC

MODEL_DIR="modelling_modules_20MAP"
res = getparams(WD, 'GNG', MODEL_DIR, MODULES_FILE = MODULES_FILE_20, IS_FC = T, bymodules=T); means.mod.20.GNG = res$means; zeros.mod.20.GNG = res$zeros 
res = getparams(WD, 'TAB', MODEL_DIR, MODULES_FILE = MODULES_FILE_20, IS_FC = T, bymodules=T); means.mod.20.TAB = res$means; zeros.mod.20.TAB = res$zeros
res = getparams(WD, 'RS', MODEL_DIR, MODULES_FILE = MODULES_FILE_20, IS_FC = T, bymodules=T);  means.mod.20.RS = res$means;  zeros.mod.20.RS = res$zeros

MODEL_DIR="modelling_modules_70MAP"
res = getparams(WD, 'GNG', MODEL_DIR, MODULES_FILE = MODULES_FILE_70, IS_FC = T, bymodules=T); means.mod.70.GNG = res$means; zeros.mod.70.GNG = res$zeros
res = getparams(WD, 'TAB', MODEL_DIR, MODULES_FILE = MODULES_FILE_70, IS_FC = T, bymodules=T); means.mod.70.TAB = res$means; zeros.mod.70.TAB = res$zeros
res = getparams(WD, 'RS', MODEL_DIR, MODULES_FILE = MODULES_FILE_70, IS_FC = T, bymodules=T);  means.mod.70.RS = res$means;  zeros.mod.70.RS = res$zeros

MODEL_DIR="modellingMAP"
res = getparams(WD, 'GNG', MODEL_DIR, MODULES_FILE = MODULES_FILE_20, IS_FC=T); means.GNG = res$means; zeros.GNG = res$zeros
FC.GNG.mu.20 = res$mu.20; FC.GNG.mu.70 = res$mu.70; FC.GNG.sigma.20 = res$sigma.20; FC.GNG.sigma.70 = res$sigma.70

res = getparams(WD, 'TAB', MODEL_DIR, MODULES_FILE = MODULES_FILE_20, IS_FC=T); means.TAB = res$means; zeros.TAB = res$zeros
FC.TAB.mu.20 = res$mu.20; FC.TAB.mu.70 = res$mu.70; FC.TAB.sigma.20 = res$sigma.20; FC.TAB.sigma.70 = res$sigma.70

res = getparams(WD, 'RS', MODEL_DIR, MODULES_FILE = MODULES_FILE_20, IS_FC=T);  means.RS = res$means;  zeros.RS = res$zeros
FC.RS.mu.20 = res$mu.20; FC.RS.mu.70 = res$mu.70; FC.RS.sigma.20 = res$sigma.20; FC.RS.sigma.70 = res$sigma.70
average.cog.RS = res$average.cog;  average.ICN.RS = res$average.ICN

save(means.PET, means.VBM, means.TAB, means.GNG, means.RS, 
     zeros.PET, zeros.VBM, zeros.TAB, zeros.GNG, zeros.RS,
     means.mod.20.TAB, means.mod.20.GNG, means.mod.20.RS, 
     zeros.mod.20.TAB, zeros.mod.20.GNG, zeros.mod.20.RS,
     means.mod.70.TAB, means.mod.70.GNG, means.mod.70.RS, 
     zeros.mod.70.TAB, zeros.mod.70.GNG, zeros.mod.70.RS,
     FC.RS.mu.20, FC.TAB.mu.20, FC.GNG.mu.20,
     FC.RS.mu.70, FC.TAB.mu.70, FC.GNG.mu.70,
     FC.RS.sigma.20, FC.TAB.sigma.20, FC.GNG.sigma.20,
     FC.RS.sigma.70, FC.TAB.sigma.70, FC.GNG.sigma.70,
     file = OUTPUT_FILE)


