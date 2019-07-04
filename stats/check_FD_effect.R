
# statistical analysis

rm(list=ls())

# RUN SETUP SCRIPTS --------------------------------------- 
setwd('~/Software/DAD/DAandConnectomes')
source("./stats/aux_funcs.R")
sync_modules = F
source("./stats/compute_modules.R")
setwd('~/Software/DAD/DAandConnectomes')

# if not run
#source("./model/run_models.R")
#source("./model/evaluate_models_MAP.R")

# PATH DEFINITIONS --------------------------------------- 

CONN_FOLDER = "" # connectomes 


params = read.table('/home/benjamin.garzon/Data/DAD/processedNOGSR_FD/TAB/modellingMAP/params.csv', sep =',', header = T)
hist(params$rho, 50)
params = read.table('/home/benjamin.garzon/Data/DAD/processed/VBM/modellingMAP_FD/params.csv', sep =',', header = T)
hist(params$rho, 50)
params = read.table('/home/benjamin.garzon/Data/DAD/processed/PET/modellingMAP_FD/params.csv', sep =',', header = T)
hist(params$rho, 50)


#NO GSR
FIGS_DIR = "~/Data/DAD/analysis_variability/figsNOGSR/"

MOTION_DIR.TAB='~/Data/DAD/processed/TAB/RealignParameter/'
MOTION_DIR.GNG='~/Data/DAD/processed/GNG/RealignParameter/'
MOTION_DIR.RS='~/Data/DAD/processed/RS/RealignParameter/'

LABELS_FILE="~/Data/DAD/parcellations/shen/fconn_150_labels.txt"

DEMO_FILE = '~/Data/DAD/behaviour/demo_data.csv'

labels = read.csv(LABELS_FILE, header=TRUE, sep='\t')


# ANALYSES --------------------------------------- 
demo = read.table(DEMO_FILE, header=TRUE)

# remove cerebellum ROIS from analysis of PET data
cerebellum_rois = which(labels$LabelMNI=='Cerebellum')
means.PET$intercept_muc[cerebellum_rois] = means.PET$slope_muc[cerebellum_rois] = means.PET$log_intercept_sigmac[cerebellum_rois] = means.PET$log_slope_sigmac[cerebellum_rois] = NA

# analyze motion
params.TAB = analyze_motion(MOTION_DIR.TAB, demo)
params.GNG = analyze_motion(MOTION_DIR.GNG, demo)
params.RS = analyze_motion(MOTION_DIR.RS, demo)

dev.off()

calcpvals = function(mydatax, mydatay, mydataz, main){
  mydata = merge(mydatax[,c("Subject.ID", "mean.FD_Power")], mydatay, by = "Subject.ID", by.y = "Subject")  
  mydata = merge(mydataz[, c("Subject", "age")] , mydata, by.x = "Subject", by.y = "Subject.ID")
  pvals.age = apply( mydata[, -seq(4)], 2, function(x) pcor.test(x, 
                                                              mydata$mean.FD_Power, 
                                                              mydata$age)$estimate) # p.value
  pvals = apply( mydata[, -seq(4)], 2, function(x) cor.test(x, mydata$mean.FD_Power)$estimate) #p.value
  hist(pvals, 20, main = main, xlab = "p-value")
  hist(pvals.age, 20, main = paste(main, "\n (age corrected)"), xlab = "p-value")
  return(pvals)
}

params.mean = merge(params.TAB[,c("Subject.ID", "mean.FD_Power")], 
                    params.GNG[,c("Subject.ID", "mean.FD_Power")], 
                    by = "Subject.ID", all = T, suffixes = c(".TAB", ".GNG"))
params.mean = merge(params.mean, params.RS[,c("Subject.ID", "mean.FD_Power")], 
                    by = "Subject.ID", all = T, suffixes = c("", ".RS"))

params.mean = params.mean %>% dplyr::rename(mean.FD_Power.RS = mean.FD_Power)

params.mean$mean.FD_Power = rowMeans(params.mean[-1], na.rm = T)

params.complete = params.mean[complete.cases(params.mean), ]
cor(params.complete[-1])

# save all these
write.table(params.TAB[,c("Subject.ID", "mean.FD_Power")] %>% 
              dplyr::rename(Subject = Subject.ID,
                            FD = mean.FD_Power),
            file = "~/Data/DAD/FD.TAB.csv")
write.table(params.GNG[,c("Subject.ID", "mean.FD_Power")] %>% 
              dplyr::rename(Subject = Subject.ID,
                            FD = mean.FD_Power),
            file = "~/Data/DAD/FD.GNG.csv")
write.table(params.RS[,c("Subject.ID", "mean.FD_Power")] %>% 
              dplyr::rename(Subject = Subject.ID,
                            FD = mean.FD_Power),
            file = "~/Data/DAD/FD.RS.csv")
write.table(params.mean[,c("Subject.ID", "mean.FD_Power")] %>% 
              dplyr::rename(Subject = Subject.ID,
                            FD = mean.FD_Power),
            file = "~/Data/DAD/FD.mean.csv")


################
# PET
################
PET_DATA = read.table('~/Data/DAD/processed/PET/pet_data_ROI.csv') %>% dplyr::rename(Subject = V279)

par(mfrow = c(2, 2))
PET.pvals = calcpvals(params.mean, PET_DATA, demo, "PET")

################
# VBM
################
VBM_DATA = read.table('~/Data/DAD/processed/VBM/vbm_data.csv')%>% dplyr::rename(Subject = V279)
VBM.pvals = calcpvals(params.mean, VBM_DATA, demo, "VBM")


################
# FC NOGSR
################

openmat = function(input_file){
  Conn_info= readMat(input_file)
  labels = unlist(Conn_info$labels)
  X = data.frame(Conn_info$merged.matrices.mat)
  X$Subject  = subjects_FC.1 = unlist(Conn_info$subjects)
  return(X)
}

TAB_DATA.NOGSR = openmat('~/Data/DAD/processedNOGSR/TAB/zFC_all_150_valid.mat')
GNG_DATA.NOGSR = openmat('~/Data/DAD/processedNOGSR/GNG/zFC_all_150_valid.mat')
RS_DATA.NOGSR = openmat('~/Data/DAD/processedNOGSR/RS/zFC_all_150_valid.mat')

par(mfrow = c(3, 2))
TAB.NOGSR.pvals = calcpvals(params.TAB,TAB_DATA.NOGSR, demo, "FC TAB - NO GSR")
GNG.NOGSR.pvals = calcpvals(params.GNG, RS_DATA.NOGSR, demo, "FC GNG - NO GSR")
RS.NOGSR.pvals = calcpvals(params.RS, RS_DATA.NOGSR, demo, "FC RS - NO GSR")


################
# FC NOGSR
################

TAB_DATA.GSR = openmat('~/Data/DAD/processedGSR/TAB/zFC_all_150_valid.mat')
GNG_DATA.GSR = openmat('~/Data/DAD/processedGSR/GNG/zFC_all_150_valid.mat')
RS_DATA.GSR = openmat('~/Data/DAD/processedGSR/RS/zFC_all_150_valid.mat')

par(mfrow = c(3, 2))
TAB.GSR.pvals = calcpvals(params.TAB,TAB_DATA.NOGSR, demo, "FC TAB - GSR")
GNG.GSR.pvals = calcpvals(params.GNG, RS_DATA.NOGSR, demo, "FC GNG - GSR")
RS.GSR.pvals = calcpvals(params.RS, RS_DATA.NOGSR, demo, "FC RS - GSR")
