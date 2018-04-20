# statistical analysis

rm(list=ls())

# RUN SETUP SCRIPTS --------------------------------------- 

setwd('~/Software/DAD/DAandConnectomes')
source("./stats/aux_funcs.R")
source("./stats/compute_modules.R")
setwd('~/Software/DAD/DAandConnectomes')

# if not run
#system("../model/run_all_models.sh", wait=FALSE)
#source("../model/evaluate_models.R")

# PATH DEFINITIONS --------------------------------------- 

INPUT_FILE.TAB='~/Data/DAD/processed/TAB/Connectome0.3/zFC_all_150_0.3_TAB_valid.mat'
INPUT_FILE.GNG='~/Data/DAD/processed/GNG/Connectome0.3/zFC_all_150_0.3_GNG_valid.mat'
INPUT_FILE.RS='~/Data/DAD/processed/RS/Connectome0.4/zFC_all_150_0.4_RS_valid.mat'

MOTION_DIR.TAB='~/Data/DAD/processed/TAB/RealignParameter/'
MOTION_DIR.GNG='~/Data/DAD/processed/GNG/RealignParameter/'
MOTION_DIR.RS='~/Data/DAD/processed/RS/RealignParameter/'

LABELS_FILE="~/Data/DAD/parcellations/shen/fconn_150_labels.txt"

DEMO_FILE = '~/Data/DAD/behaviour/demo_data.csv'
DISTANCE_FILE = '~/Data/DAD/parcellations/shen/parc_shen_150.dist.csv'

labels = read.csv(LABELS_FILE, header=TRUE, sep='\t')

doplot = T # deactivate some plots that take a long time
FIGS_DIR = "~/Data/DAD/processed/analysis_variability/figs/"
#unlink(paste0(FIGS_DIR, '/*'))
fig_count = 0


# ANALYSES --------------------------------------- 
load("~/Data/DAD/processed/analysis_variability/data.Rdata")
demo = read.table(DEMO_FILE, header=TRUE)

# remove cerebellum ROIS from analysis of PET data
cerebellum_rois = which(labels$LabelMNI=='Cerebellum')
means.PET$intercept_muc[cerebellum_rois] = means.PET$slope_muc[cerebellum_rois] = means.PET$log_intercept_sigmac[cerebellum_rois] = means.PET$log_slope_sigmac[cerebellum_rois] = NA

# analyze demographics
tabulate_demo(demo)
tabulate_demo(demo, INPUT_FILE.GNG)
tabulate_demo(demo, INPUT_FILE.TAB)
tabulate_demo(demo, INPUT_FILE.GNG)

# analyze motion
params.TAB = analyze_motion(MOTION_DIR.TAB, demo, "motion.TAB")
params.GNG = analyze_motion(MOTION_DIR.GNG, demo, "motion.GNG")
params.RS = analyze_motion(MOTION_DIR.RS, demo, "motion.RS")

t.test(params.TAB$mean.FD_Power, params.GNG$mean.FD_Power)
t.test(params.RS$mean.FD_Power, params.GNG$mean.FD_Power)
t.test(params.TAB$mean.FD_Power, params.RS$mean.FD_Power)

t.test(params.TAB$mean.FD_Power_scrubbed, params.GNG$mean.FD_Power_scrubbed)
t.test(params.RS$mean.FD_Power_scrubbed, params.GNG$mean.FD_Power_scrubbed)
t.test(params.TAB$mean.FD_Power_scrubbed, params.RS$mean.FD_Power_scrubbed)


# COMPARISON BETWEEN PARAMETERS ----------------------

save_fig(res=BWRES)
#all_params = cbind(means.GNG$intercept_muc, means.GNG$slope_muc, exp(means.GNG$log_intercept_sigmac), exp(50*means.GNG$log_slope_sigmac) ) 
all_params.GNG = cbind(means.GNG$intercept_muc, means.GNG$slope_muc, means.GNG$log_intercept_sigmac, means.GNG$log_slope_sigmac ) 

colnames(all_params.GNG) = c("alpha_mu", "beta_mu", "alpha_sigma", "beta_sigma")
pairs( all_params.GNG[seq(1, nrow(all_params.GNG), 10), ], diag.panel=panel.hist, labels=c(expression(alpha[mu]),  expression(beta[mu]),  expression(alpha[sigma]),  expression(beta[sigma])), pch=20, cex=1)
cor(all_params.GNG)

save_fig(res=BWRES)
par(mfrow=c(3, 3), mar=c(8,8,5,5), mgp=MYMGP)
plot_correl(means.GNG$slope_muc, means.GNG$intercept_muc, expression(beta[mu] ~ "for FC (GNG)"), expression(alpha[mu] ~ "for FC (GNG)"), cex = MYCEX)
plot_correl(means.TAB$slope_muc, means.TAB$intercept_muc, expression(beta[mu] ~ "for FC (TAB)"), expression(alpha[mu] ~ "for FC (TAB)"), cex = MYCEX)
plot_correl(means.RS$slope_muc, means.RS$intercept_muc, expression(beta[mu] ~ "for FC (RS)"), expression(alpha[mu] ~ "for FC (RS)"), cex = MYCEX)

plot_correl(means.GNG$log_slope_sigmac, means.GNG$log_intercept_sigmac, expression(beta[sigma] ~ "for FC (GNG)"), expression(alpha[sigma] ~ "for FC (GNG)"), cex = MYCEX)
plot_correl(means.TAB$log_slope_sigmac, means.TAB$log_intercept_sigmac, expression(beta[sigma] ~ "for FC (TAB)"), expression(alpha[sigma] ~ "for FC (TAB)"), cex = MYCEX)
plot_correl(means.RS$log_slope_sigmac, means.RS$log_intercept_sigmac, expression(beta[sigma] ~ "for FC (RS)"), expression(alpha[sigma] ~ "for FC (RS)"), cex = MYCEX)

plot_correl(exp(50*means.GNG$log_slope_sigmac), exp(means.GNG$log_intercept_sigmac),   expression("exp(50*" ~ beta[sigma] ~ ") for FC (GNG)"),  expression( "exp(50*" ~ alpha[sigma] ~ ") for FC (GNG)"), cex = MYCEX)
plot_correl(exp(50*means.TAB$log_slope_sigmac), exp(means.TAB$log_intercept_sigmac),  expression( "exp(50*" ~ beta[sigma] ~ ") for FC (TAB)"),  expression( "exp(50*" ~ alpha[sigma] ~ ") for FC (TAB)"), cex = MYCEX)
plot_correl(exp(50*means.RS$log_slope_sigmac), exp(means.RS$log_intercept_sigmac),  expression( "exp(50*" ~ beta[sigma] ~ ") for FC (RS)"),  expression("exp(50*" ~ alpha[sigma] ~ ") for FC (RS)"), cex = MYCEX)


save_fig(res=BWRES)
par(mfrow=c(1, 3), mar=c(8,8,5,5), mgp=MYMGP)
plot_correl(means.GNG$slope_muc, means.GNG$log_slope_sigmac, expression(beta[mu] ~ "for FC (GNG)"), expression(beta[sigma] ~ "for FC (GNG)"), cex = MYCEX)
plot_correl(means.TAB$slope_muc, means.TAB$log_slope_sigmac, expression(beta[mu] ~ "for FC (TAB)"), expression(beta[sigma] ~ "for FC (TAB)"), cex = MYCEX)
plot_correl(means.RS$slope_muc, means.RS$log_slope_sigmac, expression(beta[mu] ~ "for FC (RS)"), expression(beta[sigma] ~ "for FC (RS)"), cex = MYCEX)

# COMPARISON BETWEEN TASKS ----------------------

save_fig(res=BWRES)
par(mfrow=c(2, 3), mar=c(8,8,5,5), mgp=MYMGP)
plot_correl(means.GNG$slope_muc, means.TAB$slope_muc, expression(beta[mu] ~ "for FC (GNG)"), expression(beta[mu] ~ "for FC (TAB)"), cex = MYCEX, asp=1)
plot_correl(means.TAB$slope_muc, means.RS$slope_muc, expression(beta[mu] ~ "for FC (TAB)"), expression(beta[mu] ~ "for FC (RS)"), cex = MYCEX, asp=1)
plot_correl(means.RS$slope_muc, means.GNG$slope_muc, expression(beta[mu] ~ "for FC (RS)"), expression(beta[mu] ~ "for FC (GNG)"), cex = MYCEX, asp=1)

plot_correl(means.GNG$log_slope_sigmac, means.TAB$log_slope_sigmac, expression(beta[sigma] ~ "for FC (GNG)"), expression(beta[sigma] ~ "for FC (TAB)"), cex = MYCEX, asp=1)
plot_correl(means.TAB$log_slope_sigmac, means.RS$log_slope_sigmac, expression(beta[sigma] ~ "for FC (TAB)"), expression(beta[sigma] ~ "for FC (RS)"), cex = MYCEX, asp=1)
plot_correl(means.RS$log_slope_sigmac, means.GNG$log_slope_sigmac, expression(beta[sigma] ~ "for FC (RS)"), expression(beta[sigma] ~ "for FC (GNG)"), cex = MYCEX, asp=1)

save_fig(res=BWRES)

par(mfrow=c(2, 3), mar=c(8,8,5,5), mgp=MYMGP)
plot_correl(means.GNG$intercept_muc, means.TAB$intercept_muc, expression(alpha[mu] ~ "for FC (GNG)"), expression(alpha[mu] ~ "for FC (TAB)"), cex = MYCEX, asp=1)
plot_correl(means.TAB$intercept_muc, means.RS$intercept_muc, expression(alpha[mu] ~ "for FC (TAB)"), expression(alpha[mu] ~ "for FC (RS)"), cex = MYCEX, asp=1)
plot_correl(means.RS$intercept_muc, means.GNG$intercept_muc, expression(alpha[mu] ~ "for FC (RS)"), expression(alpha[mu] ~ "for FC (GNG)"), cex = MYCEX, asp=1)

plot_correl(means.GNG$log_intercept_sigmac, means.TAB$log_intercept_sigmac, expression(alpha[sigma] ~ "for FC (GNG)"), expression(alpha[sigma] ~ "for FC (TAB)"), cex = MYCEX, asp=1)
plot_correl(means.TAB$log_intercept_sigmac, means.RS$log_intercept_sigmac, expression(alpha[sigma] ~ "for FC (TAB)"), expression(alpha[sigma] ~ "for FC (RS)"), cex = MYCEX, asp=1)
plot_correl(means.RS$log_intercept_sigmac, means.GNG$log_intercept_sigmac, expression(alpha[sigma] ~ "for FC (RS)"), expression(alpha[sigma] ~ "for FC (GNG)"), cex = MYCEX, asp=1)


# REPRESENT PARAMETERS ----------------------
# Represent parameters in a matrix
slope_muc.GNG = means.GNG$slope_muc
slope_muc.GNG[zeros.GNG$slope_muc] = NA
log_slope_sigmac.GNG = means.GNG$log_slope_sigmac
log_slope_sigmac.GNG[zeros.GNG$log_slope_sigmac] = NA
slope_muc.GNG.adj = get_adj(slope_muc.GNG, MODULES_FILE_MNI, labels_file=LABELS_FILE)
log_slope_sigmac.GNG.adj = get_adj(log_slope_sigmac.GNG, MODULES_FILE_MNI, labels_file=LABELS_FILE)
if(doplot){
save_fig(figname="FigureSupp2a_GNG", res=CRES)
plot_adj(slope_muc.GNG.adj, MODULES_FILE_MNI) 
save_fig(figname="FigureSupp2b_GNG", res=CRES)
plot_adj(log_slope_sigmac.GNG.adj, MODULES_FILE_MNI, lim=-2) 
}

intercept_muc.GNG = means.GNG$intercept_muc
intercept_muc.GNG[zeros.GNG$intercept_muc] = NA
log_intercept_sigmac.GNG = means.GNG$log_intercept_sigmac
log_intercept_sigmac.GNG[zeros.GNG$log_intercept_sigmac] = NA

slope_muc.GNG.adj = get_adj(slope_muc.GNG, MODULES_FILE_20, labels_file=LABELS_FILE)
log_slope_sigmac.GNG.adj = get_adj(log_slope_sigmac.GNG, MODULES_FILE_20, labels_file=LABELS_FILE)
intercept_muc.GNG.adj = get_adj(intercept_muc.GNG, MODULES_FILE_20, labels_file=LABELS_FILE)
log_intercept_sigmac.GNG.adj = get_adj(log_intercept_sigmac.GNG, MODULES_FILE_20, labels_file=LABELS_FILE)

if(doplot){  
save_fig(figname="FigureSuppBM20_GNG_slope_mu", res=CRES)
plot_adj(slope_muc.GNG.adj, MODULES_FILE_20) 
save_fig(figname="FigureSuppBM20_GNG_slope_sigma", res=CRES)
plot_adj(log_slope_sigmac.GNG.adj, MODULES_FILE_20) 

save_fig(figname="FigureSuppBM20_GNG_intercept_mu", res=CRES)
plot_adj(intercept_muc.GNG.adj, MODULES_FILE_20) 
save_fig(figname="FigureSuppBM20_GNG_intercept_sigma", res=CRES)
plot_adj(log_intercept_sigmac.GNG.adj, MODULES_FILE_20) 
}

# Represent parameters by modules 
if(doplot){
save_fig(res=CRES)
plot_matrix(get_adj_modules(means.mod.20.GNG$intercept_muc, MODULES_FILE_20))
save_fig(figname="FigureBetaMu_20_GNG", res=CRES)
plot_matrix(get_adj_modules(means.mod.20.GNG$slope_muc, MODULES_FILE_20))
save_fig(res=CRES)
plot_matrix(get_adj_modules(exp(means.mod.20.GNG$log_intercept_sigmac), MODULES_FILE_20))
save_fig(figname="FigureBetaSigma_20_GNG", res=CRES)
plot_matrix(get_adj_modules(exp(means.mod.20.GNG$log_slope_sigmac*50), MODULES_FILE_20))
}

# average parameters over matrices 
intercept_muc.GNG.adj.20 = average.ICN.GNG$intercept_muc$means
intercept_muc.GNG.adj.20[average.ICN.GNG$intercept_muc$quantsl < 0 & average.ICN.GNG$intercept_muc$quantsh > 0] = 0
log_intercept_sigmac.GNG.adj.20 = average.ICN.GNG$log_intercept_sigmac$means
log_intercept_sigmac.GNG.adj.20[average.ICN.GNG$log_intercept_sigmac$quantsl < 0 & average.ICN.GNG$log_intercept_sigmac$quantsh > 0] = 0

if(doplot){
save_fig(figname="FigureSuppBM20_GNG_average_intercept_mu", res=CRES)
plot_matrix(intercept_muc.GNG.adj.20)
save_fig(figname="FigureSuppBM20_GNG_average_intercept_sigma", res=CRES)
plot_matrix(exp(log_intercept_sigmac.GNG.adj.20))
}

slope_muc.GNG.adj.20 = average.ICN.GNG$slope_muc$means
slope_muc.GNG.adj.20[average.ICN.GNG$slope_muc$quantsl < 0 & average.ICN.GNG$slope_muc$quantsh > 0] = 0
log_slope_sigmac.GNG.adj.20 = average.ICN.GNG$log_slope_sigmac$means
log_slope_sigmac.GNG.adj.20[average.ICN.GNG$log_slope_sigmac$quantsl < 0 & average.ICN.GNG$log_slope_sigmac$quantsh > 0] = 0

if(doplot){
save_fig(figname="FigureSuppBM20_GNG_average_slope_mu", res=CRES)
plot_matrix(slope_muc.GNG.adj.20)
save_fig(figname="FigureSuppBM20_GNG_average_slope_sigma", res=CRES)
plot_matrix(log_slope_sigmac.GNG.adj.20)
}

# plot intra correlations
dist.mat.orig = as.matrix(read.csv(DISTANCE_FILE, sep=' ', header=F))
dist.mat = order_vals(dist.mat.orig, MODULES_FILE_20, ismatrix=TRUE )
dist.ICN = average_params(dist.mat.orig, MODULES_FILE_20, ICN=TRUE)

save_fig(res=BWRES)
par(mar=c(8,8,5,5), mgp=MYMGP)
barplot(sort(diag(dist.ICN)), xlab = "Average distance (mm)", ylab="Network",  
        horiz = TRUE, cex.lab= CEX_LAB, cex.axis = CEX_AXIS, cex.names = CEX_NAMES, las = 1, 
        width = 1, space=0,  col="gray90") 


# plot intra correlations
save_fig(figname="Figure3_GNG_e", res=BWRES)
plot_intra_values.slopes(average.ICN.GNG.20, MODULES_FILE_20)
save_fig(res=BWRES)
plot_intra_values.intercept(average.ICN.GNG.20, MODULES_FILE_20)

save_fig(figname="FigureIntra", res=BWRES)
plot_intra_values.slopes(average.ICN.GNG.20, MODULES_FILE_20, plot_inter=F)


plot(diag(average.ICN.GNG$slope_muc$means), diag(dist.ICN), pch=20)
cor.test(diag(average.ICN.GNG$slope_muc$means), diag(dist.ICN))

# list strongest connections
list.mu.GNG.20 = list_connections(slope_muc.GNG.adj.20, n=NULL, incdiag=TRUE, absord=TRUE)
list.sigma.GNG.20 = list_connections(log_slope_sigmac.GNG.adj.20, n=NULL, incdiag=TRUE, absord=TRUE)
print(list.mu.GNG.20[1:10,])
print(list.sigma.GNG.20[1:10,])

# plot cognitive domains
save_fig(res=BWRES)
plot_cognitive_domains(average.cog.GNG.70, MODULES_FILE_70, bycluster=FALSE)
save_fig(figname="FigureSupp2c_GNG", res=BWRES)
plot_cognitive_domains(average.cog.GNG.70, MODULES_FILE_70, bycluster=TRUE)

save_fig(res=BWRES)
par(mfrow=c(3, 2), mar=c(8,8,5,5), mgp=MYMGP)
plot_correl(average.cog.TAB.70$slope_muc$means, average.cog.GNG.70$slope_muc$means, expression(beta[mu] ~ " cog for FC (TAB)"), expression(beta[mu] ~ " cog for FC (GNG)"), cex = MYCEX, asp=1)
plot_correl(diag(average.ICN.TAB.20$slope_muc$means), diag(average.ICN.GNG.20$slope_muc$means), expression(beta[mu] ~ " ICN for FC (TAB)"), expression(beta[mu] ~ " ICN for FC (GNG)"), cex = MYCEX, asp=1)
plot_correl(average.cog.RS.70$slope_muc$means, average.cog.GNG.70$slope_muc$means, expression(beta[mu] ~ "cog for FC (RS)"), expression(beta[mu] ~ "cog for FC (GNG)"), cex = MYCEX, asp=1)
plot_correl(diag(average.ICN.RS.20$slope_muc$means), diag(average.ICN.GNG.20$slope_muc$means), expression(beta[mu] ~ " ICN for FC (RS)"), expression(beta[mu] ~ " ICN for FC (GNG)"), cex = MYCEX, asp=1)
plot_correl(average.cog.RS.70$slope_muc$means, average.cog.TAB.70$slope_muc$means, expression(beta[mu] ~ "cog for FC (RS)"), expression(beta[mu] ~ "cog for FC (TAB)"), cex = MYCEX, asp=1)
plot_correl(diag(average.ICN.RS.20$slope_muc$means), diag(average.ICN.TAB.20$slope_muc$means), expression(beta[mu] ~ " ICN for FC (RS)"), expression(beta[mu] ~ " ICN for FC (TAB)"), cex = MYCEX, asp=1)

save_fig(res=BWRES)
par(mfrow=c(3, 2), mar=c(8,8,5,5), mgp=MYMGP)
plot_correl(average.cog.TAB.70$log_slope_sigmac$means, average.cog.GNG.70$log_slope_sigmac$means, expression(beta[sigma] ~ " cog for FC (TAB)"), expression(beta[sigma] ~ " cog for FC (GNG)"), cex = MYCEX, asp=1)
plot_correl(diag(average.ICN.TAB.20$log_slope_sigmac$means), diag(average.ICN.GNG.20$log_slope_sigmac$means), expression(beta[sigma] ~ " ICN for FC (TAB)"), expression(beta[sigma] ~ " ICN for FC (GNG)"), cex = MYCEX, asp=1)
plot_correl(average.cog.RS.70$log_slope_sigmac$means, average.cog.GNG.70$log_slope_sigmac$means, expression(beta[sigma] ~ "cog for FC (RS)"), expression(beta[sigma] ~ "cog for FC (GNG)"), cex = MYCEX, asp=1)
plot_correl(diag(average.ICN.RS.20$log_slope_sigmac$means), diag(average.ICN.GNG.20$log_slope_sigmac$means), expression(beta[sigma] ~ " ICN for FC (RS)"), expression(beta[sigma] ~ " ICN for FC (GNG)"), cex = MYCEX, asp=1)
plot_correl(average.cog.RS.70$log_slope_sigmac$means, average.cog.TAB.70$log_slope_sigmac$means, expression(beta[sigma] ~ "cog for FC (RS)"), expression(beta[sigma] ~ "cog for FC (TAB)"), cex = MYCEX, asp=1)
plot_correl(diag(average.ICN.RS.20$log_slope_sigmac$means), diag(average.ICN.TAB.20$log_slope_sigmac$means), expression(beta[sigma] ~ " ICN for FC (RS)"), expression(beta[sigma] ~ " ICN for FC (TAB)"), cex = MYCEX, asp=1)

# define parameters again
intercept_muc.GNG = means.GNG$intercept_muc
log_intercept_sigmac.GNG = means.GNG$log_intercept_sigmac
slope_muc.GNG = means.GNG$slope_muc
log_slope_sigmac.GNG = means.GNG$log_slope_sigmac

intercept_muc.TAB = means.TAB$intercept_muc
log_intercept_sigmac.TAB = means.TAB$log_intercept_sigmac
slope_muc.TAB = means.TAB$slope_muc
log_slope_sigmac.TAB = means.TAB$log_slope_sigmac

intercept_muc.RS = means.RS$intercept_muc
log_intercept_sigmac.RS = means.RS$log_intercept_sigmac
slope_muc.RS = means.RS$slope_muc
log_slope_sigmac.RS = means.RS$log_slope_sigmac

# just looking at baseline - positive FC!!!
#slope_muc.GNG[!(intercept_muc.GNG > 0)] = NA
#intercept_muc.GNG[!(intercept_muc.GNG > 0)] = NA

# DEFINE PARAMETERS FOR NETWORKS ----------------------

intercept_muc.GNG.adj = get_adj(intercept_muc.GNG, MODULES_FILE_20, sort=TRUE, labels_file=LABELS_FILE)
log_intercept_sigmac.GNG.adj = get_adj(log_intercept_sigmac.GNG, MODULES_FILE_20, sort=TRUE, labels_file=LABELS_FILE)
slope_muc.GNG.adj = get_adj(slope_muc.GNG, MODULES_FILE_20, sort=TRUE, labels_file=LABELS_FILE)
log_slope_sigmac.GNG.adj = get_adj(log_slope_sigmac.GNG, MODULES_FILE_20, sort=TRUE, labels_file=LABELS_FILE)

intercept_muc.TAB.adj = get_adj(intercept_muc.TAB, MODULES_FILE_20, sort=TRUE, labels_file=LABELS_FILE)
log_intercept_sigmac.TAB.adj = get_adj(log_intercept_sigmac.TAB, MODULES_FILE_20, sort=TRUE, labels_file=LABELS_FILE)
slope_muc.TAB.adj = get_adj(slope_muc.TAB, MODULES_FILE_20, sort=TRUE, labels_file=LABELS_FILE)
log_slope_sigmac.TAB.adj = get_adj(log_slope_sigmac.TAB, MODULES_FILE_20, sort=TRUE, labels_file=LABELS_FILE)

intercept_muc.RS.adj = get_adj(intercept_muc.RS, MODULES_FILE_20, sort=TRUE, labels_file=LABELS_FILE)
log_intercept_sigmac.RS.adj = get_adj(log_intercept_sigmac.RS, MODULES_FILE_20, sort=TRUE, labels_file=LABELS_FILE)
slope_muc.RS.adj = get_adj(slope_muc.RS, MODULES_FILE_20, sort=TRUE, labels_file=LABELS_FILE)
log_slope_sigmac.RS.adj = get_adj(log_slope_sigmac.RS, MODULES_FILE_20, sort=TRUE, labels_file=LABELS_FILE)

intercept_muc.GNG.zeros.adj = get_adj(zeros.GNG$intercept_muc, MODULES_FILE_20, sort=TRUE, labels_file=LABELS_FILE)

intercept_muc.GNG.degree = compute_degree(intercept_muc.GNG.adj)
log_intercept_sigmac.GNG.degree = compute_degree(log_intercept_sigmac.GNG.adj)
slope_muc.GNG.degree = compute_degree(slope_muc.GNG.adj)
log_slope_sigmac.GNG.degree = compute_degree(log_slope_sigmac.GNG.adj)

intercept_muc.TAB.degree = compute_degree(intercept_muc.TAB.adj)
log_intercept_sigmac.TAB.degree = compute_degree(log_intercept_sigmac.TAB.adj)
slope_muc.TAB.degree = compute_degree(slope_muc.TAB.adj)
log_slope_sigmac.TAB.degree = compute_degree(log_slope_sigmac.TAB.adj)

intercept_muc.RS.degree = compute_degree(intercept_muc.RS.adj)
log_intercept_sigmac.RS.degree = compute_degree(log_intercept_sigmac.RS.adj)
slope_muc.RS.degree = compute_degree(slope_muc.RS.adj)
log_slope_sigmac.RS.degree = compute_degree(log_slope_sigmac.RS.adj)

#slope_muc.TAB.degree = compute_degree(slope_muc.TAB.adj)
#log_slope_sigmac.TAB.degree = compute_degree(log_slope_sigmac.TAB.adj)

slope_muc.PET = order_vals(means.PET$slope_muc, MODULES_FILE_20)
slope_muc.VBM = order_vals(means.VBM$slope_muc, MODULES_FILE_20)
intercept_muc.PET = order_vals(means.PET$intercept_muc, MODULES_FILE_20)
intercept_muc.VBM = order_vals(means.VBM$intercept_muc, MODULES_FILE_20)

log_intercept_sigmac.PET = order_vals(means.PET$log_intercept_sigmac, MODULES_FILE_20)
log_intercept_sigmac.VBM = order_vals(means.VBM$log_intercept_sigmac, MODULES_FILE_20)

log_slope_sigmac.PET = order_vals(means.PET$log_slope_sigmac, MODULES_FILE_20)
log_slope_sigmac.VBM = order_vals(means.VBM$log_slope_sigmac, MODULES_FILE_20)

# CORRELATIONS BETWEEN PET AND VBM PARAMETERS
all_params.PET = data.frame(alpha_mu=intercept_muc.PET, beta_mu=slope_muc.PET, alpha_sigma=log_intercept_sigmac.PET, beta_sigma=log_slope_sigmac.PET ) 
all_params.PET = subset(all_params.PET, !is.na(slope_muc.PET))
all_params.VBM = data.frame(alpha_mu=means.VBM$intercept_muc, beta_mu=means.VBM$slope_muc, alpha_sigma=means.VBM$log_intercept_sigmac, beta_sigma=means.VBM$log_slope_sigmac ) 
cor(all_params.PET)
cor(all_params.VBM)

# ASSOCIATIONS OF FC WITH PET AND VBM

#PET
save_fig(figname="Figure3_GNG_f_BP", res=BWRES)
par(mar=c(8,9,5,5), mgp=MYMGP)
plot_correl(slope_muc.GNG.degree, slope_muc.PET, expression("Average nodal " ~ beta[mu] ~ "for FC"), expression(beta[mu] ~ "for BP"), !is.na(slope_muc.GNG.degree * slope_muc.PET) )

save_fig(figname="Figure3_GNG_f_BP_sigma", res=BWRES)
#plot_correl(log_slope_sigmac.GNG.degree, slope_muc.PET, expression("Average nodal " ~ beta[sigma] ~ "for FC"), expression(beta[mu] ~ "for BP"), !is.na(log_slope_sigmac.GNG.degree))
par(mar=c(8,9,5,5), mgp=MYMGP)
plot_correl(log_slope_sigmac.GNG.degree, log_slope_sigmac.PET, expression("Average nodal " ~ beta[sigma] ~ "for FC"), expression(beta[sigma] ~ "for BP"), !is.na(log_slope_sigmac.GNG.degree* log_slope_sigmac.PET))

# which regins have largest dopamine loss ?
print(labels[means.PET$slope_muc < -0.007, ])

# is the correlation still there if I remove those regions
slope_muc.PET.subset  = slope_muc.PET
slope_muc.PET.subset[slope_muc.PET < - 0.007] = NA

save_fig(res=BWRES)
par(mar=c(8,9,5,5), mgp=MYMGP, mfrow=c(1,3))
plot_correl(slope_muc.GNG.degree, slope_muc.PET, expression("Average nodal " ~ beta[mu] ~ "for FC"), expression(beta[mu] ~ "for BP"), !is.na(slope_muc.GNG.degree * slope_muc.PET.subset))
plot_correl(slope_muc.TAB.degree, slope_muc.PET, expression("Average nodal " ~ beta[mu] ~ "for FC"), expression(beta[mu] ~ "for BP"), !is.na(slope_muc.GNG.degree * slope_muc.PET.subset))
plot_correl(slope_muc.RS.degree, slope_muc.PET, expression("Average nodal " ~ beta[mu] ~ "for FC"), expression(beta[mu] ~ "for BP"), !is.na(slope_muc.GNG.degree * slope_muc.PET.subset))

#VBM
save_fig(figname="Figure3_GNG_f_VBM", res=BWRES)
par(mar=c(8,9,5,5), mgp=MYMGP)
plot_correl(slope_muc.GNG.degree, slope_muc.VBM, expression("Average nodal " ~ beta[mu] ~ "for FC"), expression(beta[mu] ~ "for GMV"), !is.na(slope_muc.GNG.degree), right = T)
#plot_correl(log_slope_sigmac.GNG.degree, slope_muc.VBM, expression("Average nodal " ~ beta[sigma] ~ "for FC"), expression(beta[mu] ~ "for GMV"), !is.na(log_slope_sigmac.GNG.degree))
save_fig(figname="Figure3_GNG_f_VBM_sigma", res=BWRES)
par(mar=c(8,9,5,5), mgp=MYMGP)
plot_correl(log_slope_sigmac.GNG.degree, log_slope_sigmac.VBM, expression("Average nodal " ~ beta[sigma] ~ "for FC"), expression(beta[sigma] ~ "for GMV"), !is.na(log_slope_sigmac.GNG.degree))

save_fig(figname="Figure3_PET_VBM_GNG", res=BWRES)
par(mfrow =c(2,2), mar=c(8,9,5,5), mgp=MYMGP)
plot_correl(slope_muc.GNG.degree, slope_muc.PET, expression("Average nodal " ~ beta[mu] ~ "for FC"), expression(beta[mu] ~ "for BP"), !is.na(slope_muc.GNG.degree * slope_muc.PET) )
plot_correl(log_slope_sigmac.GNG.degree, log_slope_sigmac.PET, expression("Average nodal " ~ beta[sigma] ~ "for FC"), expression(beta[sigma] ~ "for BP"), !is.na(log_slope_sigmac.GNG.degree* log_slope_sigmac.PET))
plot_correl(slope_muc.GNG.degree, slope_muc.VBM, expression("Average nodal " ~ beta[mu] ~ "for FC"), expression(beta[mu] ~ "for GMV"), !is.na(slope_muc.GNG.degree))
plot_correl(log_slope_sigmac.GNG.degree, log_slope_sigmac.VBM, expression("Average nodal " ~ beta[sigma] ~ "for FC"), expression(beta[sigma] ~ "for GMV"), !is.na(log_slope_sigmac.GNG.degree))

save_fig(figname="Figure3_PET_VBM_TAB", res=BWRES)
par(mfrow =c(2,2), mar=c(8,9,5,5), mgp=MYMGP)
plot_correl(slope_muc.TAB.degree, slope_muc.PET, expression("Average nodal " ~ beta[mu] ~ "for FC"), expression(beta[mu] ~ "for BP"), !is.na(slope_muc.TAB.degree * slope_muc.PET) )
plot_correl(log_slope_sigmac.TAB.degree, log_slope_sigmac.PET, expression("Average nodal " ~ beta[sigma] ~ "for FC"), expression(beta[sigma] ~ "for BP"), !is.na(log_slope_sigmac.TAB.degree* log_slope_sigmac.PET))
plot_correl(slope_muc.TAB.degree, slope_muc.VBM, expression("Average nodal " ~ beta[mu] ~ "for FC"), expression(beta[mu] ~ "for GMV"), !is.na(slope_muc.TAB.degree))
plot_correl(log_slope_sigmac.TAB.degree, log_slope_sigmac.VBM, expression("Average nodal " ~ beta[sigma] ~ "for FC"), expression(beta[sigma] ~ "for GMV"), !is.na(log_slope_sigmac.TAB.degree))

save_fig(figname="Figure3_PET_VBM_RS", res=BWRES)
par(mfrow =c(2,2), mar=c(8,9,5,5), mgp=MYMGP)
plot_correl(slope_muc.RS.degree, slope_muc.PET, expression("Average nodal " ~ beta[mu] ~ "for FC"), expression(beta[mu] ~ "for BP"), !is.na(slope_muc.RS.degree * slope_muc.PET) )
plot_correl(log_slope_sigmac.RS.degree, log_slope_sigmac.PET, expression("Average nodal " ~ beta[sigma] ~ "for FC"), expression(beta[sigma] ~ "for BP"), !is.na(log_slope_sigmac.RS.degree* log_slope_sigmac.PET))
plot_correl(slope_muc.RS.degree, slope_muc.VBM, expression("Average nodal " ~ beta[mu] ~ "for FC"), expression(beta[mu] ~ "for GMV"), !is.na(slope_muc.RS.degree))
plot_correl(log_slope_sigmac.RS.degree, log_slope_sigmac.VBM, expression("Average nodal " ~ beta[sigma] ~ "for FC"), expression(beta[sigma] ~ "for GMV"), !is.na(log_slope_sigmac.RS.degree))


# PET controlling for VBM
summary(lm(slope_muc.GNG.degree ~ scale(slope_muc.PET) + scale(slope_muc.VBM)))
summary(lm(slope_muc.TAB.degree ~ scale(slope_muc.PET) + scale(slope_muc.VBM)))
summary(lm(slope_muc.RS.degree ~ scale(slope_muc.PET) + scale(slope_muc.VBM)))

summary(lm(slope_muc.GNG.degree ~ scale(slope_muc.PET.subset) + scale(slope_muc.VBM)))
summary(lm(slope_muc.TAB.degree ~ scale(slope_muc.PET.subset) + scale(slope_muc.VBM)))
summary(lm(slope_muc.RS.degree ~ scale(slope_muc.PET.subset) + scale(slope_muc.VBM)))


summary(lm(log_slope_sigmac.GNG.degree ~ scale(log_slope_sigmac.PET) + scale(log_slope_sigmac.VBM)))
summary(lm(log_slope_sigmac.TAB.degree ~ scale(log_slope_sigmac.PET) + scale(log_slope_sigmac.VBM)))
summary(lm(log_slope_sigmac.RS.degree ~ scale(log_slope_sigmac.PET) + scale(log_slope_sigmac.VBM)))


# PET controlling for intercept of PET
summary(lm(slope_muc.GNG.degree ~ slope_muc.PET + intercept_muc.PET))

# REPRESENT SIMILARITIES
save_fig(figname="Figure1b", res=CRES)
plot_connectome_similarities(INPUT_FILE.GNG, demo)
save_fig()
plot_connectome_similarities(INPUT_FILE.TAB, demo)
save_fig()
plot_connectome_similarities(INPUT_FILE.RS, demo)

# MDS plots
save_fig(figname="Figure1c_GNG", res=CRES)
plot_connectome_similarities(INPUT_FILE.GNG, demo, plot_mds=TRUE)
save_fig()
plot_connectome_similarities(INPUT_FILE.TAB, demo, plot_mds=TRUE)
save_fig()
plot_connectome_similarities(INPUT_FILE.RS, demo, plot_mds=TRUE)


#COGNITIVE SCORES

# plot them
#plot(updating ~ age, data=demo, xlab = "Age", ylab = "Score", pch=20, cex.axis=CEX_AXIS, cex.lab=CEX_LAB, cex=2)

# create a WM score
estimates = p.values = coefs.Estimate = coefs.p = NULL
summary(demo[, c("updating", "n_back")])
outliers = c("D84")

demo.old = demo[demo$class=='Old', ]
pca = princomp(demo.old[!demo.old$Subject %in% outliers, c("updating", "n_back")], cor = T)

var = pca$sdev**2/sum(pca$sdev**2)
demo.clean = demo
demo.clean$combined = 0
demo.clean$combined[!demo.clean$Subject %in% outliers & demo.clean$class=='Old'] = pca$scores[,"Comp.1"]

# missing data
demo.clean[demo.clean$Subject %in% outliers, c("combined", "n_back")] = NA 
md.pattern(demo.clean[c("combined", "updating", "n_back", "age", "sex")])

# impute the combined scored base on what we have
m_imputations = 10
demo.imputed = mice(demo.clean[, c("combined", "updating", "n_back", "age", "sex")], maxit = 50, m = m_imputations, method = 'pmm', seed = 500) # impute by mean

mycombined = 0

for (m in seq(m_imputations)){  # repeat for the number of imputations
  demo.clean = demo
  demo.clean[c("combined", "updating", "n_back", "age", "sex")] = complete(demo.imputed, m)
  mycombined = demo.clean$combined + mycombined
  # model of temporal connections vs updating or nback
  print("Cognitive associations")
  test = "combined"


# ASSOCIATION SCORES AND FC
# 20 
similarity.score.mod.20.GNG = normal_pattern_sim(INPUT_FILE.GNG.mod.20, demo.clean, beta_mu= means.mod.20.GNG$slope_muc, 
                                                           beta_sigma = means.mod.20.GNG$log_slope_sigmac, 
                                                           test, modules_file = MODULES_FILE_20)

similarity.score.mod.20.TAB = normal_pattern_sim(INPUT_FILE.TAB.mod.20, demo.clean, beta_mu= means.mod.20.TAB$slope_muc, 
                                                beta_sigma = means.mod.20.TAB$log_slope_sigmac, 
                                                test, modules_file = MODULES_FILE_20)

similarity.score.mod.20.RS = normal_pattern_sim(INPUT_FILE.RS.mod.20, demo.clean, beta_mu= means.mod.20.RS$slope_muc, 
                                                beta_sigma = means.mod.20.RS$log_slope_sigmac, 
                                                test, modules_file = MODULES_FILE_20)


similarity.score.mod.70.GNG = normal_pattern_sim(INPUT_FILE.GNG.mod.70, demo.clean, beta_mu= means.mod.70.GNG$slope_muc, 
                                                 beta_sigma = means.mod.70.GNG$log_slope_sigmac, 
                                                 test, modules_file = MODULES_FILE_70)

similarity.score.mod.70.TAB = normal_pattern_sim(INPUT_FILE.TAB.mod.70, demo.clean, beta_mu= means.mod.70.TAB$slope_muc, 
                                                 beta_sigma = means.mod.70.TAB$log_slope_sigmac, 
                                                 test, modules_file = MODULES_FILE_70)

similarity.score.mod.70.RS = normal_pattern_sim(INPUT_FILE.RS.mod.70, demo.clean, beta_mu= means.mod.70.RS$slope_muc, 
                                                beta_sigma = means.mod.70.RS$log_slope_sigmac, 
                                                test, modules_file = MODULES_FILE_70)

# gather estimates
estimates = rbind(estimates, 
                  c(similarity.score.mod.20.GNG$estimate, 
                    similarity.score.mod.20.TAB$estimate,
                    similarity.score.mod.20.RS$estimate,
                    similarity.score.mod.70.GNG$estimate, 
                    similarity.score.mod.70.TAB$estimate,
                    similarity.score.mod.70.RS$estimate)
)
p.values  = rbind(p.values, 
                  c(similarity.score.mod.20.GNG$p.value, 
                    similarity.score.mod.20.TAB$p.value,
                    similarity.score.mod.20.RS$p.value,
                    similarity.score.mod.70.GNG$p.value, 
                    similarity.score.mod.70.TAB$p.value,
                    similarity.score.mod.70.RS$p.value)
)



association.score.mod.20.GNG = associate_cognition_modules(INPUT_FILE.GNG.mod.20, demo.clean, beta_mu= means.mod.20.GNG$slope_muc, 
                                                           beta_sigma = means.mod.20.GNG$log_slope_sigmac, 
                                                           test, modules_file = MODULES_FILE_20)
coefs.20.GNG = model.cognitive.association(association.score.mod.20.GNG, means.mod.20.GNG$slope_muc, means.mod.20.GNG$log_slope_sigmac, 
                            figname="Figure4_GNG_b_modules20", means.mod.20.GNG$intercept_muc, means.mod.20.GNG$log_intercept_sigmac)

association.score.mod.20.TAB = associate_cognition_modules(INPUT_FILE.TAB.mod.20, demo.clean, beta_mu= means.mod.20.TAB$slope_muc, 
                                                           beta_sigma = means.mod.20.TAB$log_slope_sigmac,
                                                           test, modules_file = MODULES_FILE_20)
coefs.20.TAB = model.cognitive.association(association.score.mod.20.TAB, means.mod.20.TAB$slope_muc, means.mod.20.TAB$log_slope_sigmac, 
                            figname="Figure4_TAB_b_modules", means.mod.20.TAB$intercept_muc, means.mod.20.TAB$log_intercept_sigmac)

association.score.mod.20.RS = associate_cognition_modules(INPUT_FILE.RS.mod.20, demo.clean, test, beta_mu= means.mod.20.RS$slope_muc, 
                                                          beta_sigma = means.mod.20.RS$log_slope_sigmac,
                                                          modules_file = MODULES_FILE_20)
coefs.20.RS = model.cognitive.association(association.score.mod.20.RS, means.mod.20.RS$slope_muc, means.mod.20.RS$log_slope_sigmac, 
                            figname="Figure4_RS_b_modules", means.mod.20.RS$intercept_muc, means.mod.20.RS$log_intercept_sigmac)


# 70 
association.score.mod.70.GNG = associate_cognition_modules(INPUT_FILE.GNG.mod.70, demo.clean, beta_mu= means.mod.70.GNG$slope_muc, 
                                                           beta_sigma = means.mod.70.GNG$log_slope_sigmac, 
                                                           test, modules_file = MODULES_FILE_70)

coefs.70.GNG = model.cognitive.association(association.score.mod.70.GNG, means.mod.70.GNG$slope_muc, means.mod.70.GNG$log_slope_sigmac, 
                            figname="Figure4_GNG_b_modules70", means.mod.70.GNG$intercept_muc, means.mod.70.GNG$log_intercept_sigmac)

association.score.mod.70.TAB = associate_cognition_modules(INPUT_FILE.TAB.mod.70, demo.clean, beta_mu= means.mod.70.TAB$slope_muc, 
                                                           beta_sigma = means.mod.70.TAB$log_slope_sigmac,
                                                           test, modules_file = MODULES_FILE_70)

coefs.70.TAB = model.cognitive.association(association.score.mod.70.TAB, means.mod.70.TAB$slope_muc, means.mod.70.TAB$log_slope_sigmac, 
                            figname="Figure4_TAB_b_modules", means.mod.70.TAB$intercept_muc, means.mod.70.TAB$log_intercept_sigmac)

association.score.mod.70.RS = associate_cognition_modules(INPUT_FILE.RS.mod.70, demo.clean, test, beta_mu= means.mod.70.RS$slope_muc, 
                                                          beta_sigma = means.mod.70.RS$log_slope_sigmac,
                                                          modules_file = MODULES_FILE_70)

coefs.70.RS = model.cognitive.association(association.score.mod.70.RS, means.mod.70.RS$slope_muc, means.mod.70.RS$log_slope_sigmac, 
                            figname="Figure4_RS_b_modules", means.mod.70.RS$intercept_muc, means.mod.70.RS$log_intercept_sigmac)


coefs.Estimate = abind(coefs.Estimate, 
                  cbind(coefs.20.GNG[c("beta_mu", "beta_sigma"), "Estimate"], 
                    coefs.20.TAB[c("beta_mu", "beta_sigma"), "Estimate"],
                    coefs.20.RS[c("beta_mu", "beta_sigma"), "Estimate"],
                    coefs.70.GNG[c("beta_mu", "beta_sigma"), "Estimate"], 
                    coefs.70.TAB[c("beta_mu", "beta_sigma"), "Estimate"],
                    coefs.70.RS[c("beta_mu", "beta_sigma"), "Estimate"]), along=3)


coefs.p = abind(coefs.p, 
                   cbind(coefs.20.GNG[c("beta_mu", "beta_sigma"), "Pr(>|t|)"], 
                         coefs.20.TAB[c("beta_mu", "beta_sigma"), "Pr(>|t|)"],
                         coefs.20.RS[c("beta_mu", "beta_sigma"), "Pr(>|t|)"],
                         coefs.70.GNG[c("beta_mu", "beta_sigma"), "Pr(>|t|)"], 
                         coefs.70.TAB[c("beta_mu", "beta_sigma"), "Pr(>|t|)"],
                         coefs.70.RS[c("beta_mu", "beta_sigma"), "Pr(>|t|)"]), along = 3)


} # end imputation


demo.clean$combined = mycombined/m_imputations

colnames(coefs.Estimate) = colnames(coefs.p) = 
  colnames(estimates) = colnames(p.values)  = c("GNG.20", "TAB.20", "RS.20", "GNG.70", "TAB.70", "RS.70")
options(digits = 3)
print(apply(estimates, 2, mean))
print(apply(p.values, 2, lichtrubin))

print(apply(coefs.Estimate, c(1, 2), mean))
print(apply(coefs.p,  c(1, 2), lichtrubin))


# 
association.score.mod.20.GNG = associate_cognition_modules(INPUT_FILE.GNG.mod.20, demo.clean, beta_mu= means.mod.20.GNG$slope_muc, 
                                                           beta_sigma = means.mod.20.GNG$log_slope_sigmac, 
                                                           test, modules_file = MODULES_FILE_20)
coefs.20.GNG = model.cognitive.association(association.score.mod.20.GNG, means.mod.20.GNG$slope_muc, means.mod.20.GNG$log_slope_sigmac, 
                                           figname="Figure4_GNG_b_modules20", means.mod.20.GNG$intercept_muc, means.mod.20.GNG$log_intercept_sigmac,
                                           plotme = T)

association.score.mod.70.GNG = associate_cognition_modules(INPUT_FILE.GNG.mod.70, demo.clean, beta_mu= means.mod.70.GNG$slope_muc, 
                                                           beta_sigma = means.mod.70.GNG$log_slope_sigmac, 
                                                           test, modules_file = MODULES_FILE_70)

coefs.70.GNG = model.cognitive.association(association.score.mod.70.GNG, means.mod.70.GNG$slope_muc, means.mod.70.GNG$log_slope_sigmac, 
                                           figname="Figure4_GNG_b_modules70", means.mod.70.GNG$intercept_muc, means.mod.70.GNG$log_intercept_sigmac,
                                           plotme =T) 


# DO PLOTS WITH AVERAGE VALUES

save_fig(figname="Figure4a", res=BWRES, height = 6.5)
par(mar=c(8,8,5,5), mgp = c(5, 2, 0))

plot(n_back ~ updating, data=demo, xlab = "Letter updating score", ylab = "3-back score", pch=c(1, 2)[class], 
     cex.axis=CEX_AXIS + .5, cex.lab=CEX_LAB + 1, cex=6, lwd = 11, bg="grey70")
points(n_back ~ updating, data=subset(demo, Subject %in% outliers), pch=4, cex=6, lwd = 10)
arrows(x0=pca$center["updating"], 
       y0=pca$center["n_back"],  
       x1=pca$center["updating"] + loadings(pca)["updating","Comp.1"]*pca$scale[1], 
       y1=pca$center["n_back"] + loadings(pca)["n_back","Comp.1"]*pca$scale[2], code=2, cex=10, lwd=20, col = 'grey40', length = .8)

legend("bottomleft", legend=c("Younger", "Older"), pch =c(2, 1), cex = CEX_LAB + 1, lwd = 11, lty = 0)
dev.off()

# CLUSTERING CONNECTIONS
#GNG
clusters = cluster_params_2points(FC.GNG.mu.20$mean, FC.GNG.mu.70$mean, means.GNG$log_slope_sigmac, FC.GNG.mu.20$zeros, FC.GNG.mu.70$zeros, zeros.GNG$log_slope_sigmac, pal= CLUSPAL9, "Figure3_GNG-")

clusters.adj = get_adj(clusters, MODULES_FILE_MNI)
clusters.adj[clusters.adj==0] = NA
save_fig(figname = "Figure3_GNG_c", res=CRES)
plot_adj(clusters.adj, MODULES_FILE_MNI, lim = -1, pal=CLUSPAL18, nolegend=TRUE) 
save_fig(figname="Figure3_GNG_d", res=CRES)
plot_pies(clusters.adj, MODULES_FILE_MNI, pal=CLUSPAL18) 
dev.off()
clusters.adj = get_adj(clusters, MODULES_FILE_MNI_RED)
clusters.adj[clusters.adj==0] = NA
save_fig(res=CRES)
plot_adj(clusters.adj, MODULES_FILE_MNI_RED, lim = -1, pal=CLUSPAL18) 

#TAB
clusters = cluster_params_2points(FC.TAB.mu.20$mean, FC.TAB.mu.70$mean, means.TAB$log_slope_sigmac, FC.TAB.mu.20$zeros, FC.TAB.mu.70$zeros, zeros.TAB$log_slope_sigmac, pal= CLUSPAL9, "Figure3_TAB-")
clusters.adj = get_adj(clusters, MODULES_FILE_MNI)
clusters.adj[clusters.adj==0] = NA
save_fig(figname = "Figure3c_TAB", res=CRES)
plot_adj(clusters.adj, MODULES_FILE_MNI, lim = -1, pal=CLUSPAL18, nolegend=TRUE) 
save_fig(figname="Figure3d_TAB", res=CRES)
plot_pies(clusters.adj, MODULES_FILE_MNI, pal=CLUSPAL18) 

#RS
clusters = cluster_params_2points(FC.RS.mu.20$mean, FC.RS.mu.70$mean, means.RS$log_slope_sigmac, FC.RS.mu.20$zeros, FC.RS.mu.70$zeros, zeros.RS$log_slope_sigmac, pal= CLUSPAL9, "Figure3_RS-")
clusters.adj = get_adj(clusters, MODULES_FILE_MNI)
clusters.adj[clusters.adj==0] = NA
save_fig(figname = "Figure3c_RS", res=CRES)
if(doplot) plot_adj(clusters.adj, MODULES_FILE_MNI, lim = -1, pal=CLUSPAL18, nolegend=TRUE) 
save_fig(figname="Figure3d_RS", res=CRES)
if(doplot) plot_pies(clusters.adj, MODULES_FILE_MNI, pal=CLUSPAL18) 


# connection removal

conn.order.mu = order(abs(means.GNG$slope_muc), decreasing=TRUE)
conn.order.sigma = order(means.GNG$log_slope_sigmac, decreasing=TRUE)
if(doplot)
plot_connectome_removal(INPUT_FILE.GNG, demo, conn.order.mu, conn.order.sigma, abs(means.GNG$slope_muc), abs(means.GNG$log_slope_sigmac), figname = "Figure2_GNG-")

dev.off()



