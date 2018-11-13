rm(list=ls())
library(rstan)
library(R.matlab)
library(ggplot2)
library(doParallel)
library(foreach)

rstan_options(auto_write = TRUE)

# args <- commandArgs(TRUE)
# input_file = args[1]
# output_file = args[2]
# stan_file = args[3]
# priors_file = args[4]
# fitFC = args[5]
# NCLUSTERS = args[6]
# 
# print(paste('input_file =', input_file))
# print(paste('output_file =', output_file))
# print(paste('stan_file =', stan_file ))
# print(paste('priors_file =', priors_file ))
# print(paste('fitFC =', fitFC))
# print(paste('NCLUSTERS =', NCLUSTERS))

fit_models_connectome = function(input_file, output_file, stan_file, priors_file, fitFC, NCLUSTERS){

print(output_file)

priors = as.matrix(read.table(file = priors_file, header = TRUE, sep='\t'))
print(priors)

data_file = '~/Data/DAD/behaviour/demo_data.csv'
data = read.table(data_file, header=TRUE)

if (fitFC == 1) {
Conn_info = readMat(input_file)
labels = unlist(Conn_info$labels)

subjects_FC.1 = unlist(Conn_info$subjects)
X = data.frame(Conn_info$merged.matrices.mat)

vars = colnames(X)
n_conn = length(vars)
X$Subject = subjects_FC.1

X = merge(X, data, by="Subject")
values = X[, vars]

} else {

values = data.frame(read.table(input_file))
ncols = ncol(values)

names(values)[ncols] = "Subject"

vars = names(values)[1:(ncols-1)]
X = merge(values, data, by="Subject")
values = X[, vars]

}

print(dim(values))

#values = values[, 1:200]

# fit stan model
age = X$age
group = X$class

n_subjects = length(age)

age.shifted = age - 20
age.centered = scale(age, center=T, scale=F)

##################################
# parallel function
##################################

fit_parallel = function (y, x, inits, priors, mymodel, parnames){
  
  stan_data = list(S = length(x),
                   values = values[, k],
                   age = as.numeric(x),
                   priors = priors)
  
  fit <- rstan::optimizing(mymodel, data = stan_data, verbose = T, init = inits)
  return(fit$par[parnames])
  
}

##################################
##################################
##################################

parnames = c('intercept_muc', 
             'slope_muc', 
             'log_intercept_sigmac', 
             'log_slope_sigmac')

inits = list(intercept_muc = 0,
             slope_muc = 0,
             log_intercept_sigmac = 0,
             log_slope_sigmac = 0)


mymodel = stan_model(stan_file)

cl <- makeCluster(NCLUSTERS)
registerDoParallel(cl)

print(cl)

mypars = foreach(k = seq(ncol(values)), 
                  .packages = c('rstan'),
                  .export = c('fit_parallel'), 
                  .combine = 'rbind') %dopar% fit_parallel(values[, k], 
                                                           age.shifted, 
                                                           inits, 
                                                           priors, 
                                                           mymodel, 
                                                           parnames)

stopCluster(cl)

colnames(mypars) <- parnames
write.csv(x = mypars, file = output_file, row.names = F)

}