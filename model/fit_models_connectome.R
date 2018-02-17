rm(list=ls())
library(rstan)
library(R.matlab)
library(ggplot2)

print(tempdir())

args <- commandArgs(TRUE)
input_file = args[1]
samples_file = args[2]
stan_file = args[3]
priors_file = args[4]
fitFC = args[5]
USEMCMC = args[6]

print(paste('input_file =', input_file))
print(paste('samples_file =', samples_file))
print(paste('stan_file =', stan_file ))
print(paste('priors_file =', priors_file ))
print(paste('fitFC =', fitFC))
print(paste('USEMCMC =', USEMCMC))

filesep='/'

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

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

print(values)

#values = values[, 1:200]

# fit stan model
age = X$age
group = X$class

n_subjects = length(age)

age.shifted = age - 20
age.centered = scale(age, center=T, scale=F)

stan_data = list(S = length(age.shifted),
                 C = dim(values)[2],
                 values = values,
                 age = as.numeric(age.shifted),
                 priors = priors)

WARMUP=35000 
ITERS=WARMUP + 500

OUTPUT_SAMPLES=1000
VBITERS=200000


if (USEMCMC==1) {
# USE MCMC
#  stan_control = list(adapt_delta=0.99,  stepsize = 0.01, max_treedepth = 15)
  stan_control = list(adapt_delta=0.99, stepsize = 0.01, max_treedepth = 15)
  
  fit <- stan(file = stan_file, data = stan_data, iter = ITERS, chains = 4, pars = c("muc", "sigmac"), 
  include = FALSE, warmup = WARMUP, sample_file=samples_file, control = stan_control, thin = 2)
  #save(fit.FC, file="fit.FC.Rdata")
} else {
# USE ADVI
  # run several times!
  fit <- vb(stan_model(stan_file), data = stan_data, output_samples=OUTPUT_SAMPLES, sample_file=samples_file, iter=VBITERS)
  #save(fit.FC, file="fit.vb.FC.Rdata")
}
