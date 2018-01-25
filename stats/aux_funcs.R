#install.packages("ppcor", repos="http://cran.us.r-project.org")

library(R.matlab)
library(pracma)
library(reshape)
library(dplyr)
#library(rstan)
library(ggplot2)
library(data.table)
library(abind)
library(ppcor)
library(yhat)
library(mclust)
library(LSD)
#library(interplot)
library(foreach)
library(doMC)
library(RColorBrewer)
#library(lmerTest)
#library(sem)
library(mvtnorm)
library(parallel)
library(mice)
#library(caret)
#library(glmnet)

N_CORES = 10
registerDoMC(cores=20)

MYPALETTE=brewer.pal(11, "Spectral")[seq(11, 1, -1)]# colorspace::diverge_hsv(15, power=2)
#CLUSPAL = c("#00FFFFFF", "#FFFF00FF", "#2A2A2A", "#FF0000FF", "#00FF00FF", "#0000FFFF")

PAL2 = c('yellow', 'red')
CLUSPAL = c("#00FFFF", "#FFFF00", "#707070", "#FF0000", "#00FF00", "#0000FF")

CLUSPAL9 = c("red", "green", "blue", 
             "yellow", "gray50", "violet",
             "black", "cyan", "orange")

CLUSPAL18 = rep(c("red", "green", "blue", 
                  "yellow", "gray50", "violet",
                  "black", "cyan", "orange"), 2)


#CLUSPAL18= c(hsv(h=seq(0, 1, 1/8)*.9, v=.7),  hsv(h=seq(0, 1, 1/8)*.9, v=1))

#CLUSPAL = c("#49ae8a", "#ab62c0", "#7ca343", "#6587cd", "#c57c3c", "#ca5670")
# right resolution
BWRES=1200
CRES=1200
POINTSIZE=35
CEX_AXIS=2
CEX_LAB=2
CEX_NAMES=1.5
CEX_LEG=2
CEX_MAIN=2

LWD=3
FINELINE=.6
GGTEXTSIZE1=20
GGTEXTSIZE2=60#50
GGTEXTSIZE3=25#20
GGTEXTSIZE4=75#60

blank_theme <- theme(
  axis.title.x = element_text(size = GGTEXTSIZE2),
  axis.title.y = element_text(size = GGTEXTSIZE2),
  panel.grid = element_blank(),
  axis.text.x = element_text(size = GGTEXTSIZE2), 
  axis.text.y = element_text(size = GGTEXTSIZE2),
  legend.title = element_blank(), 
  legend.text = element_text(size = GGTEXTSIZE2)
)


# params for connection removal 
NREPS = 20
NVALUES = 21

#NREPS = 5
#NVALUES = 5


plot_adj = function(adj, modules_file, lim=NULL, pal=MYPALETTE, nolegend=FALSE, use.valid=T){
  
  # plot adjacency matrix
  load(modules_file)
  
  if (use.valid) {
    partition = partition[valid.indices]
    adj = adj[valid.indices, valid.indices]
  
  } 

  ord = order(partition)
  partition.ord = partition[ord]
  
  adj = adj[ord, ord]
  
  adj = adj[partition.ord>0, partition.ord>0]
  partition.ord = partition.ord[partition.ord>0]
  rownames(adj) = colnames(adj) = seq(nrow(adj))
  
  #module_names = sort(unique(names(partition.ord)))
  breaks = which(diff(partition.ord)==1) + .5
  maxval = length(partition.ord) + 0.5
  midpoints = (c(0.5, breaks) + c(breaks, maxval )) / 2
  
  molten_adj = melt(adj)
  plot_data = rename(molten_adj, from = X1, to = X2) 
  
  if (nolegend) {
    
    theme_new = theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),   
      legend.text = element_text(size=GGTEXTSIZE2),
      legend.title = element_blank(), 
      aspect.ratio = 1,
      legend.position = "none"
    )
    
  } else {
    
    theme_new = theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),   
      legend.text = element_text(size=GGTEXTSIZE2 - 5),
      legend.title = element_blank(), 
      aspect.ratio = 1)
  }
  
  
  if (is.null(lim)) {
    #lim = max(abs(adj), na.rm=TRUE)
    lim = 1
    lim1 = quantile(abs(adj), .95, na.rm=TRUE)
    lim0 = -lim1 
  }
  
  if (lim == -1) {
    myplot <- ggplot(plot_data, aes(x = from, y = to, fill = factor(value))) +
      geom_raster() +
      theme_bw() +
      theme_new + 
      scale_fill_manual(values=pal) #+
    #scale_y_reverse()
  }

  if (lim == -2) {
    lim = 1
    lim0 = quantile(abs(adj), .05, na.rm=TRUE)
    lim1 = quantile(abs(adj), .95, na.rm=TRUE)
  }  
   
  if (0 < lim ) {
    myplot <- ggplot(plot_data, aes(x = from, y = to, fill = value)) +
      geom_raster() +
      theme_bw() +
      theme_new + 
      scale_fill_gradientn(colours = pal, 
                           limits=c(lim0, lim1), na.value = "white") #+
    #  scale_y_reverse()
  }
  
  
  # separate in modules
  
  for (i in seq(length(breaks))){
    myplot <- myplot + 
      geom_segment(x = 0.5, xend = maxval, y = breaks[i], yend = breaks[i],  size = FINELINE) + 
      geom_segment(y = 0.5, yend = maxval, x = breaks[i], xend = breaks[i],  size = FINELINE)
    
  }
  
  myplot <- myplot + 
    geom_segment(x = 0.5, xend = maxval, y = 0.5, yend = 0.5,  size = FINELINE) + 
    geom_segment(x = 0.5, xend = maxval, y = maxval, yend = maxval,  size = FINELINE) + 
    geom_segment(y = 0.5, yend = maxval, x = 0.5, xend = 0.5,  size = FINELINE) + 
    geom_segment(y = 0.5, yend = maxval, x = maxval, xend = maxval,  size = FINELINE) 
  
  if (nolegend) MYTEXTSIZE = GGTEXTSIZE3
  else MYTEXTSIZE = GGTEXTSIZE3 - 3 
  
    
  for (i in seq(length(midpoints))){
    myplot <- myplot + 
      annotate("text", x = -.5, y = midpoints[i], label = module_names[i], size = MYTEXTSIZE, hjust=1) +
      annotate("text", y = -.5, x = midpoints[i], label = module_names[i], size = MYTEXTSIZE, hjust=1, angle=90) 
  }
  
  myplot = myplot + ylim(-10, maxval) + xlim(-10, maxval)
  print(myplot)
}

save_fig = function(figname=NULL, width=6.5, height=6.5, res=600, jpg=F){
  #size in inches
  if (dev.cur()!=1) dev.off()
  if (is.null(figname)) {
    figname = paste0('plot', fig_count)
    
    fig_count <<- fig_count+1
  } 
  
  print(paste("Generating figure: ", figname))
  if (jpg){
    figname = paste0(FIGS_DIR, figname, '.jpg')
    jpeg(figname, width = floor(res/2.54*width), height = floor(res/2.54*height), pointsize=POINTSIZE)
  } else {
    figname = paste0(FIGS_DIR, figname, '.png')
    png(figname, width = floor(res/2.54*width), height = floor(res/2.54*height), pointsize=POINTSIZE)
  }
  
  
}

extract_means = function(X){
  slope_muc = colMeans(X[, grep("slope_muc\\.", colnames(X)), with=FALSE])
  intercept_muc = colMeans(X[, grep("intercept_muc\\.", colnames(X)), with=FALSE])
  log_slope_sigmac = colMeans(X[, grep("log_slope_sigmac\\.", colnames(X)), with=FALSE])
  log_intercept_sigmac = colMeans(X[, grep("log_intercept_sigmac\\.", colnames(X)), with=FALSE])
  return(list(slope_muc = slope_muc, intercept_muc = intercept_muc, 
              log_slope_sigmac = log_slope_sigmac, log_intercept_sigmac = log_intercept_sigmac ))
}

iszero = function(X){
  
  quants = apply(X[, , with=FALSE], 2, quantile,  c(.025, .975))
  zeros = (quants[1, ] < 0) & (quants[2, ] > 0)
  
  slope_muc = zeros[grep("slope_muc\\.", names(zeros))]
  intercept_muc = zeros[grep("intercept_muc\\.", names(zeros))]
  log_slope_sigmac = zeros[grep("log_slope_sigmac\\.", names(zeros))]
  log_intercept_sigmac = zeros[grep("log_intercept_sigmac\\.", names(zeros))]
  
  return(list(slope_muc = slope_muc, intercept_muc = intercept_muc, 
              log_slope_sigmac = log_slope_sigmac, log_intercept_sigmac = log_intercept_sigmac ))
}

get_adj = function(values, modules_file, sort=FALSE, labels_file=NULL){
  
  load(modules_file)
  
  n = n.rois
  adj = matrix(NA, n, n)
  adj[valid.indices, valid.indices] = squareform2(values)
  
  if (!is.null(labels_file)){
    labels = read.csv(labels_file, header=TRUE, sep='\t')
    rownames(adj) = colnames(adj) = labels$Label
    
  }
  
  # order according to modules 
  if (sort) 
  {
    ord = order(partition)
    adj = adj[ord, ord]
  }
  
  return(adj)
}

get_adj_sparse = function(values, modules_file, sort=FALSE, labels_file=NULL, valid=NULL){
  
  load(modules_file)
  
  n = n.rois
  
  adj = matrix(NA, n, n)
  adj.s = squareform2(adj)
  adj.s[valid] = values
  adj = squareform2(adj.s)
  
  if (!is.null(labels_file)){
    labels = read.csv(labels_file, header=TRUE, sep='\t')
    rownames(adj) = colnames(adj) = labels$Label
    
  }
  
  # order according to modules 
  if (sort) 
  {
    ord = order(partition)
    adj = adj[ord, ord]
  }
  
  return(adj)
}


get_adj_modules = function(values, modules_file){
  
  load(modules_file)
  
  n = length(ICNnames)
  n_connections = n*(n-1)/2
  
  adj = squareform(values[seq(n_connections)])
  diag(adj) = values[(n_connections + 1): length(values)]
  rownames(adj) = colnames(adj) = ICNnames
  
  return(adj)
}

order_vals = function(values, modules_file, ismatrix=FALSE){
  
  load(modules_file)
  ord = order(partition)
  if (ismatrix) {
    values = values[ord, ord]
  }
  
  else {
    if (is.null(dim(values))) values = values[ord]
    else values = values[, ord]
  }
  return(values)
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, 100, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

diagnostics = function(vals, meas, params){
  pairs(vals, diag.panel=panel.hist, main = meas, labels=params, pch=20, cex=1)
  
}

evaluate = function(FIGS_DIR, PARAMS_FILE, title){
  print(paste("Evaluating model in file: ", PARAMS_FILE))
  dir.create(FIGS_DIR)
  fit = fread(PARAMS_FILE, header=TRUE)
  save_fig(figname="scatterplot", jpg=T)
  
  outlier_rows = find_outliers(fit[, params, with=FALSE])
  #X = fit[!outlier_rows, ]
  X = fit
  means = extract_means(X)
  zeros = iszero(X)
  
  diagnostics(X[, params, with=FALSE], title, params)
  
  means.params = colMeans(X[, params, with=FALSE])
  quants.params = apply(X[, params, with=FALSE], 2, quantile,  c(.025, .975))
  print(means.params)
  print(quants.params)
  return(list(means = means, zeros = zeros))
}

check_chains = function(FIGS_DIR, PARAMS_FILE_LIST, title){
  print(paste("Checking chains for model in files: ", PARAMS_FILE_LIST))
  
  fit = read_stan_csv(PARAMS_FILE_LIST)
  save_fig(figname="traceplot", jpg=T)
  print(rstan::traceplot(fit, pars = params))
  Rhat = summary(fit)$summary[, "Rhat"] 
  n_eff = summary(fit)$summary[, "n_eff"] 
  
  print("min n_eff")
  min_neff = min(n_eff, na.rm=TRUE)
  print(n_eff[n_eff = min_neff])

  print("Rhat above 1.1 : ")
  print(Rhat[Rhat > 1.1])
  
}

evaluate_averages = function(PARAMS_FILE, MODULES_FILE, ICN, bymodules = FALSE){
  print(paste("Evaluating averages for model in file: ", PARAMS_FILE))
  
  fit = fread(PARAMS_FILE, header=TRUE)
  
  outlier_rows = find_outliers(fit[, params, with=FALSE])
  
  #X = fit[!outlier_rows, ]
  X = fit
  
  slope_muc = X[, grep("slope_muc\\.", names(X)), with=FALSE]
  intercept_muc = X[, grep("intercept_muc\\.", names(X)), with=FALSE]
  log_slope_sigmac = X[, grep("log_slope_sigmac\\.", names(X)), with=FALSE]
  log_intercept_sigmac = X[, grep("log_intercept_sigmac\\.", names(X)), with=FALSE]
  
  if (!bymodules) {
    slope_muc.average = average_params_variability(slope_muc, MODULES_FILE, ICN=ICN)
    intercept_muc.average = average_params_variability(intercept_muc, MODULES_FILE, ICN=ICN)
    log_slope_sigmac.average = average_params_variability(log_slope_sigmac, MODULES_FILE, ICN=ICN)
    log_intercept_sigmac.average = average_params_variability(log_intercept_sigmac, MODULES_FILE, ICN=ICN)
  } else {
    slope_muc.average = average_params_variability_modules(slope_muc, MODULES_FILE, ICN=ICN)
    intercept_muc.average = average_params_variability_modules(intercept_muc, MODULES_FILE, ICN=ICN)
    log_slope_sigmac.average = average_params_variability_modules(log_slope_sigmac, MODULES_FILE, ICN=ICN)
    log_intercept_sigmac.average = average_params_variability_modules(log_intercept_sigmac, MODULES_FILE, ICN=ICN)
    
  }
  
  return(list(
    slope_muc = slope_muc.average, 
    intercept_muc = intercept_muc.average,
    log_slope_sigmac = log_slope_sigmac.average, 
    log_intercept_sigmac = log_intercept_sigmac.average)
  )
}

predict_value = function(age.ref, PARAMS_FILE, musigma="mu"){
  print(paste("Predicting value for model in file: ", PARAMS_FILE))
  
  X = fread(PARAMS_FILE, header=TRUE)
  
  slope_muc = X[, grep("slope_muc\\.", names(X)), with=FALSE]
  intercept_muc = X[, grep("intercept_muc\\.", names(X)), with=FALSE]
  log_slope_sigmac = X[, grep("log_slope_sigmac\\.", names(X)), with=FALSE]
  log_intercept_sigmac = X[, grep("log_intercept_sigmac\\.", names(X)), with=FALSE]
  
  
  if (musigma == "mu") {
    m = intercept_muc + slope_muc * age.ref
  } else {
    m = log_intercept_sigmac + log_slope_sigmac * age.ref
  }
  
  value = list()
  value$mean = colMeans(m)
  
  quants = apply(m, 2, quantile,  c(.025, .975))
  value$zeros = (quants[1, ] < 0) & (quants[2, ] > 0)
  
  return(value)
}

get_samples = function(param, PARAMS_FILE){
  print(paste("Getting Sample from file: ", PARAMS_FILE))
  
  X = fread(PARAMS_FILE, header=TRUE)
  
  samples = as.matrix(X[, grep(param, names(X)), with=FALSE])
  
  return(samples)
}

correlate_samples = function(X, Y, Z=NULL){
  c = list(samples=rep(0, nrow(X)))
  
  if (is.null(Z)){
    for (i in 1:nrow(X)){
      s = sample(ncol(X), replace=T)
      c$samples[i] = cor.test(X[i, s], Y[i, s])$estimate
    }
  } else {
    for (i in 1:nrow(X)){
      s = sample(ncol(X), replace=T)
      c$samples[i] = pcor.test(X[i, s], Y[i, s], Z[i, s])$estimate
    }
  }
  
  
  c$CI = quantile(c$samples, c(.025, .975))
  c$mean = mean(c$samples)
  return(c)
}

find_outliers = function(data){
  #print(colMeans(data))
  
  outliers = apply(data, 2, function(x){x %in% boxplot.stats(x)$out})
  if (0 < length(outliers)) { 
    outlier_rows = apply(outliers, 1, any)
    print(paste("Found ", sum(outlier_rows), "outliers"))
  }
  else outlier_rows=NULL
  
  #print(colMeans(data[!outlier_rows, ]))
  
  return(outlier_rows)
}

compute_degree = function(X, correct=T){
  
  c = ifelse(correct, 1, 0)
  deg = function(y){ 
    valid = !is.na(y)  
    if (sum(valid)==0){
      deg = NA
    } else {
      deg = sum(y[valid])/(sum(valid)-c)
    }
    
  }
  
  degree = apply(X, 2, deg)
  
}

compute_degree_distance = function(adj, dist.matrix, d, short=T){
  if (short) {
    adj[dist.matrix > d] = NA
    correct = TRUE
  }
  else 
  {
    adj[dist.matrix < d] = NA
    correct = FALSE
  }
  
  
  degree = compute_degree(adj, correct)
  
  
}

compute_degree_kernel = function(adj, dist.matrix, center, FWHM){
  
  K = (2*sqrt(2*log(2)))
  sigma=(FWHM/K)
  
  adj = adj * dnorm(dist.matrix, center, sigma) # exp(-K*(dist.matrix-center)**2/FWHM2)
  
  degree = compute_degree(adj, TRUE)
  
  
}

compute_degree_samples = function(X, modules_file, labels_file, modular=F){
  
  labels = read.csv(labels_file, header=TRUE, sep='\t')
  
  degree = matrix(0, nrow(X), nrow(labels))
  colnames(degree) = labels$Label
  
  for (i in 1:nrow(X)){
    adj = get_adj(X[i, ], modules_file, sort=TRUE)
    
    if (modular) {
      adj[outer(partition, partition, "!=")] = NA
    }
    degree[i, ] = compute_degree(adj)
  }
  
  return(degree)
}

plot_correl = function(y, x, ylab, xlab, valid=abs(y)>=0, cex=1, asp = -1, right = F ){
  c = cor.test(x[valid], y[valid])
  r = format(round(c$estimate, 3), nsmall=3)
  plot(x[valid], y[valid], xlab=xlab, ylab=ylab, pch=20, cex=cex, cex.axis=CEX_AXIS, cex.lab=CEX_LAB, asp = asp, main="")
  
  if (c$p.value<=0.05) {
    model = lm(y[valid] ~ x[valid])
    abline(model, lwd=LWD)
  }
  
  if (asp==1) {
    abline(0, 1, lwd=LWD, col="red")
  }
  
  
  options(scipen=3)
  if (c$p.value>=0.001) {
    p = format(round(c$p.value, 3), nsmall=3)
  } else if (c$p.value>0){
    p = format(c$p.value, scientific=TRUE, digits=3)
  }
  else p = 0
  
  side = ifelse(right, "topright", "topleft")   
  
  legend(side, legend = c( paste0('r=', r), paste0('p=', p)), bty="n", cex=CEX_LEG )
  
}

plot_correl_repeated = function(y, x, ylab, xlab, cex=1, asp = -1 ){
  
  plot_data = melt(x)
  molten.y = melt(y)
  plot_data$value.y = molten.y$value
  
  plot_data = subset(plot_data, !is.na(value.y) & (value.y!=0)) 
  model = lmer( value.y ~ value + (1 | X2), data = plot_data)
  p.value = coef(summary(model))["value", "Pr(>|t|)"] #Satterwaithe approx
  intercept = coef(summary(model))["(Intercept)", "Estimate"]
  slope = coef(summary(model))["value", "Estimate"]
  r = format(round(slope, 3), nsmall=3)
  
  heatscatter(plot_data$value, plot_data$value.y, xlab=xlab, ylab=ylab, pch=20, cex=cex, cex.axis=CEX_AXIS, cex.lab=CEX_LAB, asp = asp, main ="")
  
  if (p.value<=0.5) {
    abline(intercept, slope, lwd=LWD)
  }
  
  if (asp==1) {
    abline(0, 1, lwd=LWD, col="red")
  }
  
  options(scipen=3)
  if (p.value>=0.001) {
    p = format(round(p.value, 3), nsmall=3)
  } else if (p.value>0){
    p = format(p.value, scientific=TRUE, digits=3)
  }
  else p = 0
  
  #legend("topleft", legend = c( paste0('r=', r), paste0('p=', p)), bty="n", cex=CEX_LEG )
  legend("topleft", legend = paste0('p=', p), bty="n", cex=CEX_LEG )
  
}

adj_cor <- function(X){
  n = dim(X)[1] 
  r = matrix(NA, n, n)
  for (i in seq(n)){
    for (j in seq(n)){
      
      r[i, j] = cor.test(X[i, ], X[j, ])$estimate
      
    }
  }
  return(r)  
}

plot_similarities = function(FC.ordered, sorted.age, line=NULL){
  D = adj_cor(FC.ordered)
  
  molten_D = melt(D)
  plot_data = rename(molten_D, from = X1, to = X2)
  plot_data$age.1 = sorted.age[plot_data$from]
  plot_data$age.2 = sorted.age[plot_data$to]
  
  myplot <- ggplot(plot_data, aes(x = from, y = to, fill = value)) +
    geom_raster() +
    theme_bw() +
    #  scale_y_reverse() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),       
      legend.text = element_text(size=GGTEXTSIZE2 + 15),
      legend.title = element_blank(), 
      aspect.ratio = 1) +
    scale_fill_gradientn(colours = colorspace::diverge_hsv(15), 
                         limits=c(0, 1), na.value = "white")
  
  if (!is.null(line)){
    
    myplot <- myplot + 
      geom_hline(yintercept = line, size = FINELINE) + 
      geom_vline(xintercept = line, size = FINELINE) 
    
    myplot <- myplot + 
      annotate("text", x = line/2, y = -2, label = "Younger", size = GGTEXTSIZE1 + 10) +
      annotate("text", x = (line + dim(D)[1])/2, y = -2, label = "Older", size = GGTEXTSIZE1 + 10) + 
      annotate("text", y = line/2, x = -2, label = "Younger", size = GGTEXTSIZE1 +10, angle = 90) +
      annotate("text", y = (line + dim(D)[1])/2, x = -2, label = "Older", size = GGTEXTSIZE1 + 10, angle = 90) 
  }
  
  print(myplot)
}

plot_mds = function(FC, group){
  D = adj_cor(FC)
  D.old = squareform2(D[group=='Old', group=='Old'])
  D.young = squareform2(D[group=='Young', group=='Young'])
  D.oldyoung = D[group=='Old', group=='Young']

  print(t.test(D.old, D.young))
  print(t.test(D.old, D.oldyoung))

  mds = cmdscale(as.dist(1 - D), eig=T, k=2)
  x1 = mds$points[, 1]
  x2 = mds$points[, 2]
  par(mar=c(8,8,5,5), mgp = c(4, 2, 0))
  
#  plot(x1, x2, pch = ifelse(group=="Young", 21, 24), xaxt='n', yaxt='n', ann=FALSE, bty='n', cex=4)
#  plot(x1, x2, pch = ifelse(group=="Young", 21, 24), xlab='Dimension 1', ylab='Dimension 2', cex=4, bg= 'grey50', lwd= 3)
  plot(x1, x2, pch = ifelse(group=="Young", 2, 1), xlab='Dimension 1', ylab='Dimension 2', cex = 4, cex.lab = 3, cex.axis = 3, lwd = 11)
  
#  legend("bottomleft", legend=c("Younger", "Older"), pch =c(21, 24), pt.bg= 'grey50', cex = CEX_LAB)
  legend("bottomleft", legend=c("Younger", "Older"), pch = c(2, 1), cex = CEX_LAB + 2, lwd = 11, lty = 0)
  
}

plot_connectome_similarities = function(input_file, demo, plot_mds=F){
  Conn_info = readMat(input_file)
  labels = unlist(Conn_info$labels)
  
  subjects_FC = unlist(Conn_info$subjects)
  X = data.frame(Conn_info$merged.matrices.mat)
  
  vars = colnames(X)
  n_conn = length(vars)
  X$Subject = subjects_FC
  
  X = merge(X, demo, by="Subject")
  FC = X[, vars]
  
  age = X$age
  group = X$class
  
  FC.ordered = as.matrix(arrange(data.frame(FC), age))
  sorted.age = sort(age)
  lim = which(sorted.age>50)[1]-0.5
  if (plot_mds) plot_mds(as.matrix(FC), group)
  else plot_similarities(FC.ordered, sorted.age, lim)
  
}

plot_connection_removal = function(FC.ordered, values, conn.order, group1, group2, main, group.ordered, ylabel, lim=length(conn.order)){
  print("Doing connection removal")
  
  var.ordered = var.ordered.sampled = c()
  
  d.ordered.sampled = c()
  removed = round(seq(0, lim, length.out = NVALUES))[1:(NVALUES-1)]
  #removed = seq(0, dim(FC.ordered)[2], length.out = NVALUES)[1:(NVALUES-1)]
  
  sampled_connections = c()
  for (k in seq(NREPS)){ 
    sampled_connections = rbind(sampled_connections, 
                                sample(conn.order, replace=F))
  }
  
  
  for (rem_cons in removed) {
    print(rem_cons)
    sel = conn.order[(rem_cons + 1):length(conn.order)]
    
    d.ordered = adj_cor(FC.ordered[, sel])
    d.ordered = d.ordered - diag(diag(d.ordered))
    d.ordered = d.ordered[group.ordered==group1, group.ordered==group2]
    
    if (group1==group2) myvals = squareform(d.ordered)
    else myvals = d.ordered
    
    var.ordered = c(var.ordered, mean(myvals))
    #  image(d.ordered, main=paste("Ordered", p), zlim = conlims, asp = 1)  
    
    myvar = c()
    for (k in 1:NREPS){  
      mymat = adj_cor(FC.ordered[, sampled_connections[k, 1:(length(conn.order)-rem_cons)]])
      d.ordered.sampled = abind(d.ordered.sampled, mymat, along = 3)
      mymat = mymat - diag(diag(mymat))
      mymat = mymat[group.ordered==group1, group.ordered==group2]
      
      if (group1==group2) myvals = squareform(mymat)
      else myvals = mymat
      
      myvar = c(myvar,  mean(myvals))
    }
    var.ordered.sampled = cbind(var.ordered.sampled, myvar)  
    #image(apply(d.ordered.sampled, c(1, 2), mean), main=paste("Shuffled", p), zlim = conlims, asp = 1)
  }
  
  plot(removed, var.ordered, type="b", col="red", lwd=3, pch=20, cex.axis=CEX_AXIS, cex.lab=CEX_LAB, cex.main = CEX_MAIN, ylim=c(0.3, 0.5),
       ylab="Average similarity between subjects", xlab = "Number of removed connections", main = main)
  
  for (i in 1:NREPS) lines(removed, var.ordered.sampled[i,], type="b", col="blue", pch=20, cex = 0.5)
  lines(removed, var.ordered, type="b", col="red", lwd=3, pch=20, cex = 0.5)
  par(new=T)
  plot(values[conn.order[1:removed[length(removed)]]], pch=20, cex = 0.5, col="green", xaxt="n",yaxt="n",xlab="",ylab="")
  axis(4, cex.axis=CEX_AXIS) 
  mtext(ylabel,side=4,line=4, col="green",  cex=CEX_LAB)
  
  #legend("topright", legend=c("Ordered","Shuffled"), fill=c("red", "blue"), cex=CEX_LEG)
}

plot_connectome_removal = function(input_file, demo, conn.order.mu, conn.order.sigma, values_mu, values_sigma, lim.muc = length(conn.order.mu), lim.sigmac = length(conn.order.sigma), figname="FigureRemoval") {
  Conn_info = readMat(input_file)
  labels = unlist(Conn_info$labels)
  
  subjects_FC = unlist(Conn_info$subjects)
  X = data.frame(Conn_info$merged.matrices.mat)
  
  vars = colnames(X)
  n_conn = length(vars)
  X$Subject = subjects_FC
  
  X = merge(X, demo, by="Subject")
  FC = X[, vars]
  
  age = X$age
  group = X$class
  FC.ordered = as.matrix(arrange(data.frame(FC), age))
  sorted.age = sort(age)
  group.ordered = as.matrix(arrange(data.frame(group), age))
#  par(mfrow=c(2, 2), mar=c(8,8,5,5), mgp=c(4,1,0))
  
  save_fig(figname=paste0(figname,"a"), res=BWRES)
  par(mar=c(8,8,5,5), mgp=c(4,1,0))
  plot_connection_removal(FC.ordered, values_mu, conn.order.mu, 'Young', 'Old', expression("Ordering by |" * beta[mu] * "|, Younger-Older"), group.ordered, expression("|"*beta[mu]*"|"), lim=lim.muc)

  save_fig(figname=paste0(figname,"b"), res=BWRES)
  par(mar=c(8,8,5,5), mgp=c(4,1,0))
  plot_connection_removal(FC.ordered, values_mu, conn.order.mu, 'Old', 'Old', expression("Ordering by |" * beta[mu] * "|, Older-Older"), group.ordered, expression("|"*beta[mu]*"|"), lim=lim.muc)

  save_fig(figname=paste0(figname,"c"), res=BWRES)
  par(mar=c(8,8,5,5), mgp=c(4,1,0))
  plot_connection_removal(FC.ordered, values_sigma, conn.order.sigma, 'Young', 'Old', expression("Ordering by " * beta[sigma] * ", Young-Older"), group.ordered, expression(""*beta[sigma]*""), lim = lim.sigmac)

  save_fig(figname=paste0(figname,"d"), res=BWRES)
  par(mar=c(8,8,5,5), mgp=c(4,1,0))
  plot_connection_removal(FC.ordered, values_sigma, conn.order.sigma, 'Old', 'Old', expression("Ordering by " * beta[sigma] * ", Older-Older"), group.ordered, expression(""*beta[sigma]*""), lim= lim.sigmac)
  
}

associate_cognition = function(input_file, demo, MODULES_FILE_NET, MODULES_FILE_MNI, test, permute=FALSE, beta_mu=NULL, beta_sigma=NULL, alpha_mu=NULL){
  
  Conn_info = readMat(input_file)
  labels = unlist(Conn_info$labels)
  
  subjects_FC = unlist(Conn_info$subjects)
  FC.all = data.frame(Conn_info$merged.matrices.mat)
  
  vars = colnames(FC.all)
  n_conn = length(vars)
  FC.all$Subject = subjects_FC
  
  data = merge(FC.all, demo, by="Subject")
  X = data
  
  FC.all = data[, vars]
  FC = X[,vars]
  
  age = X$age
  group = X$class
  cogscore = X[, test]
  
  if(permute) {
    set.seed(42)
    cogscore = sample(cogscore)
    
  }
  
  # find WM network 
  load(MODULES_FILE_NET)
  WMROIs = cogweights.all[rownames(cogweights.all)=='Cognition.Memory.Working', ] > 0
  
  WMnetwork = outer(WMROIs, WMROIs, FUN="|")
  WMnetwork = WMnetwork[valid.indices, valid.indices]
  WMnetwork = squareform2(WMnetwork)
  
  sel = (WMnetwork > 0)
  
  cors.young = atanh(apply(FC[group=="Young", sel], 2, function(x) pcor.test(x, cogscore[group=="Young"], age[group=="Young"])$estimate))
  cors.old = atanh(apply(FC[group=="Old", sel], 2, function(x) pcor.test(x, cogscore[group=="Old"], age[group=="Old"])$estimate))
  
  cors = c(cors.young, cors.old)
  
  association.score.old = association.score.young = sel
  association.score.young[sel] = cors.young
  association.score.old[sel] = cors.old
  association.score = c(association.score.young, association.score.old)
  q.beta_mu = quantile(beta_mu, c(.2, .8))
  q.beta_sigma = quantile(beta_sigma, c(.2, .8))
  
  ex1 = which(beta_mu > q.beta_mu[2] & beta_sigma > q.beta_sigma[2] & alpha_mu > 0 & sel) 
  ex2 = which(beta_mu < q.beta_mu[1] & beta_sigma > q.beta_sigma[2] & alpha_mu > 0 & sel)
  ex3 = which(beta_mu > q.beta_mu[2] & beta_sigma < q.beta_sigma[1] & alpha_mu > 0 & sel)
  ex4 = which(beta_mu < q.beta_mu[1] & beta_sigma < q.beta_sigma[1] & alpha_mu > 0 & sel)
  
  indices = c(ex1, ex2, ex3, ex4)
  clus = c(rep(1, length(ex1)), rep(2, length(ex2)), rep(3, length(ex3)), rep(4, length(ex4)))
  colors = c("red", "blue", "green", "orange")
  young.means = aggregate(association.score.young[indices], by = list(clus), FUN=mean)
  old.means = aggregate(association.score.old[indices], by = list(clus), FUN=mean)
  
  #print(young.means)
  
  save_fig()
  par(mfrow=c(1,1))
  plot(jitter(clus, factor=.7)+.25, association.score.young[indices], pch=20, col = colors[clus], xlab="", ylab="Association with WM score", xlim=c(0.5, 4.5), ylim = c(-.7, .7) , cex=0.8) 
  points(jitter(clus, factor= .7)-.25, association.score.old[indices], pch=17, col = colors[clus] , cex=0.8) 
  points(seq(4) - .25, young.means$x, col = "black", cex = 4, pch=3)
  points(seq(4) + .25, old.means$x, col = "black", cex = 4, pch=3)
  
  print(sum(sel))
  association.score[!sel] = NA
  
  #save_fig(res=CRES)
  #plot_adj(get_adj(WMnetwork, MODULES_FILE_MNI), MODULES_FILE_MNI, lim = 1)
  
  #save_fig()
  #association.score.adj = get_adj(association.score, MODULES_FILE_MNI)
  #plot_adj(association.score.adj, MODULES_FILE_MNI) 
  
  
  return(association.score)
  
}

associate_cognition_component = function(input_file, demo, MODULES_FILE_NET, MODULES_FILE_MNI, test, permute=FALSE, beta_mu=NULL, beta_sigma=NULL, comp, thr){
  
  Conn_info = readMat(input_file)
  labels = unlist(Conn_info$labels)
  
  subjects_FC = unlist(Conn_info$subjects)
  FC.all = data.frame(Conn_info$merged.matrices.mat)
  
  vars = colnames(FC.all)
  n_conn = length(vars)
  FC.all$Subject = subjects_FC
  
  data = merge(FC.all, demo, by="Subject")
  X = data
  
  FC.all = data[, vars]
  FC = X[,vars]
  
  age = X$age
  group = X$class
  cogscore = X[, test]
  
  if(permute) {
    set.seed(42)
    cogscore = sample(cogscore)
    
  }
  
  #comp[is.na(comp)]=0
  comp = comp[valid.indices, valid.indices]
  
  sel = (squareform2(comp) > thr)
  
  print(sum(sel, na.rm=T))
  
  cors = c(atanh(apply(FC[group=="Young", sel], 2, function(x) pcor.test(x, cogscore[group=="Young"], age[group=="Young"])$estimate)),
           atanh(apply(FC[group=="Old", sel], 2, function(x) pcor.test(x, cogscore[group=="Old"], age[group=="Old"])$estimate))
  )
  sel = c(sel, sel)
  association.score = sel
  association.score[sel] = cors
  
  association.score[!sel] = NA
  
  return(association.score)
  
}

associate_cognition_modules_2groups = function(input_file, demo, test, permute=FALSE, beta_mu=NULL, beta_sigma=NULL, alpha_mu=NULL, modules_file){
  
  Conn_info = readMat(input_file)
  labels = unlist(Conn_info$labels)
  
  subjects_FC = unlist(Conn_info$subjects)
  FC.all = data.frame(Conn_info$merged.matrices.mat)
  
  vars = colnames(FC.all)
  n_conn = length(vars)
  FC.all$Subject = subjects_FC
  
  data = merge(FC.all, demo, by="Subject")
  
  X = data
  
  # balance the groups
  min_n = min(table(X$class))
  which.old = which(X$class=='Old')
  which.young = which(X$class=='Young')
  which.balanced = sort(c(sample(which.old, min_n), sample(which.young, min_n)))
#  X = X[which.balanced, ]

  FC.all = data[, vars]
  FC = X[,vars]
  
  age = X$age
  group = X$class
  cogscore = X[, test]

  if(permute) {
    set.seed(42)
    cogscore = sample(cogscore)
    
  }
  
  load(modules_file)
  WMnets = (cogICNs[rownames(cogICNs)=='Cognition.Memory.Working', ] > 0)
  
  WMnetwork = drop(outer(WMnets, WMnets, FUN="|"))

  sel = c(squareform2(WMnetwork), diag(WMnetwork))

  cors.young = atanh(apply(FC[group=="Young", ], 2, function(x) pcor.test(x, cogscore[group=="Young"], age[group=="Young"])$estimate))
  cors.old = atanh(apply(FC[group=="Old", ], 2, function(x) pcor.test(x, cogscore[group=="Old"], age[group=="Old"])$estimate))
  
  cors.young[!sel] = NA
  cors.old[!sel] = NA
  
  cors = c(cors.young, cors.old)
  association.score = c(cors.young, cors.old)
  
  cors.young.adj = get_adj_modules(cors.young, modules_file)
  cors.old.adj = get_adj_modules(cors.old, modules_file)
  
  save_fig(res=CRES)
  plot_matrix(cors.young.adj)
  save_fig(res=CRES)
  plot_matrix(cors.old.adj)
  
  q.beta_mu = quantile(beta_mu, c(.2, .8))
  q.beta_sigma = quantile(beta_sigma, c(.2, .8))
  
  ex1 = which(beta_mu > q.beta_mu[2] & beta_sigma > q.beta_sigma[2] & sel) # & alpha_mu > 0 
  ex2 = which(beta_mu < q.beta_mu[1] & beta_sigma > q.beta_sigma[2] & sel)
  ex3 = which(beta_mu > q.beta_mu[2] & beta_sigma < q.beta_sigma[1] & sel)
  ex4 = which(beta_mu < q.beta_mu[1] & beta_sigma < q.beta_sigma[1] & sel)
  
  indices = c(ex1, ex2, ex3, ex4)
  clus = c(rep(1, length(ex1)), rep(2, length(ex2)), rep(3, length(ex3)), rep(4, length(ex4)))
  colors = c("red", "blue", "green", "orange")
  young.means = aggregate(cors.young[indices], by = list(clus), FUN=mean)
  old.means = aggregate(cors.old[indices], by = list(clus), FUN=mean)
  
  print(young.means)
  
  save_fig()
  par(mfrow=c(1,1))
  plot(jitter(clus, factor=.7) - .25, cors.young[indices], pch=20, col = colors[clus], xlab="", ylab="Association with WM score", xlim=c(0.5, 4.5), ylim = c(-.7, .7) , cex=0.8) 
  points(jitter(clus, factor= .7) + .25, cors.old[indices], pch=17, col = colors[clus] , cex=0.8) 
  points(seq(4) - .25, young.means$x, col = "black", cex = 4, pch=3)
  points(seq(4) + .25, old.means$x, col = "black", cex = 4, pch=3)
  
  #print(sum(sel))
  #association.score[!sel] = NA
  
  #save_fig(res=CRES)
  #plot_adj(get_adj(WMnetwork, MODULES_FILE_MNI), MODULES_FILE_MNI, lim = 1)
  
  #save_fig()
  #association.score.adj = get_adj(association.score, MODULES_FILE_MNI)
  #plot_adj(association.score.adj, MODULES_FILE_MNI) 
  
  
  return(association.score)
  
}

normal_pattern_sim = function(input_file, demo, test, permute=FALSE, beta_mu=NULL, beta_sigma=NULL, alpha_mu=NULL, modules_file){
  
  Conn_info = readMat(input_file)
  labels = unlist(Conn_info$labels)
  
  subjects_FC = unlist(Conn_info$subjects)
  FC.all = data.frame(Conn_info$merged.matrices.mat)
  
  vars = colnames(FC.all)
  n_conn = length(vars)
  FC.all$Subject = subjects_FC
  
  data = merge(FC.all, demo, by="Subject")
  
  X = data
  
  # balance the groups
  min_n = min(table(X$class))
  which.old = which(X$class=='Old')
  which.young = which(X$class=='Young')
  which.balanced = sort(c(sample(which.old, min_n), sample(which.young, min_n)))
  #  X = X[which.balanced, ]
  
  FC.all = data[, vars]
  FC = X[,vars]
  
  age = X$age
  group = X$class
  cogscore = X[, test]
  
  if(permute) {
    set.seed(42)
    cogscore = sample(cogscore)
    
  }
  
  load(modules_file)
  WMnets = (cogICNs[rownames(cogICNs)=='Cognition.Memory.Working', ] > 0)
  
  WMnetwork = drop(outer(WMnets, WMnets, FUN="|"))
  
  sel = c(squareform2(WMnetwork), diag(WMnetwork))
  sel2 = beta_mu < median(beta_mu, na.rm=T)
  #sel = sel & sel2
  mean.young = colMeans(FC[group=="Young", sel])
  weights = beta_mu[sel] # - max(beta_mu, na.rm = T)
  similarity = apply(FC[group=="Old", sel], 1, function(x) cor.test(weights*x, weights*mean.young)$estimate )
#  deviation = sqrt(rowSums((beta_mu[sel]*diffs)^2) )
  cortest = pcor.test(similarity, cogscore[group=="Old"], age[group=="Old"])
  print(cortest)
  plot(similarity, cogscore[group=="Old"], pch = 20, cex = 0.5, ylab = "WM score", xlab = "Similarity" ) 

#   return(data.frame(similarity=similarity, Subject=data$Subject[group=="Old"], 
#                     age=age[group=="Old"], combined=cogscore[group=="Old"]))
    return(cortest)
}


associate_cognition_modules = function(input_file, demo, test, permute=FALSE, beta_mu=NULL, beta_sigma=NULL, alpha_mu=NULL, modules_file){
  
  Conn_info = readMat(input_file)
  labels = unlist(Conn_info$labels)
  
  subjects_FC = unlist(Conn_info$subjects)
  FC.all = data.frame(Conn_info$merged.matrices.mat)
  
  vars = colnames(FC.all)
  n_conn = length(vars)
  FC.all$Subject = subjects_FC
  
  data = merge(FC.all, demo, by="Subject")
  
  X = data
  
  # balance the groups
  min_n = min(table(X$class))
  which.old = which(X$class=='Old')
  which.young = which(X$class=='Young')
  which.balanced = sort(c(sample(which.old, min_n), sample(which.young, min_n)))
  #  X = X[which.balanced, ]
  
  FC.all = data[, vars]
  FC = X[,vars]
  
  age = X$age
  group = X$class
  cogscore = X[, test]
  
  if(permute) {
    set.seed(42)
    cogscore = sample(cogscore)
    
  }
  
  load(modules_file)
  WMnets = (cogICNs[rownames(cogICNs)=='Cognition.Memory.Working', ] > 0)
  
  WMnetwork = drop(outer(WMnets, WMnets, FUN="|"))
  
  sel = c(squareform2(WMnetwork), diag(WMnetwork))
  
  cors.old = atanh(apply(FC[group=="Old", ], 2, function(x) pcor.test(x, cogscore[group=="Old"], age[group=="Old"])$estimate))
  cors.old[!sel] = NA
  
  association.score = cors.old
  
  cors.old.adj = get_adj_modules(cors.old, modules_file)
  
  save_fig(res=CRES)
  plot_matrix(cors.old.adj)
  
  q.beta_mu = quantile(beta_mu, c(.3, .7))
  q.beta_sigma = quantile(beta_sigma, c(.3, .7))
  
  ex1 = which(beta_mu > q.beta_mu[2] & beta_sigma > q.beta_sigma[2] & sel) # & alpha_mu > 0 
  ex2 = which(beta_mu < q.beta_mu[1] & beta_sigma > q.beta_sigma[2] & sel)
  ex3 = which(beta_mu > q.beta_mu[2] & beta_sigma < q.beta_sigma[1] & sel)
  ex4 = which(beta_mu < q.beta_mu[1] & beta_sigma < q.beta_sigma[1] & sel)
  
  indices = c(ex1, ex2, ex3, ex4)
  
  clus = c(rep(1, length(ex1)), rep(2, length(ex2)), rep(3, length(ex3)), rep(4, length(ex4)))
  colors = c("red", "blue", "green", "orange")
  old.means = aggregate(cors.old[indices], by = list(clus), FUN=mean)

  save_fig()
  par(mfrow=c(1,1))
  plot(jitter(clus, factor=.7), cors.old[indices], pch=20, col = colors[clus], xlab="", ylab="Association with WM score", xlim=c(0.5, 4.5), ylim = c(-.7, .7) , cex=0.8) 
  points(seq(4), old.means$x, col = "black", cex = 4, pch=3)
  
  return(association.score)
  
}


model.cognitive.association = function(association.score, slope_muc, log_slope_sigmac, figname, intercept_muc=NULL, log_intercept_sigmac=NULL, plotme = F){
  
  log_slope_sigmac.scaled = scale(log_slope_sigmac, scale=TRUE)
  slope_muc.scaled = scale(slope_muc, scale=TRUE)
  log_intercept_sigmac.scaled = scale(log_intercept_sigmac, scale=TRUE)
  intercept_muc.scaled = scale(intercept_muc, scale=TRUE)
  

  beta_sigma = as.numeric(log_slope_sigmac.scaled)
  beta_mu = as.numeric(slope_muc.scaled)
  alpha_sigma = as.numeric(log_intercept_sigmac.scaled)
  alpha_mu = as.numeric(intercept_muc.scaled)
  
#   save_fig()
#   par(mfrow=c(3,1))
# #  plot(as.factor(ifelse(group==1, "Younger", "Older")), association.score, ylim=c(-.3, .3))
#   hist(association.score, 100)
#   heatscatter(slope_muc, association.score, cex = .5, main="Older")
#   heatscatter(log_slope_sigmac, association.score, cex = .5, main="Older")
  print(sum(!is.na(association.score)))
  model = lm(association.score ~ beta_mu + beta_sigma + alpha_mu + alpha_sigma )
#  model = lm(association.score ~ beta_mu + beta_sigma + beta_mu*beta_sigma )
  
  print(summary(model))
  coefs = coef(summary(model))
  #model = lm(association.score ~ beta_mu + beta_sigma )
  
  
#   N=100
#   seq_mu.scaled = seq(min(beta_mu, na.rm=TRUE), max(beta_mu, na.rm=TRUE), length.out=N)
#   seq_sigma.scaled = seq(min(beta_sigma, na.rm=TRUE), max(beta_sigma, na.rm=TRUE), length.out=N)
#   seq_mu  = seq_mu.scaled*attr(slope_muc.scaled,"scaled:scale") + attr(slope_muc.scaled,"scaled:center")
#   seq_sigma  = seq_sigma.scaled*attr(log_slope_sigmac.scaled,"scaled:scale") + attr(log_slope_sigmac.scaled,"scaled:center")
#   
#   x = expand.grid(beta_mu=seq_mu.scaled, beta_sigma=seq_sigma.scaled)
#   z = predict(model, x)
#   x = expand.grid(beta_mu=seq_mu, beta_sigma=seq_sigma)
#   d = data.frame(z=z, beta_mu = x$beta_mu, beta_sigma=x$beta_sigma)
if (plotme){
  
  save_fig(figname  = figname, res=CRES)
  
#   blank_theme <- theme_minimal() +
#     theme(
#       axis.title.x = element_text(size = GGTEXTSIZE2),
#       axis.title.y = element_text(size = GGTEXTSIZE2),
#       panel.grid = element_blank(),
#       axis.text.x = element_text(size = GGTEXTSIZE2), 
#       axis.text.y = element_text(size = GGTEXTSIZE2),
#       legend.title = element_blank(), 
#       legend.text = element_text(size = GGTEXTSIZE2)
#     )
#   
#   slopes = data.frame(slope_muc = slope_muc, log_slope_sigmac=log_slope_sigmac) 
#   
#   myplot <- ggplot(d, aes(x=beta_mu, y=beta_sigma, z=z)) + 
#     geom_tile(aes(fill=z)) + scale_fill_gradientn(limits=c(-.3, .3), colours=MYPALETTE) + #scale_fill_gradient2(limits=c(-.3, .3), low="blue", high="red") +
#     geom_point(data = slopes, aes(x=slope_muc, y= log_slope_sigmac, z=NULL), alpha=0.3, size=0.5) +
#     xlab(expression(beta[mu])) + ylab(expression(beta[sigma])) + blank_theme + ylim(.001, .008) +  xlim(-.004, .002) 
# 
#   print(myplot)
  #model.diff = lm(association.score ~ slope_muc.scaled + intercept_muc.scaled + log_intercept_sigmac.scaled )
  #print(summary(model.diff))
  par(mar=c(8,8,5,5), mgp=c(4,1,0))
  plot_correl(association.score[!is.na(association.score)], 
              slope_muc[!is.na(association.score)],
              "z-score", 
              expression(bar(beta)[mu]), 
              cex = 1,
              right = T)
  }            
  #plot_correl(model.diff$residuals, 
  #            log_slope_sigmac[!is.na(association.score)],
  #            "Difference (Younger - Older) in association score with WM (residuals)", 
  #            expression(bar(beta)[sigma]), 
  #            cex = 0.5)
  
  #print(cor.test(log_slope_sigmac.scaled[!is.na(association.score)], model.diff$residuals ))
return(coefs)
}

model.cognitive.association_2groups = function(association.score, slope_muc, log_slope_sigmac, figname, intercept_muc=NULL, log_intercept_sigmac=NULL){
  
  log_slope_sigmac.scaled = scale(log_slope_sigmac, scale=TRUE)
  slope_muc.scaled = scale(slope_muc, scale=TRUE)
  log_intercept_sigmac.scaled = scale(log_intercept_sigmac, scale=TRUE)
  intercept_muc.scaled = scale(intercept_muc, scale=TRUE)
  
  group = c(rep(1, length(slope_muc)), rep(0, length(slope_muc))) # young = 1, old = 0
  
  beta_sigma = as.numeric(rep(log_slope_sigmac.scaled, 2))
  beta_mu = as.numeric(rep(slope_muc.scaled, 2))
  alpha_sigma = as.numeric(rep(log_intercept_sigmac.scaled, 2))
  alpha_mu = as.numeric(rep(intercept_muc.scaled, 2))
  
  save_fig()
  par(mfrow=c(3,2))
  plot(as.factor(ifelse(group==1, "Younger", "Older")), association.score, ylim=c(-.3, .3))
  hist(association.score, 100)
  heatscatter(slope_muc, association.score[group==0], cex = .5, main="Older")
  heatscatter(log_slope_sigmac, association.score[group==0], cex = .5, main="Older")
  heatscatter(slope_muc, association.score[group==1], cex = .5, main="Younger")
  heatscatter(log_slope_sigmac, association.score[group==1], cex = .5, main="Younger")
  
  
  #model = lm(association.score ~ beta_mu * beta_sigma * as.factor(group))
  #model = lm(association.score ~ beta_mu * as.factor(group) + beta_sigma* as.factor(group) + alpha_mu * as.factor(group) + alpha_sigma* as.factor(group))
  #model = lm(association.score ~ 0 + alpha_mu + alpha_sigma)
  
  model = lm(association.score ~ beta_mu * as.factor(group) + beta_sigma* as.factor(group) + alpha_mu + alpha_sigma + beta_mu*beta_sigma* as.factor(group) )
  
  print(summary(model))
  

  N=100
  seq_mu.scaled = seq(min(beta_mu, na.rm=TRUE), max(beta_mu, na.rm=TRUE), length.out=N)
  seq_sigma.scaled = seq(min(beta_sigma, na.rm=TRUE), max(beta_sigma, na.rm=TRUE), length.out=N)
  seq_mu  = seq_mu.scaled*attr(slope_muc.scaled,"scaled:scale") + attr(slope_muc.scaled,"scaled:center")
  seq_sigma  = seq_sigma.scaled*attr(log_slope_sigmac.scaled,"scaled:scale") + attr(log_slope_sigmac.scaled,"scaled:center")
  
  x = expand.grid(beta_mu=seq_mu.scaled, beta_sigma=seq_sigma.scaled, group = c(1, 0), alpha_mu = mean(alpha_mu), alpha_sigma = mean(alpha_sigma))
  z = predict(model, x)
  x = expand.grid(beta_mu=seq_mu, beta_sigma=seq_sigma, group = c(1, 0))
  d = data.frame(z=z, beta_mu = x$beta_mu, beta_sigma=x$beta_sigma, group = 0)#group = factor(x$group))
  
  save_fig(figname  = figname, res=CRES)
  
  blank_theme <- theme_minimal() +
    theme(
      axis.title.x = element_text(size = GGTEXTSIZE2),
      axis.title.y = element_text(size = GGTEXTSIZE2),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = GGTEXTSIZE2), 
      axis.text.y = element_text(size = GGTEXTSIZE2),
      legend.title = element_blank(), 
      legend.text = element_text(size = GGTEXTSIZE2)
    )
  
  slopes = data.frame(slope_muc = slope_muc, log_slope_sigmac=log_slope_sigmac) 
  
  myplot <- ggplot(d, aes(x=beta_mu, y=beta_sigma, z=z)) + 
    geom_tile(aes(fill=z)) + scale_fill_gradientn(limits=c(-.3, .3), colours=MYPALETTE) + #scale_fill_gradient2(limits=c(-.3, .3), low="blue", high="red") +
    geom_point(data = slopes, aes(x=slope_muc, y= log_slope_sigmac, z=NULL), alpha=0.3, size=0.5) +
    xlab(expression(beta[mu])) + ylab(expression(beta[sigma])) + blank_theme + ylim(.001, .008) +  xlim(-.004, .002) 
  #facet_grid(. ~ group) +
  
  #print(myplot)
  diff.association.score = association.score[group==1] - association.score[group==0] 
  model.diff = lm(diff.association.score ~ slope_muc.scaled + intercept_muc.scaled + log_intercept_sigmac.scaled )
  print(summary(model.diff))
  par(mar=c(8,8,5,5), mgp=c(4,1,0))
  
  plot_correl(model.diff$residuals, 
              log_slope_sigmac[!is.na(diff.association.score)],
              "Difference (Younger - Older) in association score with WM (residuals)", 
              expression(bar(beta)[sigma]), 
              cex = 0.5)

    print(cor.test(log_slope_sigmac.scaled[!is.na(diff.association.score)], model.diff$residuals ))
}

cluster_params_manual_inter = function(alpha_mu, beta_mu, beta_sigma, alpha_mu.zeros, beta_mu.zeros, beta_sigma.zeros, pal, figname){
  
  print("Doing clustering")
  clusters = 0*beta_mu + 2
  
  clusters = ifelse((!beta_mu.zeros)&(beta_mu < 0), 1, clusters)
  clusters = ifelse((!beta_mu.zeros)&(beta_mu > 0), 3, clusters)
  
  clusters = ifelse(alpha_mu.zeros, clusters + 3, clusters)
  clusters = ifelse((!alpha_mu.zeros)&(alpha_mu > 0), clusters + 6, clusters)
  
  clusters = ifelse(!beta_sigma.zeros, clusters + 9, clusters)
  #median_sigma = median(beta_sigma)
  #clusters = ifelse(beta_sigma > median_sigma, clusters + 9, clusters)
  
  save_fig(figname=figname, res=CRES, width = 13)
  #par(mar=c(8,8,5,5), mfrow=c(1,2)), mgp=c(4,1,0))
  par(mar=c(5,6,4,1), mfrow=c(1,2))
  
  plot(alpha_mu, beta_mu, col=pal[clusters], pch = 20, xlab=expression(alpha[mu]), ylab=expression(beta[mu]), cex=0.5, cex.axis=CEX_AXIS, cex.lab=CEX_LAB)
  
  t.clusters = table(clusters)
  print(t.clusters)
  barplot(rbind(t.clusters[1:9], t.clusters[10:18]), beside=T, ylab = "Number of connections", col = rbind(pal[1:9], pal[10:18]), axisnames=F, cex.lab= CEX_LAB, cex.axis = CEX_AXIS, cex.names = CEX_NAMES)
  
  return(clusters)
}

cluster_params_2points = function(mu.1, mu.2, beta_sigma, mu.1.zeros, mu.2.zeros, beta_sigma.zeros, pal, figname){
  
  print("Doing clustering")
  clusters = 0*mu.1 + 2
  
  clusters = ifelse((!mu.1.zeros)&(mu.1 < 0), 1, clusters)
  clusters = ifelse((!mu.1.zeros)&(mu.1 > 0), 3, clusters)
  
  clusters = ifelse(mu.2.zeros, clusters + 3, clusters)
  clusters = ifelse((!mu.2.zeros)&(mu.2 > 0), clusters + 6, clusters)
  
  #clusters = ifelse(!beta_sigma.zeros, clusters + 9, clusters)
  #median_sigma = median(beta_sigma)
  #clusters = ifelse(beta_sigma > median_sigma, clusters + 9, clusters)
  
#  save_fig(figname=figname, res=CRES, width= 13)
  #par(mar=c(8,8,5,5), mfrow=c(1,2)), mgp=c(4,1,0))
#  par(mar=c(5,6,4,1), mfrow=c(1,2))

  save_fig(figname=paste0(figname, "a"), res=CRES)
  par(mar=c(8,8,5,5), mgp=c(4,1,0))
  plot(mu.1, mu.2, col=pal[clusters], pch = 20, xlab="FC 20 years", ylab="FC 70 years", cex=0.5, cex.axis=CEX_AXIS, cex.lab=CEX_LAB)
  abline(0, 1, lwd=4)
  
  abline(lm(mu.2 ~ mu.1), lty=2, lwd=2)
  
  t.clusters = table(clusters)
  print(t.clusters)
  save_fig(figname=paste0(figname, "b"), res=CRES)
  par(mar=c(8,8,5,5), mgp=c(4,1,0))
    #barplot(rbind(t.clusters[1:9], t.clusters[10:18]), beside=T, ylab = "Number of connections", col = rbind(pal[1:9], pal[10:18]), axisnames=F)
  barplot(t.clusters, ylab = "Number of connections", col = pal, axisnames=F, cex.axis=CEX_AXIS, cex.lab=CEX_LAB)
  
  return(clusters)
}

cluster_params_manual = function(beta_mu, beta_sigma, beta_mu.zeros, beta_sigma.zeros, pal, figname){
  
  print("Doing clustering")
  median_sigma = median(beta_sigma)
  clusters = 0*beta_mu + 3
  
  clusters = ifelse((!beta_mu.zeros)&(beta_mu < 0), 1, clusters)
  clusters = ifelse((!beta_mu.zeros)&(beta_mu > 0), 5, clusters)
  clusters = ifelse(!beta_sigma.zeros, clusters + 1, clusters)
  
  save_fig(figname=figname, res=CRES)
  par(mar=c(8,8,5,5), mgp=c(4,1,0))
  plot(beta_mu, beta_sigma, col=pal[clusters], pch = 20, xlab=expression(beta[mu]), ylab=expression(beta[sigma]), cex=0.5, cex.axis=CEX_AXIS, cex.lab=CEX_LAB)
  
  return(clusters)
}

cluster_params = function(slope_muc, log_slope_sigmac){
  
  print("Doing clustering")
  
  X = data.frame(beta_mu=slope_muc, beta_sigma=log_slope_sigmac)
  print(dim(X))
  #BIC = mclustBIC(X)
  save_fig()
  par(mfrow=c(1,3))
  
  molten_data = melt(X)
  heatscatter(slope_muc, log_slope_sigmac, xlab="beta_mu", ylab = "beta_sigma", cex.main= 0.1)
  
  #model = Mclust(X, x = BIC)
  #plot(model, what = "classification")
  #plot(BIC)
  #print(summary(BIC))
  
}

filled.contour2 <- function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
            col = color.palette(length(levels) - 1), plot.title, plot.axes, 
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
            axes = TRUE, frame.plot = axes,mar, ...) {
    # modification by Ian Taylor of the filled.contour function
    # to remove the key and facilitate overplotting with contour()
    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
      stop("increasing 'x' and 'y' values expected")
    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    on.exit(par(par.orig))
    w <- (3 + mar.orig[2]) * par("csi") * 2.54
    par(las = las)
    mar <- mar.orig
    plot.new()
    par(mar=mar)
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
      stop("no proper 'z' matrix specified")
    if (!is.double(z)) 
      storage.mode(z) <- "double"
    filled.contour(as.double(x), as.double(y), z, as.double(levels), 
                   col = col)
    if (missing(plot.axes)) {
      if (axes) {
        title(main = "", xlab = "", ylab = "")
        Axis(x, side = 1)
        Axis(y, side = 2)
      }
    }
    else plot.axes
    if (frame.plot) 
      box()
    if (missing(plot.title)) 
      title(...)
    else plot.title
    invisible()
  }

average_over_modules = function(input_file, modules_file, output_file){
  
  Conn_info = readMat(input_file)
  labels = unlist(Conn_info$labels)
  
  subjects = unlist(Conn_info$subjects)
  X = Conn_info$merged.matrices
  merged.matrices.mat = NULL
  n_subjects = length(subjects)
  for (i in seq(n_subjects)){
    FC= X[ , , i]
    average.FC = average_params(FC, modules_file, ICN = T)
    merged.matrices.mat = rbind(merged.matrices.mat, c(squareform2(average.FC), diag(average.FC)) )
    
  }
  
  valid.indices = seq(nrow(average.FC))
  writeMat(con = output_file, subjects=subjects, labels=labels, merged.matrices.mat=merged.matrices.mat, valid.indices = valid.indices)
  
}

average_params = function(mat, modules_file, ICN){
  
  load(modules_file)
  if (ICN) weights = as.matrix(ICNweights)
  else weights = as.matrix(cogweights.all)
  
  n = dim(weights)[1]
  if (!ICN) {
    average = c()
    for (j in seq(n)){
      Id = matrix(1, length(valid.indices), length(valid.indices)) 
      Id = Id - diag(diag(Id))
      average[j] = t(weights[j, valid.indices]) %*% mat[valid.indices, valid.indices] %*% weights[j, valid.indices] /(t(weights[j, valid.indices]) %*% Id %*% weights[j, valid.indices])
    }
    
    names(average) = cluster_names
  } else {
    average = matrix(NA, n, n)
    for (i in seq(n)){
      for (j in seq(n)){
        Id = matrix(1, length(valid.indices), length(valid.indices))
        if (i == j) Id = Id - diag(diag(Id))
        average[i, j] = t(weights[i, valid.indices]) %*% mat[valid.indices, valid.indices] %*% weights[j, valid.indices] /(t(weights[i, valid.indices]) %*% Id %*% weights[j, valid.indices])
      }
    }
    rownames(average) = colnames(average) = rownames(weights)
  }
  
  return(average)
}

average_params_variability = function(values, modules_file, ICN){
  
  nsamples = nrow(values)
  load(modules_file)
  
  if (ICN) weights = as.matrix(ICNweights)
  else weights = as.matrix(cogweights.all)
  
  n = nrow(weights)
  
  doit.1 = function(values, valid.indices, modules_file, weights){
    n = nrow(weights)
    
    mat = get_adj(as.numeric(values), modules_file, sort=FALSE)
    average = rep(NA, n)
    
    for (j in seq(n)){
      
      Id = matrix(1, length(valid.indices), length(valid.indices)) 
      Id = Id - diag(diag(Id))
      average[j] = t(weights[j, valid.indices]) %*% mat[valid.indices, valid.indices] %*% weights[j, valid.indices] /(t(weights[j, valid.indices]) %*% Id %*% weights[j, valid.indices])
    }
    return(average)
  }
  
  
  doit.2 = function(values, valid.indices, modules_file, weights){
    n = nrow(weights)
    mat = get_adj(as.numeric(values), modules_file, sort=FALSE)
    average = matrix(NA, n, n)
    
    for (i in seq(n)){
      for (j in seq(n)){
        Id = matrix(1, length(valid.indices), length(valid.indices))
        if (i == j) Id = Id - diag(diag(Id))
        average[i, j] = t(weights[i, valid.indices]) %*% mat[valid.indices, valid.indices] %*% weights[j, valid.indices] /(t(weights[i, valid.indices]) %*% Id %*% weights[j, valid.indices])
      }
    }
    diag.average = diag(average)
    return(c(squareform2(average), diag.average))
  }
  
  if (!ICN) {
    averages = foreach(i = seq(nsamples), .combine=rbind) %dopar% doit.1(values[i, ], valid.indices, modules_file, weights)
    
    means = colMeans(averages)
    quants = apply(averages, 2, quantile,  c(.025, .975))
    quantsl = quants[1, ]
    quantsh = quants[2, ]
    names(means) = rownames(weights)
    names(quantsl) = rownames(weights)
    names(quantsh) = rownames(weights)
    
  } else {
    averages = foreach(i = seq(nsamples), .combine=rbind) %dopar% doit.2(values[i, ], valid.indices, modules_file, weights)
    
    l = n*(n-1)*0.5
    m = colMeans(averages)
    means = squareform(m[1:l])
    diag(means) = m[(l+1):length(m)] 
    
    quants = apply(averages, 2, quantile,  c(.025, .975))
    
    quantsl = squareform(quants[1, 1:l])
    quantsh = squareform(quants[2, 1:l])
    diag(quantsl) = quants[1, (l+1):length(m)]
    diag(quantsh) = quants[2, (l+1):length(m)]
    
    rownames(means) = colnames(means) = rownames(weights)
    rownames(quantsl) = colnames(quantsl) = rownames(weights)
    rownames(quantsh) = colnames(quantsh) = rownames(weights)
    
  }
  
  return(list(means = means, quantsl = quantsl, quantsh = quantsh))
}

average_params_variability_modules = function(values, modules_file, ICN){
  
  nsamples = nrow(values)
  load(modules_file)
  
  weights = as.matrix(cogICNs)
  
  
  doit.cog = function(values, modules_file, weights){
    nrows = nrow(weights)
    ncols = ncol(weights)
    mat = get_adj_modules(as.numeric(values), modules_file)
    average = rep(NA, nrows)
    
    for (j in seq(nrows)){
      
      Id = matrix(1, ncols, ncols)
      Id = Id - diag(diag(Id))
      average[j] = t(weights[j, ]) %*% mat %*% weights[j, ] /(t(weights[j, ]) %*% Id %*% weights[j, ])
    }
    return(average)
  }
  
  
  if (!ICN) {
    averages = foreach(i = seq(nsamples), .combine=rbind) %dopar% doit.cog(values[i, ], modules_file, weights)
    
    means = colMeans(averages)
    quants = apply(averages, 2, quantile,  c(.025, .975))
    quantsl = quants[1, ]
    quantsh = quants[2, ]
    names(means) = rownames(weights)
    names(quantsl) = rownames(weights)
    names(quantsh) = rownames(weights)
    
  } else {
    # the average in this case was already done prior to fitting the model
    averages = values 
    n = nrow(ICNweights)
    l = n*(n-1)*0.5
    m = colMeans(averages)
    means = squareform(m[1:l])
    diag(means) = m[(l+1):length(m)] 
    
    quants = apply(averages, 2, quantile,  c(.025, .975))
    
    quantsl = squareform(quants[1, 1:l])
    quantsh = squareform(quants[2, 1:l])
    diag(quantsl) = quants[1, (l+1):length(m)]
    diag(quantsh) = quants[2, (l+1):length(m)]
    
    rownames(means) = colnames(means) = rownames(ICNweights)
    rownames(quantsl) = colnames(quantsl) = rownames(ICNweights)
    rownames(quantsh) = colnames(quantsh) = rownames(ICNweights)
    
  }
  
  return(list(means = means, quantsl = quantsl, quantsh = quantsh))
}

plot_matrix = function(adj, lim=NULL, pal=MYPALETTE){ #colorspace::diverge_hsv(15, power=2)
  
  if (is.null(rownames(adj))) adj_names = seq(nrow(adj))
  else adj_names = rownames(adj)
  
  
  seq_names = seq(length(adj_names))
  rownames(adj) = colnames(adj) = seq_names
  molten_adj = melt(adj)
  plot_data = rename(molten_adj, from = X1, to = X2) 
  
  theme_new = theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(), 
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),       
    legend.text = element_text(size=GGTEXTSIZE2),
    legend.title = element_blank(), 
    aspect.ratio = 1)
  
  if (is.null(lim)) {
    lim = max(abs(adj), na.rm=TRUE)
  }
  
  if (lim == -1) {
    myplot <- ggplot(plot_data, aes(x = from, y = to, fill = factor(value))) +
      geom_raster() +
      theme_bw() +
      theme_new + 
      # scale_y_reverse() +
      scale_fill_manual(values=pal) 
  }
  else 
  {
    myplot <- ggplot(plot_data, aes(x = from, y = to, fill = value)) +
      geom_raster() +
      theme_bw() +
      theme_new + 
      # scale_y_reverse() +
      scale_fill_gradientn(colours = pal, 
                           limits=c(-lim, lim), na.value = "white")
  }
  
  for (i in seq_names){
    myplot <- myplot + 
      annotate("text", x = 0, y = i, label = adj_names[i], size = GGTEXTSIZE3) +
      annotate("text", y = 0, x = i, label = adj_names[i], size = GGTEXTSIZE3) +
      annotate("text", x = length(seq_names) + 1, y = i, label = adj_names[i], size = GGTEXTSIZE3) +
      annotate("text", y = length(seq_names) + 1, x = i, label = adj_names[i], size = GGTEXTSIZE3) 
    
  }
  
  print(myplot)
}

list_connections = function(adj, n=NULL, incdiag=TRUE, absord=FALSE){
  rows = nrow(adj)
  d = diag(adj)
  vals = squareform2(adj)
  
  rnames = rownames(adj)
  if (is.null(rnames)) rnames = seq(rows)
  
  tags = matrix(rep(rnames, rows), rows, rows)
  from = squareform2(t(tags))
  to = squareform2(tags)
  
  if (incdiag) {
    vals = c(d, vals)
    from = c(rnames, from)
    to = c(rnames, to)
  }
  
  if (absord) ord = order(abs(vals), decreasing=TRUE) 
  else ord = order(vals, decreasing=TRUE)
  
  from = from[ord]
  to = to[ord]
  vals = vals[ord]
  if (is.null(n)) n = length(vals)
  
  return(data.frame(from=from[1:n], to=to[1:n], value=vals[1:n]))
  
}

compute_modules = function(WD, artifacts, ICNnames, figname=NULL){
  setwd(WD)
  
  # get ROI sizes
  
  roi_volumes = as.numeric(read.table(ROI_VOL_FILE))
  
  # get BRAINMAP ICN weights
  ICNs = read.table("ICNs.csv")
  rownames(ICNs) = ICNnames
  colnames(ICNs) = labels$Label
  ICNs = ICNs[-artifacts, ]
  ICNnames = ICNnames[-artifacts]
  ICNweights = ICNs  
  save_fig()
  image(as.matrix(ICNweights))
  
  #compute partition by ICNs
  partition = apply(ICNweights, 2, which.max)
  partition[colSums(ICNweights)==0] = 0 
  names(partition) = formatC(partition, width=2, flag="0")
  #print(partition)
  
  # get cognitive BRAINMAP weights and select behavioural domains
  cogICNs  = read.table("cognitive.txt", sep="\t", header=F, row.names=1)
  if(ONLY_BEHAV) {
    inc = grep("Action\\.|Cognition\\.|Emotion\\.|Interoception\\.|Perception\\.", rownames(cogICNs))
    cogICNs = cogICNs[inc, ]
  }
  
  cogICNs = cogICNs[, -artifacts]
  
  # hierarchical clustering of behavioural domains
#  cogICNs.aux = cogICNs 
#  rownames(cogICNs.aux) = paste0(seq(nrow(cogICNs.aux)), ' - ', rownames(cogICNs.aux) )
  d.cor = as.dist(1-cor(t(cogICNs)))
  
  dendro.cor = hclust(d.cor, method="average")
  
  clusters = cutree(dendro.cor, h = HCTHRESH) 
  
  cluster_labels = list()
  cogICNs.red = matrix(0, max(clusters), ncol(cogICNs))
  colnames(cogICNs) = colnames(cogICNs.red) = ICNnames
  
  for (i in unique(clusters)){
    cogICNs.red[i, ] = colMeans(cogICNs[clusters == i, ])
    cluster_labels[i] = paste(rownames(cogICNs)[clusters == i], collapse=", ")
  }
  cluster_names = paste0('BEHAV_', unique(clusters))
  
  clusters = data.frame(Name=rownames(cogICNs), Cluster = clusters, ClusterName=cluster_names[clusters])
  cogICNs = cogICNs[order(clusters$Cluster), ]
  
  clusters = clusters[order(clusters$Cluster), ]
  
  # plot dendrogram 
  cogICNs.aux = cogICNs[order(clusters$Cluster), ] 
  rownames(cogICNs.aux) = paste0(seq(nrow(cogICNs.aux)), '(', clusters$Cluster, ') - ',rownames(cogICNs.aux) )
  
  d.cor = as.dist(1-cor(t(cogICNs.aux)))
  dendro.cor = hclust(d.cor, method="average")
  
  save_fig(figname=figname, res=BWRES)
  plot(dendro.cor, cex=CEX_LAB, ylab="", xlab="", sub="", lwd=3, main ="")
  abline(HCTHRESH, 0, lty = 2, col="red", lwd=3)
  
  write.csv(clusters, file = "cog_clusters.txt", row.names=FALSE)
  
  # find cognitive correspondence of ROIs
  rownames(cogICNs.red) = cluster_names
  cogweights = as.matrix(cogICNs.red) %*% as.matrix(ICNs)
  
  # clean weights
  y = x = unlist(ICNweights)
  x = x[x!=0]
  ICNweights.clust = Mclust(x, G=2)
  
  #print(length(ICNweights.clust$classification))
  y[y!=0] = ICNweights.clust$classification
  
  breaks= seq(-100, 100, .1)
  save_fig(figname = 'ICN_weights')
  par(mfrow=c(1, 2))
  plot(unlist(ICNweights), y)
  hist(unlist(ICNweights), breaks=breaks, col=rgb(1,0,0,.5), xlim=c(-10,10), ylim=c(0,1000))
  ICNweights[ matrix(y==which.min(ICNweights.clust$parameters$mean), nrow(ICNweights), ncol(ICNweights))] = 0
  ICNweights[ICNweights<0] = 0
  print(ICNweights.clust$parameters)
  hist(unlist(ICNweights), breaks=breaks, col=rgb(0,0,1,.5), add=T)
  
  
  # finally, normalize the weights to take into account ROI size
  ICNweights = t(t(ICNweights)*roi_volumes/sum(roi_volumes))
  
  
  save_fig()
  par(mfrow = c(2, 2))
  image(as.matrix(t(cogICNs[dendro.cor$order, ])), main="Ordered weights")
  image(cor(t(cogICNs[dendro.cor$order, ])), , main="Ordered weight correlations")
  image(as.matrix(cogweights), main = "Cognitive ROI weights")
  image(as.matrix(cor(t(cogICNs.red))), main = "Cluster correlations")
  
  # clean cogweights
  cogweights.all = as.matrix(cogICNs) %*% as.matrix(ICNs)
  y = x = as.vector(cogweights.all)
  x = x[x!=0]
  cogweights.clust = Mclust(x, G=2)
  y[y!=0] = cogweights.clust$classification
  
  save_fig(figname = 'CogROI_weights')
  par(mfrow=c(1, 2))
  plot(as.numeric(cogweights.all), y)
  hist(as.numeric(cogweights.all), breaks=breaks, col=rgb(1,0,0,.5), xlim=c(-10,10))
  
  cogweights.all[y==which.min(cogweights.clust$parameters$mean)] = 0
  cogweights.all[cogweights.all<0] = 0
  print(cogweights.clust$parameters)
  hist(unlist(cogweights.all), breaks=breaks, col=rgb(0,0,1,.5), add=T)
  
  # clean cogICNs
  save_fig(figname = 'CogICN_weights')
  y = x = unlist(cogICNs)
  x = x[x!=0]
  cogICNs.clust = Mclust(x, G=2)
  y[y!=0] = cogICNs.clust$classification
  par(mfrow=c(1, 2))
  plot(unlist(cogICNs), y)
  hist(unlist(cogICNs), breaks=breaks, col=rgb(1,0,0,.5), xlim=c(-1,1))
  
  cogICNs[ matrix(y==which.min(cogICNs.clust$parameters$mean), nrow(cogICNs), ncol(cogICNs))] = 0
  cogICNs[cogICNs<0] = 0
  hist(unlist(cogICNs), breaks=breaks, col=rgb(0,0,1,.5), add=T)
  
  module_names = sort(unique(names(partition)))
  save(partition, cogweights, n.rois, ICNweights, valid.indices, ICNnames, cluster_labels, cluster_names, cogICNs, cogICNs.red, cogweights.all, clusters, module_names, file = MODULES_FILE )
}

squareform2 = function (x) 
{
  if (is.vector(x)) {
    n <- length(x)
    m <- floor(sqrt(2 * n))
    if (m * (m + 1) != 2 * n) 
      stop("Argument 'x' does not correspond to a distance matrix.")
    inds <- c()
    k <- m + 1
    for (i in 1:(k - 1)) inds <- c(inds, (1 + i + (i - 1) * 
                                            k):(i * k))
    y <- numeric(k * k)
    y[inds] <- x
    y <- matrix(y, k, k) + t(matrix(y, k, k))
  }
  else if (is.matrix(x)) {
    m <- nrow(x)
    n <- ncol(x)
    if (m != n) 
      stop("Argument 'x' must be a vector or a square matrix.")
    if (n == 1) 
      return(c())
    inds <- c()
    for (i in 1:(n - 1)) inds <- c(inds, (1 + i + (i - 1) * 
                                            n):(i * n))
    y <- x[inds]
  }
  return(y)
}

wideScreen <- function(howWide=Sys.getenv("COLUMNS")) {
  options(width=as.integer(howWide))
}

flatten  = function(arr, indices) {
  n = dim(arr)[3]
  l = length(indices)
  mat = matrix(NA, n, l*(l-1)/2)
  for (i in seq(n)) mat[i, ] = squareform(arr[indices, indices, i])
  
  return(mat)
}

plot_intra_values.slopes = function(average.ICN, MODULES_FILE, scatter=TRUE, plot_inter=TRUE){
  
  load(MODULES_FILE)

  slope_muc.means = diag(average.ICN$slope_muc$means)
  log_slope_sigmac.means = diag(average.ICN$log_slope_sigmac$means)
  slope_muc.quantsl = diag(average.ICN$slope_muc$quantsl)
  log_slope_sigmac.quantsl = diag(average.ICN$log_slope_sigmac$quantsl)
  slope_muc.quantsh = diag(average.ICN$slope_muc$quantsh)
  log_slope_sigmac.quantsh = diag(average.ICN$log_slope_sigmac$quantsh)
  N = length(slope_muc.means)

  
  print(cor.test( squareform2(average.ICN$slope_muc$means), squareform2(average.ICN$log_slope_sigmac$means)))  
  print(cor.test( slope_muc.means, log_slope_sigmac.means))  
  
    
  muc_order = order(slope_muc.means, decreasing=FALSE)
  sigmac_order = order(log_slope_sigmac.means, decreasing=TRUE)

    ylimits = range(c(log_slope_sigmac.quantsl, log_slope_sigmac.quantsh))
    xlimits = range(c(slope_muc.quantsl, slope_muc.quantsh))

  if (scatter){
    
    par(mfrow=c(1, 1), mar=c(8,8,5,5), mgp=c(4,1,0))
    plot(0, xlab = expression(bar(beta)[mu]), ylab = expression(bar(beta)[sigma]), xlim = xlimits, ylim = ylimits, 
         cex.lab= CEX_LAB, cex.axis = CEX_AXIS, type="n")

  if (plot_inter){    
    points(squareform2(average.ICN$slope_muc$means), squareform2(average.ICN$log_slope_sigmac$means), pch = 20, cex = 1)
    }
    
    segments(slope_muc.quantsl, log_slope_sigmac.means, slope_muc.quantsh, log_slope_sigmac.means, lwd=2) 
    segments(slope_muc.means, log_slope_sigmac.quantsl, slope_muc.means, log_slope_sigmac.quantsh, lwd=2) 

    points(jitter(slope_muc.means), jitter(log_slope_sigmac.means), pch = 22, cex = 7, bg="white")
    text(slope_muc.means, log_slope_sigmac.means, labels = ICNnames, cex = CEX_NAMES)
    
  } else {
    par(mfrow=c(2, 1), mar=c(8,20,5,5), mgp=c(4,1,0))
    barplot(slope_muc.means[muc_order], ylab = "", xlab = expression("Average " ~ beta[mu] ), 
            names.arg=ICNnames[muc_order], horiz = TRUE, cex.lab= CEX_LAB, cex.axis = CEX_AXIS, cex.names = CEX_NAMES, las = 1, width = 1, space=0, xlim=c(-0.003, 0.004), col="gray90") 
    segments(slope_muc.quantsl[muc_order], seq(N)-0.5, slope_muc.quantsh[muc_order],  seq(N)-0.5, lwd=2) 
    
    barplot(log_slope_sigmac.means[sigmac_order], ylab = "", xlab = expression("Average " ~ beta[sigma]),  
            names.arg=ICNnames[sigmac_order], horiz = TRUE, cex.lab= CEX_LAB, cex.axis = CEX_AXIS, cex.names = CEX_NAMES, las = 1, width = 1, space=0, xlim=c(0.004, 0.0065), xpd=FALSE, col="gray90") 
    segments(log_slope_sigmac.quantsl[sigmac_order], seq(N)-0.5, log_slope_sigmac.quantsh[sigmac_order],  seq(N)-0.5, lwd=2) 
  }
}

plot_intra_values.intercept = function(average.ICN, MODULES_FILE, scatter=TRUE){
  
  load(MODULES_FILE)
  
  slope_muc.means = diag(average.ICN$slope_muc$means)
  intercept_muc.means = diag(average.ICN$intercept_muc$means)
  slope_muc.quantsl = diag(average.ICN$slope_muc$quantsl)
  intercept_muc.quantsl = diag(average.ICN$intercept_muc$quantsl)
  slope_muc.quantsh = diag(average.ICN$slope_muc$quantsh)
  intercept_muc.quantsh = diag(average.ICN$intercept_muc$quantsh)
  N = length(slope_muc.means)
  
  muc_order = order(slope_muc.means, decreasing=FALSE)
  sigmac_order = order(intercept_muc.means, decreasing=TRUE)
  
  if (scatter){
    
    par(mfrow=c(1, 1), mar=c(8,8,5,5), mgp=c(4,1,0))
    plot(0, xlab = expression("Average " ~ beta[mu]), ylab = expression("Average " ~ alpha[mu]), xlim = c(-0.0025, 0.002), ylim = c(0.2, 0.5), 
         cex.lab= CEX_LAB, cex.axis = CEX_AXIS, type="n")
    segments(slope_muc.quantsl, intercept_muc.means, slope_muc.quantsh, intercept_muc.means, lwd=2) 
    segments(slope_muc.means, intercept_muc.quantsl, slope_muc.means, intercept_muc.quantsh, lwd=2) 
    
    points(slope_muc.means, intercept_muc.means, pch = 22, cex = 5, bg="white")
    
    
    text(slope_muc.means, intercept_muc.means, labels = ICNnames, cex = CEX_NAMES)
    
  } else {
    par(mfrow=c(2, 1), mar=c(8,20,5,5), mgp=c(4,1,0))
    barplot(slope_muc.means[muc_order], ylab = "", xlab = expression("Average " ~ beta[mu] ), 
            names.arg=ICNnames[muc_order], horiz = TRUE, cex.lab= CEX_LAB, cex.axis = CEX_AXIS, cex.names = CEX_NAMES, las = 1, width = 1, space=0, xlim=c(-0.003, 0.004), col="gray90") 
    segments(slope_muc.quantsl[muc_order], seq(N)-0.5, slope_muc.quantsh[muc_order],  seq(N)-0.5, lwd=2) 
    
    barplot(intercept_muc.means[sigmac_order], ylab = "", xlab = expression("Average " ~ alpha[mu]),  
            names.arg=ICNnames[sigmac_order], horiz = TRUE, cex.lab= CEX_LAB, cex.axis = CEX_AXIS, cex.names = CEX_NAMES, las = 1, width = 1, space=0, xlim=c(0.004, 0.0065), xpd=FALSE, col="gray90") 
    segments(intercept_muc.quantsl[sigmac_order], seq(N)-0.5, intercept_muc.quantsh[sigmac_order],  seq(N)-0.5, lwd=2) 
  }
}

plot_cognitive_domains = function(average.cog, MODULES_FILE, bycluster=FALSE, scatter=TRUE){
  load(MODULES_FILE)
  
  N = length(average.cog$slope_muc$means)
  clusters$slope_muc = average.cog$slope_muc$means
  clusters$log_slope_sigmac = average.cog$log_slope_sigmac$means
  
  slope_muc_clusters = aggregate(slope_muc ~ Cluster, FUN=mean, data=clusters)$slope_muc
  log_slope_sigmac_clusters = aggregate(log_slope_sigmac ~ Cluster, FUN=mean, data=clusters)$log_slope_sigmac
  
  clusters$cluster_slope_muc = slope_muc_clusters[clusters$Cluster]
  clusters$cluster_log_slope_sigmac = log_slope_sigmac_clusters[clusters$Cluster]
  
  if (bycluster){
    muc_order = order(clusters$cluster_slope_muc, decreasing=FALSE)
    sigmac_order = order(clusters$cluster_log_slope_sigmac, decreasing=TRUE)
    cognames = paste0(clusters$Cluster, ". ", clusters$Name)
    set.seed(1)
    colors = rep(c("red", "cyan", "green", "yellow", "orange"), 4)[clusters$Cluster] #sample(rainbow(max(clusters$Cluster)))[clusters$Cluster]
    pchs = c(rep(21, 5), rep(22, 5), rep(23, 5), rep(25, 5))[clusters$Cluster] #sample(seq(21, 25), max(clusters$Cluster), replace=T)[clusters$Cluster]
    colors_muc = colors[muc_order]
    colors_sigmac = colors[sigmac_order]
    
  } else {
    muc_order = order(clusters$slope_muc, decreasing=FALSE)
    sigmac_order = order(clusters$log_slope_sigmac, decreasing=TRUE)
    cognames = clusters$Name
    colors = "white"
    colors_muc = colors_sigmac = "gray90"
    pchs = 22
  }
  
  
  if (scatter){
    ylimits = range(c(average.cog$log_slope_sigmac$quantsl, average.cog$log_slope_sigmac$quantsh))
    xlimits = range(c(average.cog$slope_muc$quantsl, average.cog$slope_muc$quantsh))
    
    par(mfrow=c(1, 1), mar=c(8,8,5,5), mgp=c(4,1,0))
    plot(0, xlab = expression("Average " ~ hat(beta)[mu]), ylab = expression("Average " ~ hat(beta)[sigma]),  xlim = xlimits, ylim = ylimits, 
         cex.lab= CEX_LAB, cex.axis = CEX_AXIS, type="n")
    segments(average.cog$slope_muc$quantsl, average.cog$log_slope_sigmac$means, average.cog$slope_muc$quantsh, average.cog$log_slope_sigmac$means, lwd=2) 
    segments(average.cog$slope_muc$means, average.cog$log_slope_sigmac$quantsl, average.cog$slope_muc$means, average.cog$log_slope_sigmac$quantsh, lwd=2)
    
    points(average.cog$slope_muc$means, average.cog$log_slope_sigmac$means, cex = 5, bg = colors, pch = pchs)
    text(average.cog$slope_muc$means, average.cog$log_slope_sigmac$means, labels = seq(length(cognames)), cex = CEX_NAMES)
    
  } else {
    par(mfrow=c(2, 1), mar=c(8,20,5,5), mgp=c(4,1,0))
    
    barplot(average.cog$slope_muc$means[muc_order], ylab = "", xlab = expression("Average " ~ beta[mu] ), 
            names.arg=cognames[muc_order], horiz = TRUE, cex.lab= CEX_LAB, cex.axis = CEX_AXIS, las = 1, width = 1, space=0, xlim=c(-0.002, 0.001), col=colors_muc) 
    segments(average.cog$slope_muc$quantsl[muc_order], seq(N)-0.5, average.cog$slope_muc$quantsh[muc_order],  seq(N)-0.5, lwd=2) 
    
    barplot(average.cog$log_slope_sigmac$means[sigmac_order], ylab = "", xlab = expression("Average " ~ beta[sigma]),  
            names.arg=cognames[sigmac_order], horiz = TRUE, cex.lab= CEX_LAB, cex.axis = CEX_AXIS, las = 1, width = 1, space=0, xlim=c(0.004, 0.006), xpd=FALSE, col=colors_sigmac) 
    segments(average.cog$log_slope_sigmac$quantsl[sigmac_order], seq(N)-0.5, average.cog$log_slope_sigmac$quantsh[sigmac_order],  seq(N)-0.5, lwd=2) 
  }
}

tabulate_demo = function(demo, INPUT_FILE=NULL){
  
  print(INPUT_FILE)
  
  if (!is.null(INPUT_FILE)) {
    info = readMat(INPUT_FILE)
    indices = match(unlist(info$subjects), demo$Subject)
    demo = demo[indices, ]
  }
  print(summarise(group_by(demo, class), n()))
  print(summarise(group_by(demo, class), mean(age)))
  print(summarise(group_by(demo, class), std(age)))
  print(summarise(group_by(demo, class, sex), n()))
}

plot_pies = function(adj, modules_file, pal=MYPALETTE){
  
  load(modules_file)
  
  partition = partition[valid.indices]
  adj = adj[valid.indices, valid.indices]
  
  ord = order(partition)
  partition.ord = partition[ord]
  
  adj = adj[ord, ord]
  
  adj = adj[partition.ord>0, partition.ord>0]
  partition.ord = partition.ord[partition.ord>0]
  rownames(adj) = colnames(adj) = seq(nrow(adj))
  
  n_modules = length(unique(partition.ord))
  # count connections
  plot_data = c()
  for (r in seq(length(pal))){
    value = matrix(0, n_modules, n_modules)
    
    for (i in seq(n_modules)){
      for (j in seq(n_modules)){
        
        if (i==j) {
          x = squareform2(adj[partition.ord == i, partition.ord == j])==r
        } else {
          x = adj[partition.ord == i, partition.ord == j]==r
        }
        
        value[i, j] = sum(x)/length(x)
        
      }
    }
    
    molten = melt(value)
    molten$group = as.factor(r)
    
    plot_data = rbind(plot_data, molten)
  }
  plot_data$from = factor(module_names[plot_data$X1], levels = module_names)
  plot_data$to = factor(module_names[plot_data$X2], levels = module_names)
  
  
  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      axis.text.x=element_blank(), 
      axis.text.y=element_blank(),
      legend.position="none", 
      strip.text.x = element_text(size = GGTEXTSIZE4, angle = 90),
      strip.text.y = element_text(size = GGTEXTSIZE4, angle = 180),
      plot.margin = unit(c(2.5,.5,.5,.5),"cm")
#      plot.margin = unit(3, "points")
      
    )
  
  pies = ggplot(plot_data, aes(x="", y=value, fill=group)) +
    geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0) + blank_theme + 
    scale_fill_manual(values=pal) + facet_grid(from ~ to, as.table=FALSE, switch = "both") 
     #+ margin(8, 8, 1, 1)
  
  print(pies)
  
}

analyze_motion = function(motion_dir, demo, figname){
  
  params = read.table(paste0(motion_dir, 'HeadMotion.csv'), sep=',', header =T)
  
  params$mean.FD_Power_scrubbed = t(read.table(paste0(motion_dir, 'FD_Power_scrubbed.csv'), sep=','))
  params = merge(params, demo, by.x="Subject.ID", by.y="Subject" )
  
  save_fig(figname)
  par(mfrow=c(1, 3))
  boxplot( mean.FD_Power ~ class, data = params, ylim= c(0, 0.5))
  legend("topright", legend = t.test( mean.FD_Power ~ class, data = params)$p.value)
  legend("bottomright", legend = t.test( mean.FD_Power ~ class, data = params)$statistic)
  
  params = subset(params, !is.nan(mean.FD_Power_scrubbed))
  boxplot( mean.FD_Power ~ class, data = params, ylim= c(0, 0.5))
  legend("topright", legend = t.test( mean.FD_Power ~ class, data = params)$p.value)
  legend("bottomright", legend = t.test( mean.FD_Power ~ class, data = params)$statistic)
  
    boxplot( mean.FD_Power_scrubbed ~ class, data = params, ylim= c(0, 0.5))
  legend("topright", legend = t.test( mean.FD_Power_scrubbed ~ class, data = params)$p.value)
  legend("bottomright", legend = t.test( mean.FD_Power_scrubbed ~ class, data = params)$statistic)
  
  print(dim(params))
  
  return( params )
  
}

do_crossvalidate = function(fold, data){
  X = data$X
  y = data$y
  X.train = X[-fold, , drop=FALSE]
  X.test = X[fold, , drop=FALSE]
  y.train = y[-fold]  
  y.test = y[fold]
  
  mod = train(X.train, y.train,  method = method, trControl = ctrl, tuneGrid=tuneGrid) 
  #print(summary(mod))
  pred = predict(mod, X.test)
  RMSE = sqrt(mean((pred - y.test)^2 ))
  
  return(list(pred=pred, RMSE=RMSE, y.test = y.test, model=mod))
}

get_nspn_data = function(NSPN_FILE, modules_file){
  nspn_info = readMat(NSPN_FILE)
  load(modules_file)
  FC.nspn.orig = data.frame(nspn_info$merged.matrices.mat)
  FC.nspn = matrix(NA, nrow(FC.nspn.orig), (length(valid.indices)-1)*length(valid.indices)/2)
  n_rois = dim(nspn_info$merged.matrices)[1]
  
  # sync nspn data with DAD's
  for (i in seq(nrow(FC.nspn))){
    print(i)
    
    M = matrix(NA, n_rois, n_rois)
    M[nspn_info$valid.indices, nspn_info$valid.indices] = squareform(as.numeric(FC.nspn.orig[i, ]))
    FC.nspn[i, ] = squareform2(M[valid.indices, valid.indices])
  }
  
  return(FC.nspn)
}

corr_slope_distance = function(x, adj, dist.mat){ 
  slope_muc.degree = compute_degree_kernel(adj, dist.mat, x[1], x[2])
  #print(x)
  cor.test(slope_muc.degree, slope_muc.PET)$estimate
}






