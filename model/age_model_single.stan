// non-hierarchical, single region/edge at a time 

data {
   int<lower=0> S;          // number of subjects
   vector[S] values; // values
   vector[S] age; // age
   matrix[4, 2] priors;

}

parameters {

   real intercept_muc;
   real slope_muc;
   real log_intercept_sigmac;
   real log_slope_sigmac; 

}


transformed parameters {

   vector[S] muc;
   vector<lower=0>[S] sigmac;

   muc = intercept_muc + slope_muc*age;
   sigmac = exp(log_intercept_sigmac + log_slope_sigmac*age);
  
}


model {
  // priors to have a MAP model 
  
   intercept_muc ~ normal(priors[1, 1], priors[1, 2]);
   slope_muc ~ normal(priors[2, 1], priors[2, 2]);
   log_intercept_sigmac ~ normal(priors[3, 1], priors[3, 2]);
   log_slope_sigmac ~ normal(priors[4, 1], priors[4, 2]);
   values ~ normal(muc, sigmac);
   
}
