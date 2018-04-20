// non-centered parametrization

data {
   int<lower=0> S;          // number of subjects
   int<lower=0> C;          // number of regions
   vector[C] values[S]; // values
   vector[S] age; // age
   matrix[8, 2] priors;
}

parameters {

   real raw_intercept_muc_mu;
   real<lower=0> raw_intercept_muc_sigma;
   real raw_slope_muc_mu;
   real<lower=0> raw_slope_muc_sigma;

   real raw_log_intercept_sigmac_mu;
   real<lower=0> raw_log_intercept_sigmac_sigma;
   real raw_log_slope_sigmac_mu;
   real<lower=0> raw_log_slope_sigmac_sigma;

   vector[C] raw_intercept_muc;
   vector[C] raw_slope_muc;
   vector[C] raw_log_intercept_sigmac;
   vector[C] raw_log_slope_sigmac; 

}


transformed parameters {

   // global hyperparameters
   real intercept_muc_mu;
   real<lower=0> intercept_muc_sigma;
   real slope_muc_mu;
   real<lower=0> slope_muc_sigma;

   real log_intercept_sigmac_mu;
   real<lower=0> log_intercept_sigmac_sigma;
   real log_slope_sigmac_mu;
   real<lower=0> log_slope_sigmac_sigma;

   // connection-specific intercepts and slopes
   vector[C] intercept_muc;
   vector[C] slope_muc;
   vector[C] log_intercept_sigmac;
   vector[C] log_slope_sigmac; 

   vector[C] muc[S];
   vector<lower=0>[C] sigmac[S];

   intercept_muc_mu = priors[1, 1] + raw_intercept_muc_mu * priors[1, 2];
   intercept_muc_sigma = priors[2, 1] + raw_intercept_muc_sigma * priors[2, 2]; 
   slope_muc_mu = priors[3, 1] + raw_slope_muc_mu * priors[3, 2];
   slope_muc_sigma = priors[4, 1] + raw_slope_muc_sigma * priors[4, 2];

   log_intercept_sigmac_mu = priors[5, 1] + raw_log_intercept_sigmac_mu * priors[5, 2]; 
   log_intercept_sigmac_sigma = priors[6, 1] + raw_log_intercept_sigmac_sigma * priors[6, 2];
   log_slope_sigmac_mu = priors[7, 1] + raw_log_slope_sigmac_mu * priors[7, 2];
   log_slope_sigmac_sigma = priors[8, 1] + raw_log_slope_sigmac_sigma * priors[8, 2];

   intercept_muc = intercept_muc_mu + raw_intercept_muc * intercept_muc_sigma;
   slope_muc = slope_muc_mu + raw_slope_muc * slope_muc_sigma;
   log_intercept_sigmac = log_intercept_sigmac_mu + raw_log_intercept_sigmac * log_intercept_sigmac_sigma;
   log_slope_sigmac = log_slope_sigmac_mu + raw_log_slope_sigmac * log_slope_sigmac_sigma;

   for (s in 1:S){
     muc[s] = intercept_muc + slope_muc*age[s];
     sigmac[s] = exp(log_intercept_sigmac + log_slope_sigmac*age[s]);
   }

}


model {

   // priors
   raw_intercept_muc_mu ~ normal(0, 1);
   raw_intercept_muc_sigma ~ normal(0, 1);
   raw_slope_muc_mu ~ normal(0, 1);
   raw_slope_muc_sigma ~ normal(0, 1);

   raw_log_intercept_sigmac_mu ~ normal(0, 1);   
   raw_log_intercept_sigmac_sigma ~ normal(0, 1); 
   raw_log_slope_sigmac_mu ~ normal(0, 1);
   raw_log_slope_sigmac_sigma ~ normal(0, 1); 

   raw_intercept_muc ~ normal(0, 1);
   raw_slope_muc ~ normal(0, 1);
   raw_log_intercept_sigmac ~ normal(0, 1);
   raw_log_slope_sigmac ~ normal(0, 1);

   for (s in 1:S)
     values[s] ~ normal(muc[s], sigmac[s]);
}
