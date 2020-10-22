//Pella-Tomlinson surplus production model 
// CJ Brown 2020-04-08

data{
  int<lower=1> N; 
  vector[N] lnC; 
  vector[N] lnCPUE; 
  vector<lower=0>[N] catches; 
  real logK_mean;
  real<lower=0> logK_sd;
  real<lower=0> sigma_shape;
  real<lower=0> init_fraction;
  real logr_mean;
  real<lower=0> logr_sd;
}

parameters{
  real lnK;
  real lnr;
  real lnq1;
  real lnq2;
  vector[N] mu_raw; //put here so we can use in transformed params
  real<lower=0> sigma;
  real<lower=0> tau;
}

transformed parameters{
 real<lower=0> K;
 real<lower=0> r;
 real<lower=0> Btemp;
 vector<lower=0>[N] B; 
 vector[N] lnCPUE_hat;
 vector[N] mu;
 
 K = exp(lnK);
 r = exp(lnr);
 mu = mu_raw*sigma; // implies mu ~ normal(0, sigma);

 // Population model. But putting it in transformed params we 
 // speed up computations
 B[1] = K * init_fraction;
 
 for (t in 2:N){
   //Assuming m = 1.19 (use 0.19 here as in formula it is m-1)
      Btemp = (B[t-1] + (r/0.19)*B[t-1]*(1-(B[t-1]/K)^0.19) - catches[t-1])*exp(mu[t]);
  
   if (Btemp < 0.001){
      B[t] = 0.001;
   } else {
    B[t] = Btemp;
   }
 }
 
  lnCPUE_hat[1:16] = lnq1 + log(B[1:16]);
  lnCPUE_hat[17:N] = lnq2 + log(B[17:N]);
 
}

model{ 
 // It is most efficient if this section has 
 // only sampling statements. 
  // Priors
  lnK ~ normal(logK_mean, logK_sd);
  lnr ~ normal(logr_mean, logr_sd);
  lnq1 ~ uniform(0.01, 5);
  lnq2 ~ uniform(0.01, 5);
  // sigma ~ exponential(10);
  sigma ~ gamma(sigma_shape, 0.15);
  tau ~ gamma(6, 0.5);
  // tau ~ exponential(0.05);

  // process errors
  mu_raw ~ std_normal();

  // sampling 
  lnCPUE ~ normal(lnCPUE_hat, tau);

}

generated quantities{
  real<lower = 0> MSY;
  real<lower = 0> FMSY;
  real<lower = 0> EMSY;
  real<lower = 0> BMSY;
  vector<lower=0>[N] Brel; 

  FMSY =  (r/0.19) * (1 - (1/1.19));
  EMSY =  FMSY * exp(lnq2);
  BMSY = K * (1.19^(-1/(0.19)));
  MSY = FMSY * BMSY;
  Brel = B/K;
  
}

  

