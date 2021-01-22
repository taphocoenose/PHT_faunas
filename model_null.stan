data{
  
  int<lower=0> n; // Total n with holdout sample
  int<lower=0> tn; // Training n
  int Faunas[n]; // Count of total faunas
  int Archaeol[n]; // Count of archaeological faunas
  int I[n]; // Index for taxon

  int ll_n; // n holdout observations to evaluate log lik
  int ll_i[ll_n]; // Indices for log lik
  
  // Values for priors
  real<lower=0> phi_scale_prior;
  real<lower=0> mu_alpha_prior;
  real<lower=0> mu_beta_prior;
  
}
parameters{
  
  real<lower=0, upper=1> mu;
  real<lower=2> phi;
  vector<lower=0, upper=1>[tn] Arch_prop;
  
}
model{
  
  // Priors for archaeological proportions
  mu ~ beta(mu_alpha_prior, mu_beta_prior);
  phi ~ normal(2, phi_scale_prior) T[2, ];
  Arch_prop ~ beta(mu * phi, (1 - mu) * phi);
  
  // Likelihood for archaeological proportions, training data
  for(i in 1:tn){
    // Update posterior when observations occur
    if(Faunas[i] > 0){
      target += binomial_lpmf(Archaeol[i] | Faunas[i], Arch_prop[i]);
    }
  }
  
}
generated quantities{
  
  // Evaluate log lik for withheld data
  vector[ll_n] log_lik;
  
  for(i in 1:ll_n){
    log_lik[i] = binomial_lpmf(Archaeol[ll_i[i]] | Faunas[ll_i[i]], 
                                                   Arch_prop[I[ll_i[i]]]);
  }
}
