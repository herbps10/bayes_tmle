data {
  int<lower=0, upper=1> hierarchical;
  int<lower=0> N;
  int<lower=0> G;
  array[N] int<lower=1, upper=G> group;
  array[G] int<lower=1, upper=N> group_start;
  array[G] int<lower=1, upper=N> group_end;
  array[N] int<lower=0, upper=1> y;
  array[G] int<lower=1, upper=N> N_group;
  vector[N] Qbar;
  vector[N] Qbar0;
  vector[N] Qbar1;
  
  vector[N] H1;
  vector[N] H1_0;
  vector[N] H1_1;
  
  vector[N] H2;
}
transformed data {
  vector[N] Qw;
  for(g in 1:G) {
    Qw[group_start[g]:group_end[g]] = rep_vector(1.0 / N_group[g], N_group[g]);
  }
}
parameters {
  vector[G] epsilon;
  array[hierarchical] real<lower=-1,upper=1> mu;
  array[hierarchical] real<lower=0> sigma;
}
transformed parameters {
  vector[N] Qbar_fluctuated_logit;
  vector[N] Qw_fluctuated;
  vector[N] Qbar0_fluctuated;
  vector[N] Qbar1_fluctuated;
  vector<lower=-1, upper=1>[G] psi;
  
  vector[N] dQbar;
  vector[N] dQw;
  vector[G] dpsi;
  
  for(g in 1:G) {
    Qw_fluctuated[group_start[g]:group_end[g]] = exp(epsilon[g] .* H2[group_start[g]:group_end[g]]) .* Qw[group_start[g]:group_end[g]];
    Qw_fluctuated[group_start[g]:group_end[g]] = Qw_fluctuated[group_start[g]:group_end[g]] / sum(Qw_fluctuated[group_start[g]:group_end[g]]);
    
    Qbar_fluctuated_logit[group_start[g]:group_end[g]] = logit(Qbar[group_start[g]:group_end[g]]) + epsilon[g] * H1[group_start[g]:group_end[g]];
    
    Qbar0_fluctuated[group_start[g]:group_end[g]] = inv_logit(logit(Qbar0[group_start[g]:group_end[g]]) + epsilon[g] * H1_0[group_start[g]:group_end[g]]);
    Qbar1_fluctuated[group_start[g]:group_end[g]] = inv_logit(logit(Qbar1[group_start[g]:group_end[g]]) + epsilon[g] * H1_1[group_start[g]:group_end[g]]);
    
    psi[g] = sum(Qw_fluctuated[group_start[g]:group_end[g]] .* (Qbar1_fluctuated[group_start[g]:group_end[g]] - Qbar0_fluctuated[group_start[g]:group_end[g]]));
    
    dQbar[group_start[g]:group_end[g]] = H1_1[group_start[g]:group_end[g]] .* Qbar1_fluctuated[group_start[g]:group_end[g]] .* (1 - Qbar1_fluctuated[group_start[g]:group_end[g]]) 
      - H1_0[group_start[g]:group_end[g]] .* Qbar0_fluctuated[group_start[g]:group_end[g]] .* (1 - Qbar0_fluctuated[group_start[g]:group_end[g]]);
      
    dQw[group_start[g]:group_end[g]] = Qw_fluctuated[group_start[g]:group_end[g]] .* (H2[group_start[g]:group_end[g]] - Qw[group_start[g]:group_end[g]] .* H2[group_start[g]:group_end[g]] .* exp(epsilon[g] .* H2[group_start[g]:group_end[g]]) / sum(Qw[group_start[g]:group_end[g]] .* exp(epsilon[g] .* H2[group_start[g]:group_end[g]])));
    
    dpsi[g] = sum(dQw[group_start[g]:group_end[g]] .* (Qbar1_fluctuated[group_start[g]:group_end[g]] - Qbar0_fluctuated[group_start[g]:group_end[g]]) + Qw_fluctuated[group_start[g]:group_end[g]] .* dQbar[group_start[g]:group_end[g]]);
  }
  
}
model {
  // Jacobian adjustment
  for(g in 1:G) {
    target += log(dpsi[g]);
  }
    
  // Hierarchical prior
  if(hierarchical == 1) {
    psi ~ normal(mu[1], sigma[1]);
    mu[1] ~ uniform(-1, 1);
    //sigma[1] ~ inv_gamma(0.001, 0.001);
    //sigma[1] ~ uniform(0, 10);
  }
  else {
    psi ~ uniform(-1, 1);
  }
  
  // Likelihood
  y ~ bernoulli_logit(Qbar_fluctuated_logit);
  target += sum(log(Qw_fluctuated));
}
