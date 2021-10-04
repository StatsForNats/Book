data{
  int<lower=1> N;
  int<lower=1> N_pj;
  real lfr[N];
  real lgh[N];
  int pj[N];
}
parameters{
  real a;
  vector[N_pj] p;
  real b;
  real<lower=0> sigmap;
  real<lower=0> sigma;
}
model{
  vector[N] mu;
  sigma ~ cauchy( 0 , 1 );
  sigmap ~ cauchy( 0 , 1 );
  b ~ normal( 0 , 10 );
  p ~ normal( 0 , sigmap );
  a ~ normal( 0 , 100 );
  for ( i in 1:N ) {
    mu[i] = a + p[pj[i]] + b * lgh[i];
  }
  lfr ~ normal( mu , sigma );
}
generated quantities{
  vector[N] mu;
  real dev;
  vector[N] log_lik;
  dev = 0;
    for ( i in 1:N ) {
    mu[i] = a + p[pj[i]] + b * lgh[i];
  }
  dev = dev + (-2)*normal_lpdf( lfr | mu , sigma );
  for ( i in 1:N ) {
    mu[i] = a + p[pj[i]] + b * lgh[i];
    log_lik[i]= normal_lpdf( lfr[i] | mu[i] , sigma);
  }
}
