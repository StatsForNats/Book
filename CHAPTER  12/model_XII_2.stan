data{
  int<lower=1> N;
  int<lower=1> N_pj;
  real lfr[N];
  real lgh_s[N];
  real tsf_s[N];
  int pj[N];
}
parameters{
  real a;
  real b;
  real c;
  real cc;
  vector[N_pj] p;
  real<lower=0> sigmap;
  real<lower=0> sigma;
}
model{
  vector[N] mu;
  sigma ~ cauchy( 0 , 1 );
  sigmap ~ cauchy( 0 , 1 );
  p ~ normal( 0 , sigmap );
  cc ~ normal( 0 , 1 );
  c ~ normal( 0 , 1 );
  b ~ normal( 0 , 1 );
  a ~ normal( 0 , 50 );
  for ( i in 1:N ) {
    mu[i] = a + p[pj[i]] + b * lgh_s[i] + c * tsf_s[i] + cc * tsf_s[i] * lgh_s[i];
  }
  lfr ~ normal( mu , sigma );
}
generated quantities{
  vector[N] mu;
  vector[N] log_lik;
  real dev;
  dev = 0;
    for ( i in 1:N ) {
    mu[i] = a + p[pj[i]] + b * lgh_s[i] + c * tsf_s[i] + cc * tsf_s[i] * lgh_s[i];
  }
  dev = dev + (-2)*normal_lpdf( lfr | mu , sigma );
  for ( i in 1:N ) {
    mu[i] = a + p[pj[i]] + b * lgh_s[i] + c * tsf_s[i] + cc * tsf_s[i] * lgh_s[i];
        log_lik[i]= normal_lpdf( lfr[i] | mu[i] , sigma);
  }
}
