data{
  int<lower=1> N;
  real lfr[N];
  real lgh_s[N];
  real tsf_s[N];
}
parameters{
  real a;
  real b;
  real c;
  real cc;
  real<lower=0,upper=10> sigma;
}
model{
  vector[N] mu;
  // sigma ~ uniform( 0 , 1 );
  cc ~ normal( 0 , 1 );
  c ~ normal( 0 , 1 );
  b ~ normal( 0 , 1 );
  a ~ normal( 0 , 50 );
  for ( i in 1:N ) {
    mu[i] = a + b * lgh_s[i] + c * tsf_s[i] + cc * tsf_s[i] * lgh_s[i];
  }
  lfr ~ normal( mu , sigma );
}
generated quantities{
  vector[N] mu;
  real dev;
  vector[N] log_lik;
  dev = 0;
  for ( i in 1:N ) {
    mu[i] = a + b * lgh_s[i] + c * tsf_s[i] + cc * tsf_s[i] * lgh_s[i];
  }
  dev = dev + (-2)*normal_lpdf( lfr | mu , sigma );
      for ( i in 1:N ) {
    mu[i] = a + b * lgh_s[i] + c * tsf_s[i] + cc * tsf_s[i] * lgh_s[i];
    log_lik[i]= normal_lpdf( lfr[i] | mu[i] , sigma);
  }
}
