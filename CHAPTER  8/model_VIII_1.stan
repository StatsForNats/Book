data{
  int<lower=1> N;
  real lfr[N];
  real lgh[N];
}
parameters{
  real a;
  real b;
  real<lower=0,upper=10> sigma;
}
model{
  vector[N] mu;
  // sigma ~ uniform( 0 , 10 );
  b ~ normal( 0 , 10 );
  a ~ normal( 0 , 100 );
  for ( i in 1:N ) {
    mu[i] = a + b * lgh[i];
  }
  lfr ~ normal( mu , sigma );
}
generated quantities{
  vector[N] mu;
  real dev;
  vector[N] log_lik;
  dev = 0;
  for ( i in 1:N ) {
    mu[i] = a + b * lgh[i];
  }
  dev = dev + (-2)*normal_lpdf( lfr | mu , sigma );
  for ( i in 1:N ) {
    mu[i] = a + b * lgh[i];
    log_lik[i]= normal_lpdf( lfr[i] | mu[i] , sigma);
  }
}
