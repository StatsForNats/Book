data{
  int<lower=1> N;
  real lnrep[N];
  real lnht[N];
}
parameters{
  real a;
  real b;
  real<lower=0,upper=10> sigma;
}
model{
  vector[N] mu;
  sigma ~ uniform( 0 , 10 );
  b ~ normal( 0 , 1 );
  a ~ normal( 0 , 1 );
  for ( i in 1:N ) {
    mu[i] = a + b * lnht[i];
  }
  lnrep ~ normal( mu , sigma );
}
generated quantities{
  vector[N] mu;
  real dev;
  dev = 0;
  for ( i in 1:N ) {
    mu[i] = a + b * lnht[i];
  }
  dev = dev + (-2)*normal_lpdf( lnrep | mu , sigma );
}
