data{
  int<lower=1> N;
  real height[N];
}
parameters{
  real a;
  real<lower=0> sigma;
}
model{
  vector[N] mu;
  sigma ~ lognormal( 2.245235  , 0.06751135 );
  a ~ normal( 34.39132 , 3.927421 );
  for ( i in 1:N ) {
    mu[i] = a;
  }
  height ~ normal( mu , sigma );
}
generated quantities{
  vector[N] mu;
  real dev;
  dev = 0;
  for ( i in 1:N ) {
    mu[i] = a;
  }
  dev = dev + (-2)*normal_lpdf( height | mu , sigma );
}
