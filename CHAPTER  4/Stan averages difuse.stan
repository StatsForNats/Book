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
  sigma ~ cauchy( 0 , 1 );
  a ~ normal( 0 , 100 );
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
