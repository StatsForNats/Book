data{
  int<lower=1> N;
  real rep_structures[N];
  real height[N];
}
parameters{
  real a;
  real b;
  real<lower=0,upper=400> sigma;
}
model{
  vector[N] mu;
  sigma ~ uniform( 0 , 400 );
  b ~ normal( 0 , 100 );
  a ~ normal( 0 , 500 );
  for ( i in 1:N ) {
    mu[i] = a + b * height[i];
  }
  rep_structures ~ normal( mu , sigma );
}
generated quantities{
  vector[N] mu;
  real dev;
  dev = 0;
  for ( i in 1:N ) {
    mu[i] = a + b * height[i];
  }
  dev = dev + (-2)*normal_lpdf( rep_structures | mu , sigma );
}
