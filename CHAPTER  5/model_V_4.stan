data{
  int<lower=1> N;
  real rep_structures[N];
  real height[N];
}
parameters{
  real a;
  real b;
  real<lower=0,upper=100> sigma;
}
model{
  vector[N] mu;
  // sigma ~ uniform( 0 , 100 );
  b ~ normal( 0 , 1 );
  a ~ normal( 0 , 1 );
  for ( i in 1:N ) {
    mu[i] = a * height[i]^b;
  }
  rep_structures ~ lognormal( mu , sigma );
}
generated quantities{
  vector[N] mu;
  real dev;
  dev = 0;
  for ( i in 1:N ) {
    mu[i] = a * height[i]^b;
  }
  dev = dev + (-2)*lognormal_lpdf( rep_structures | mu , sigma );
}
