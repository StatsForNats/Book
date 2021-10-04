data{
  int<lower=1> N;
  real y[N];
  real x[N];
  real x2[N];
}
parameters{
  real a;
  real b1;
  real b2;
  real<lower=0,upper=10> sigma;
}
model{
  vector[N] mu;
  // sigma ~ uniform( 0 , 10 );
  b2 ~ normal( -0.16 , 0.13 );
  b1 ~ normal( 3.37 , 1.85 );
  a ~ normal( 29.32 , 11.59 );
  for ( i in 1:N ) {
    mu[i] = a + b1 * x[i] + b2 * x2[i];
  }
  y ~ normal( mu , sigma );
}
generated quantities{
  vector[N] mu;
  real dev;
  dev = 0;
  for ( i in 1:N ) {
    mu[i] = a + b1 * x[i] + b2 * x2[i];
  }
  dev = dev + (-2)*normal_lpdf( y | mu , sigma );
}
