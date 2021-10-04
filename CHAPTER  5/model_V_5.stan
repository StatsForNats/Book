data{
  int<lower=1> N;
  int rep_structures[N];
  real height[N];
}
parameters{
  real a;
  real b;
  real<lower=0> scale;
}
model{
  vector[N] pbar;
  scale ~ cauchy( 0 , 10 );
  b ~ normal( 0 , 10 );
  a ~ normal( 0 , 10 );
  for ( i in 1:N ) {
    pbar[i] = a + b * height[i];
    pbar[i] = exp(pbar[i]);
  }
  rep_structures ~ neg_binomial_2( pbar , scale );
}
generated quantities{
  vector[N] pbar;
  real dev;
  dev = 0;
  for ( i in 1:N ) {
    pbar[i] = a + b * height[i];
    pbar[i] = exp(pbar[i]);
  }
  dev = dev + (-2)*neg_binomial_2_lpmf( rep_structures | pbar , scale );
}