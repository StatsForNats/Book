data{
  int<lower=1> N;
  int<lower=1> N_population;
  int<lower=1> N_yrs;
  int seeds[N];
  int population[N];
  int yrs[N];
}
parameters{
  real b0;
  real<lower=0> scale;
  vector[N_population] p;
  vector[N_yrs] y;
  real<lower=0> sigmap;
  real<lower=0> sigmay;
}
model{
  vector[N] pbar;
  scale ~ cauchy( 0 , 10 );
  sigmay ~ cauchy( 0 , 1 );
  sigmap ~ cauchy( 0 , 1 );
  y ~ normal( 0 , sigmay );
  p ~ normal( 0 , sigmap );
  b0 ~ normal( 0 , 100 );
  for ( i in 1:N ) {
    pbar[i] = b0 + p[population[i]] + y[yrs[i]];
    pbar[i] = exp(pbar[i]);
  }
  seeds ~ neg_binomial_2( pbar , scale );
}
generated quantities{
   vector[N] pbar;
    real dev;
    real xsq;
   dev = 0;
   xsq = 0;
   for ( i in 1:N ) {
      pbar[i] = b0 + p[population[i]] + y[yrs[i]];
      pbar[i] = exp(pbar[i]);
      xsq = xsq + pow(seeds[i]-exp(b0), 2)/(exp(b0)+(pow(exp(b0), 2)/scale));
     }
  dev = dev + (-2)*neg_binomial_2_lpmf( seeds | pbar , scale );
}
