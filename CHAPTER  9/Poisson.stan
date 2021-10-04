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
  vector[N_population] p;
  vector[N_yrs] y;
  real<lower=0> sigmap;
  real<lower=0> sigmay;
}
model{
  vector[N] lambda;
  sigmay ~ cauchy( 0 , 1 );
  sigmap ~ cauchy( 0 , 1 );
  y ~ normal( 0 , sigmay );
  p ~ normal( 0 , sigmap );
  b0 ~ normal( 0 , 10 );
  for ( i in 1:N ) {
    lambda[i] = b0 + p[population[i]] + y[yrs[i]];
  }
  seeds ~ poisson_log( lambda );
}
generated quantities{
  vector[N] lambda;
  vector[N] pseeds;
  real dev;
  real xsqp;
  dev = 0;
  xsqp = 0;
  for ( i in 1:N ) {
    lambda[i] = b0 + p[population[i]] + y[yrs[i]];
    pseeds[i] = exp(lambda[i]);
  }
  dev = dev + (-2)*poisson_log_lpmf( seeds | lambda );
  for ( i in 1:N ) {
    xsqp = xsqp + pow(seeds[i]-mean(pseeds), 2)/mean(pseeds);
  }
}
