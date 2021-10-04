data{
  int<lower=1> N;
  int<lower=1> N_plant;
  int seed[N];
  int plant[N];
  real self[N];
  real cross[N];
  real control[N];
}
parameters{
  real a;
  real b1;
  real b2;
  real b3;
  vector[N_plant] p;
  real<lower=0> sigmap;
}
model{
  vector[N] lambda;
  sigmap ~ cauchy( 0 , 10 );
  p ~ normal( 0 , sigmap );
  b3 ~ normal( 0 , 10 );
  b2 ~ normal( 0 , 10 );
  b1 ~ normal( 0 , 10 );
  a ~ normal( 0 , 100 );
  for ( i in 1:N ) {
    lambda[i] = a + b1 * self[i] + b2 * cross[i] + b3 * control[i] + p[plant[i]];
  }
  seed ~ poisson_log( lambda );
}
generated quantities{
  vector[N] lambda;
  vector[N] log_lik;
  real dev;
  dev = 0;
  for ( i in 1:N ) {
    lambda[i] = a + b1 * self[i] + b2 * cross[i] + b3 * control[i] + p[plant[i]];
  }
  dev = dev + (-2)*poisson_log_lpmf( seed | lambda );
   for ( i in 1:N ) {
    lambda[i] = a + b1 * self[i] + b2 * cross[i] + b3 * control[i] + p[plant[i]];
    log_lik[i]= poisson_log_lpmf( seed[i] | lambda ); 
  }
}
