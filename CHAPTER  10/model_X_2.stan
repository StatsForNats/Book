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
  real z1;
  real z2;
  real z3;
  real q;
  vector[N_plant] p1;
  real<lower=0> sigmap1;
  vector[N_plant] p2;
  real<lower=0> sigmap2;
}
model{
  vector[N] lambda;
  vector[N] pbar;
  sigmap2 ~ cauchy( 0 , 10 );
  p2 ~ normal( 0 , sigmap2 );
  sigmap1 ~ cauchy( 0 , 10 );
  p1 ~ normal( 0 , sigmap1 );
  b3 ~ normal( 0 , 10 );
  b2 ~ normal( 0 , 10 );
  b1 ~ normal( 0 , 10 );
  q ~ normal( 0 , 100 );
  z3 ~ normal( 0 , 10 );
  z2 ~ normal( 0 , 10 );
  z1 ~ normal( 0 , 10 );
  a ~ normal( 0 , 100 );
  for ( i in 1:N ) {
    lambda[i] = q + z1 * self[i] + z2 * cross[i] + z3 * control[i] + p2[plant[i]];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:N ) {
    pbar[i] = a + b1 * self[i] + b2 * cross[i] + b3 * control[i] + p1[plant[i]];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:N )
    if (seed[i] == 0)
      target += (log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                             bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(seed[i]|lambda[i])));
      else
        target += (bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(seed[i]|lambda[i]));
}
generated quantities{
  vector[N] lambda;
  vector[N] pbar;
  real dev;
  dev = 0;
  for ( i in 1:N ) {
    lambda[i] = q + z1 * self[i] + z2 * cross[i] + z3 * control[i] + p2[plant[i]];
    lambda[i] = exp(lambda[i]);
  }
  for ( i in 1:N ) {
    pbar[i] = a + b1 * self[i] + b2 * cross[i] + b3 * control[i] + p1[plant[i]];
    pbar[i] = inv_logit(pbar[i]);
  }
  for ( i in 1:N )
    if (seed[i] == 0)
      dev = dev + (-2)*(log_sum_exp(bernoulli_lpmf(1|pbar[i]),
                                    bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(seed[i]|lambda[i])));
    else
      dev = dev + (-2)*(bernoulli_lpmf(0|pbar[i]) + poisson_lpmf(seed[i]|lambda[i]));
}
