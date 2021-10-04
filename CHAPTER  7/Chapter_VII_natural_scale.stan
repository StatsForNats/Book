data{
  int<lower=1> N;
  int surv_fin[N];
  int rep[N];
}
parameters{
  real a;
  real b;
}
model{
  vector[N] p;
  b ~ normal( 0 , 10 );
  a ~ normal( 0 , 10 );
  for ( i in 1:N ) {
    p[i] = a + b * rep[i];
  }
  surv_fin ~ binomial_logit( 1 , p );
}
generated quantities{
  vector[N] p;
  vector[N] log_lik;
  real dev;
  dev = 0;
  for ( i in 1:N ) {
    p[i] = a + b * rep[i];
  }
  dev = dev + (-2)*binomial_logit_lpmf( surv_fin | 1 , p );
  for ( i in 1:N ) {
    p[i] = a + b * rep[i];
    log_lik[i]= binomial_logit_lpmf( surv_fin[i] | 1 , p[i]);
  }
}
