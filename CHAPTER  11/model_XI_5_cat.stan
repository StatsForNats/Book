data{
  int<lower=1> N;
  int surv_fin[N];
  real lgh[N];
  int TSF2[N];
  int TSF3[N];
  real lgh2[N];
}
parameters{
  real a;
  real b;
  real bb;
  real c2;
  real c3;
}
model{
  vector[N] p;
  c2 ~ normal( 0 , 1 );
  c3 ~ normal( 0 , 1 );
  b ~ normal( 0 , 1 );
  a ~ normal( 0 , 100 );
  for ( i in 1:N ) {
    p[i] = a + b * lgh[i] + bb * lgh2[i] + c2 * TSF2[i] + c3 * TSF3[i];
  }
  surv_fin ~ binomial_logit( 1 , p );
}
generated quantities{
  vector[N] p;
  vector[N] log_lik;
  real dev;
  dev = 0;
  for ( i in 1:N ) {
    p[i] = a + b * lgh[i] + bb * lgh2[i] + c2 * TSF2[i] + c3 * TSF3[i];
  }
  dev = dev + (-2)*binomial_logit_lpmf( surv_fin | 1 , p );
  for ( i in 1:N ) {
    p[i] = a + b * lgh[i] + bb * lgh2[i] + c2 * TSF2[i] + c3 * TSF3[i];
    
    log_lik[i]= binomial_logit_lpmf( surv_fin[i] | 1 , p[i]);
  }
}
