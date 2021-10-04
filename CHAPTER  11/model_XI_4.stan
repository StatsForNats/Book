data{
  int<lower=1> N;
  int surv_fin[N];
  real stems[N];
  real lgh[N];
  int TSF[N];
  real lgh2[N];
}
parameters{
  real a;
  real b;
  real bb;
  real c;
  real d;
}
model{
  vector[N] p;
  d ~ normal( 0 , 1 );
  c ~ normal( 0 , 1 );
  bb ~ normal( 0 , 1 );
  b ~ normal( 0 , 1 );
  a ~ normal( 0 , 100 );
  for ( i in 1:N ) {
    p[i] = a + b * lgh[i] + bb * lgh2[i] + c * TSF[i] + d * stems[i];
  }
  surv_fin ~ binomial_logit( 1 , p );
}
generated quantities{
  vector[N] p;
  vector[N] log_lik;
  real dev;
  dev = 0;
  for ( i in 1:N ) {
    p[i] = a + b * lgh[i] + bb * lgh2[i] + c * TSF[i] + d * stems[i];
  }
  dev = dev + (-2)*binomial_logit_lpmf( surv_fin | 1 , p );
  for ( i in 1:N ) {
    p[i] = a + b * lgh[i] + bb * lgh2[i] + c * TSF[i] + d * stems[i];
    log_lik[i]= binomial_logit_lpmf( surv_fin[i] | 1 , p[i]);
  }
}
