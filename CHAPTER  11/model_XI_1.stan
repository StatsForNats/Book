data{
  int<lower=1> N;
  int surv_fin[N];
  real stems[N];
  real lgh[N];
  int TSF[N];
}
parameters{
  real a;
  real b;
  real bc;
  real bd;
  real c;
  real d;
  real cd;
}
model{
  vector[N] p;
  cd ~ normal( 0 , 1);
  bd ~ normal( 0 , 1 );
  bc ~ normal( 0 , 1 );
  d ~ normal( 0 , 1 );
  c ~ normal( 0 , 1 );
  b ~ normal( 0 , 1 );
  a ~ normal( 0 , 100 );
  for ( i in 1:N ) {
    p[i] = a + b * lgh[i] + c * TSF[i] + bc * lgh[i] * TSF[i] + d * stems[i] + bd * lgh[i] *  stems[i] + cd * TSF[i] *  stems[i];
  }
  surv_fin ~ binomial_logit( 1 , p );
}
generated quantities{
  vector[N] p;
  vector[N] log_lik;
  real dev;
  dev = 0;
  for ( i in 1:N ) {
    p[i] = a + b * lgh[i] + c * TSF[i] + bc * lgh[i] * TSF[i] + d * stems[i] + bd * lgh[i] *  stems[i] + cd * TSF[i] *  stems[i];
  }
  dev = dev + (-2)*binomial_logit_lpmf( surv_fin | 1 , p );
  for ( i in 1:N ) {
    p[i] = a + b * lgh[i] + c * TSF[i] + bc * lgh[i] * TSF[i] + d * stems[i] + bd * lgh[i] *      stems[i] + cd * TSF[i] *  stems[i];
    log_lik[i]= binomial_logit_lpmf( surv_fin[i] | 1 , p[i]);
  }
}
