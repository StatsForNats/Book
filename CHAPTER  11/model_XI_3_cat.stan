data{
  int<lower=1> N;
  int surv_fin[N];
  real stems[N];
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
  real d;
  real c2d;
  real c3d;
  real bc2;
  real bc3;
  real bbc2;
  real bbc3;
  real bd;
  real bbd;
}
model{
  vector[N] p;
  bbd ~ normal( 0 , 1 );
  bd ~ normal( 0 , 1 );
  bbc2 ~ normal( 0 , 1 ); 
  bbc3 ~ normal( 0 , 1 );
  bc2 ~ normal( 0 , 1 );
  bc3 ~ normal( 0 , 1 );
  bb ~ normal( 0 , 1 );
  c2d ~ normal(0 , 1 );
  c3d ~ normal(0 , 1 );
  d ~ normal( 0 , 1 );
  c2 ~ normal( 0 , 1 );
  c3 ~ normal( 0 , 1 );
  b ~ normal( 0 , 1 );
  a ~ normal( 0 , 100 );
  for ( i in 1:N ) {
    p[i] = a + b * lgh[i] + bb * lgh2[i] + c2 * TSF2[i] + c3 * TSF3[i] + 
      bc2 * lgh[i] * TSF2[i] + bc3 * lgh[i] * TSF3[i] +  bbc2 * lgh2[i] * TSF2[i] +
      bbc3 * lgh2[i] * TSF3[i] +d * stems[i] + bd * lgh[i] * stems[i] +
      bbd * lgh2[i] * stems[i] + c2d*TSF2[i]* stems[i]+ c3d*TSF3[i]* stems[i];
  }
  surv_fin ~ binomial_logit( 1 , p );
}
generated quantities{
  vector[N] p;
  vector[N] log_lik;
  real dev;
  dev = 0;
  for ( i in 1:N ) {
    p[i] = a + b * lgh[i] + bb * lgh2[i] + c2 * TSF2[i] + c3 * TSF3[i] + 
      bc2 * lgh[i] * TSF2[i] + bc3 * lgh[i] * TSF3[i] +  bbc2 * lgh2[i] *TSF2[i] +
      bbc3 * lgh2[i] *TSF3[i] +d * stems[i] + bd * lgh[i] * stems[i] +
      bbd * lgh2[i] * stems[i] + c2d*TSF2[i]* stems[i] + c3d*TSF3[i]* stems[i];
  }
  dev = dev + (-2)*binomial_logit_lpmf( surv_fin | 1 , p );
  for ( i in 1:N ) {
    p[i] = a + b * lgh[i] + bb * lgh2[i] + c2 * TSF2[i] + c3 * TSF3[i] + 
      bc2 * lgh[i] * TSF2[i] + bc3 * lgh[i] * TSF3[i] +  bbc2 * lgh2[i] *TSF2[i] +
      bbc3 * lgh2[i] *TSF3[i] +d * stems[i] + bd * lgh[i] * stems[i] +
      bbd * lgh2[i] * stems[i] + c2d*TSF2[i]* stems[i]+ c3d*TSF3[i]* stems[i];
    
    log_lik[i]= binomial_logit_lpmf( surv_fin[i] | 1 , p[i]);
  }
}
