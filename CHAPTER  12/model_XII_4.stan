data{
  int<lower=1> N;
  int<lower=1> N_pj;
  real lfr[N];
  real lgh_s[N];
  real tsf_s[N];
  int pj[N];
  matrix[N_pj,N_pj] Dmat;
}
parameters{
  vector[N_pj] p;
  real a;
  real b;
  real c;
  real<lower=0> etasq;
  real<lower=0> rhosq;
  real<lower=0> sigmap;
  real<lower=0> sigma;
}
model{
  matrix[N_pj,N_pj] SIGMA_Dmat;
  vector[N] mu;
  sigma ~ cauchy( 0 , 1 );
  sigmap ~ exponential( 1 );
  rhosq ~ cauchy( 0 , 1 );
  etasq ~ cauchy( 0 , 1 );
  c ~ normal( 0 , 1 );
  b ~ normal( 0 , 1 );
  a ~ normal( 0 , 50 );
  for ( i in 1:(N_pj-1) )
    for ( j in (i+1):N_pj ) {
      SIGMA_Dmat[i,j] = etasq*exp(-rhosq*pow(Dmat[i,j],2));
      SIGMA_Dmat[j,i] = SIGMA_Dmat[i,j];
    }
  for ( k in 1:N_pj )
    SIGMA_Dmat[k,k] = etasq + sigmap;
  p ~ multi_normal( rep_vector(0,N_pj) , SIGMA_Dmat );
  for ( i in 1:N ) {
    mu[i] = a + p[pj[i]] + b * lgh_s[i] + c * tsf_s[i];
  }
  lfr ~ normal( mu , sigma );
}
generated quantities{
  matrix[N_pj,N_pj] SIGMA_Dmat;
  vector[N] mu;
  vector[N] log_lik;
  real dev;
  dev = 0;
  for ( i in 1:N ) {
    mu[i] = a + p[pj[i]] + b * lgh_s[i] + c * tsf_s[i];
  }
  dev = dev + (-2)*normal_lpdf( lfr | mu , sigma );
   for ( i in 1:N ) {
    mu[i] = a + p[pj[i]] + b * lgh_s[i] + c * tsf_s[i];
        log_lik[i]= normal_lpdf( lfr[i] | mu[i] , sigma);
  }
}
