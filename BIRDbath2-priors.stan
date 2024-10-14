functions {
   real PPNB_lpmf(int y,int x,real alpha,real beta,real ratio) {
      // This is a custom negative-binomial implementation that makes use
      // of a prior with parameters alpha and beta.  Ratio is RNA/DNA.
      real numer=lgamma(y+x+alpha)+y*log(ratio)+(x+alpha)*log(beta+1);
      real denom=lgamma(y+1)+lgamma(x+alpha)+(y+x+alpha)*log(beta+ratio+1);
      return numer-denom;
   }}

data {
   int<lower=0> N_VARIANTS;    // number of variants
   int<lower=0> N_POOLS;       // number of pools per variant
   int<lower=0> b[N_VARIANTS,N_POOLS]; // DNA ref read counts
   int<lower=0> m[N_VARIANTS,N_POOLS]; // RNA ref read counts
}
parameters {
   real<lower=0> mu;        // mean of lognormal prior on r_ref
   real<lower=0> sigma;    // variance of lognormal prior on r_ref
   real<lower=0> alpha;     // shape parameter of gamma prior in NB
   real<lower=0> beta;      // rate parameter of gamma prior in NB
   real<lower=0.000001> r_ref[N_VARIANTS]; // ratio RNA/DNA for ref allele
   real<lower=0.000001> lambda[N_VARIANTS]; // poison rate parameter
}
model {
   // Parameters:
   alpha ~ lognormal(3,1);
   beta ~ lognormal(3,1);
   mu ~ lognormal(3,1);
   sigma ~ lognormal(3,1);
   for(i in 1:N_VARIANTS) {
      r_ref[i] ~ lognormal(mu,sigma);
      lambda[i] ~ gamma(alpha,beta);
      for(j in 1:N_POOLS) {
         b[i][j] ~ poisson(lambda);
         m[i][j] ~ PPNB(b[i][j],alpha,beta,r_ref[i]);
      }
   }
}




