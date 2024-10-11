functions {
   real PPNB_lpmf(int y,int x,real alpha,real beta,real ratio) {
      // This is a custom negative-binomial implementation that makes use
      // of a prior with parameters alpha and beta.  Ratio is RNA/DNA.
      real numer=lgamma(y+x+alpha)+y*log(ratio)+(x+alpha)*log(beta+1);
      real denom=lgamma(y+1)+lgamma(x+alpha)+(y+x+alpha)*log(beta+ratio+1);
      return numer-denom;
   }}

data {
   int<lower=0> N_VARIANTS;    // number of pools
   int<lower=0> b[N_VARIANTS]; // DNA ref read counts
   int<lower=0> m[N_VARIANTS]; // RNA ref read counts
}
parameters {
   real<lower=0> mu;        // mean of lognormal prior on r_ref
   real<lower=0> sigma2;    // variance of lognormal prior on r_ref
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
   sigma2 ~ lognormal(3,1);
   for(i in 1:N_VARIANTS) {
      r_ref[i] ~ lognormal(mu,sigma2);
      lambda[i] ~ gamma(alpha,beta);
      b[i] ~ poisson(lambda);
      m[i] ~ PPNB(b[i],alpha,beta,r_ref[i]);
   }
}




