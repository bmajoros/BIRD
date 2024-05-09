functions {
   // This function can be used via: p~betaModeConc(mode,concentration);
   real betaModeConc_lpdf(real parm,real m,real c) {
      return beta_lpdf(parm|m*(c-2)+1, (1-m)*(c-2)+1);
   }}

// TYPES OF POOLS: 1=HET, 2=REF, 3=ALT

data {
   int<lower=0> N_POOLS;          // number of pools
   int<lower=1,upper=3> POOL_TYPE[N_POOLS]; // 1=het, 2=ref, 3=alt
   int<lower=0> a[N_POOLS]; // DNA alt read counts
   int<lower=0> b[N_POOLS]; // DNA ref read counts
   int<lower=0> k[N_POOLS]; // RNA alt read counts
   int<lower=0> m[N_POOLS]; // RNA ref read counts
   real<lower=0> mu; // mean of lognormal prior on r_ref and r_alt
   real<lower=0> sigma2; // variance of lognormal prior on r_ref and r_alt
}

parameters {
   real<lower=0.001,upper=1000> d; // Dispersion parameter for negative binomials
   real<lower=0.001,upper=1000> r_ref; // ratio RNA/DNA for reference allele
   real<lower=0.001,upper=1000> r_alt; // ratio RNA/DNA for alternate allele
}

transformed parameters {
}

model {
   // Parameters:
   r_ref ~ lognormal(mu,sigma2);
   r_ref ~ lognormal(mu,sigma2);
   for(i in 1:N_POOLS) {
      if(POOL_TYPE[i]==1) { // HETEROZYGOUS POOL
	 m[i] ~ neg_binomial_2(b[i]*r_ref,1.0/d);
	 k[i] ~ neg_binomial_2(a[i]*r_alt,1.0/d);
      }
      else if(POOL_TYPE[i]==2) { // REF POOL
	 m[i] ~ neg_binomial_2(b[i]*r_ref,1.0/d);
      }
      else if(POOL_TYPE[i]==3) { // ALT POOL
	 k[i] ~ neg_binomial_2(a[i]*r_alt,1.0/d);
      }
   }
}



