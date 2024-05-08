functions {
   // This function can be used via: p~betaModeConc(mode,concentration);
   real betaModeConc_lpdf(real parm,real m,real c) {
      return beta_lpdf(parm|m*(c-2)+1, (1-m)*(c-2)+1);
   }}

// TYPES OF POOLS: 1=HET, 2=REF, 3=ALT

data {
   int<lower=0> N_POOLS;          // number of pools
   int<lower=1,upper=3> POOL_TYPE[N_POOLS]; // 1=het, 2=ref, 3=alt
   real<lower=0,upper=1> pop_freq[N_POOLS]; // Population alt allele freqs
   real<lower=2> pop_conc;        // Concentration of beta prior on p
   real<lower=0> d; // Dispersion parameter for negative binomials
   int<lower=0> a[N_POOLS]; // DNA alt read counts
   int<lower=0> b[N_POOLS]; // DNA ref read counts
   int<lower=0> k[N_POOLS]; // RNA alt read counts
   int<lower=0> m[N_POOLS]; // RNA ref read counts
   real<lower=0> mu; // mean of lognormal prior on r_ref
   real<lower=0> sigma2; // variance of lognormal prior on r_ref
}

parameters {
   real<lower=0> theta; // effect size (odds ratio)
   real<lower=0> s; // variance parameter of lognormal prior for theta
   real<lower=0,upper=1> p[N_POOLS]; // alt allele freq in DNA library
   real<lower=0> r_ref; // ratio RNA/DNA for reference allele
}

transformed parameters {
   real<lower=0,upper=1> q[N_POOLS]; // alt allele freq in RNA
   real<lower=0> r_alt=theta*r_ref; // ratio RNA/DNA for alternate allele
   for(i in 1:N_POOLS) {
      if(POOL_TYPE[i]==1) q[i]=theta*p[i]/(1.0-p[i]+theta*p[i]);
   }
}

model {
   // Parameters:
   s ~ gamma(1.1,3);       // variance parameter for prior on theta
   //log(theta)/s ~ normal(0,1); // prior on theta
   //target+=-log(theta)-log(s); // Jacobian for lognormal theta prior
   theta ~ lognormal(0,s);
   r_ref ~ lognormal(mu,sigma2);
   for(i in 1:N_POOLS) {
      if(POOL_TYPE[i]==1) { // HETEROZYGOUS POOL
         p[i] ~ betaModeConc(pop_freq[i],pop_conc);
         a[i] ~ binomial(a[i]+b[i],p[i]);
         k[i] ~ binomial(k[i]+m[i],q[i]);
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



