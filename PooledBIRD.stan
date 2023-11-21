functions {
   // This function can be used via: p~betaModeConc(mode,concentration);
   real betaModeConc_lpdf(real parm,real m,real c) {
      return beta_lpdf(parm|m*(c-2)+1, (1-m)*(c-2)+1);
   }}

data {
   int<lower=0> N_POOLS;          // number of pools
   int<lower=0> MAX_DNA;          // maximum # of DNA reps
   int<lower=0> MAX_RNA;          // maximum # of RNA reps
   int<lower=0> N_DNA[N_POOLS];   // # DNA replicates
   int<lower=0> N_RNA[N_POOLS];   // # RNA replicates
   real<lower=0,upper=1> pop_freq[N_POOLS]; // Population alt allele freqs
   real<lower=0> pop_conc;        // Concentration of beta prior on p
   int<lower=0> a[N_POOLS,MAX_DNA]; // DNA alt read counts
   int<lower=0> b[N_POOLS,MAX_DNA]; // DNA ref read counts
   int<lower=0> k[N_POOLS,MAX_RNA]; // RNA alt read counts
   int<lower=0> m[N_POOLS,MAX_RNA]; // RNA ref read counts
}

parameters {
   real<lower=0> theta; // effect size (odds ratio)
   real<lower=2> c; // concentration parameter of beta prior for qi
   real<lower=0> s; // variance parameter of lognormal prior for theta
   real<lower=0,upper=1> p[N_POOLS]; // alt allele freq in DNA library
   real<lower=0,upper=1> qi[N_POOLS,MAX_RNA]; // alt allele freqs in RNA reps
}

transformed parameters {
   real<lower=0,upper=1> q[N_POOLS]; // alt allele freq in RNA
   for(j in 1:N_POOLS) q[j]=theta*p[j]/(1.0-p[j]+theta*p[j]);
}

model {
   // Parameters:
   c ~ gamma(1.1, 0.0005); // concentration parameter for prior on qi
   s ~ gamma(1.1,3);       // variance parameter for prior on theta
   log(theta)/s ~ normal(0,1); // prior on theta
   target+=-log(theta)-log(s); // Jacobian for lognormal theta prior
   for(j in 1:N_POOLS) {
      p[j] ~ betaModeConc(pop_freq[j],pop_conc);
      for(i in 1:N_RNA[j])
         qi[j,i] ~ betaModeConc(q[j],c);
      for(i in 1:N_DNA[j])
         a[j,i] ~ binomial(a[j,i]+b[j,i],p[j]);
      for(i in 1:N_RNA[j])
         k[j,i] ~ binomial(k[j,i]+m[j,i],qi[j,i]);
   }
}



