functions {
   // This function can be used via: p~betaModeConc(mode,concentration);
   real betaModeConc_lpdf(real parm,real m,real c) {
      return beta_lpdf(parm|m*(c-2)+1, (1-m)*(c-2)+1);
   }  
}

data {
   int<lower=0> N_DNA;      // number of DNA replicates
   int<lower=0> a[N_DNA];   // DNA alt read counts
   int<lower=0> b[N_DNA];   // DNA ref read counts
   real<lower=0,upper=1> pi;// Theoretical allele freq from pool genotypes
}

parameters {
   real<lower=0,upper=1> p; // alt allele freq in DNA library
   real<lower=2> conc; // concentration parameter of beta prior on p
}

model {
   // Parameters:
   conc ~ gamma(1.1, 0.0005);
   p ~ betaModeConc(pi,conc);

   // Likelihoods:
   for(i in 1:N_DNA)
      a[i] ~ binomial(a[i]+b[i],p);
}



