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
   int<lower=0> N_RNA;      // number of RNA replicates
   int<lower=0> k[N_RNA];   // RNA alt read counts
   int<lower=0> m[N_RNA];   // RNA ref read counts
}

parameters {
   real<lower=0,upper=1> p; // alt allele freq in DNA library
   real<lower=0,upper=1> q; // alt allele freq in RNA library
   real<lower=0,upper=1> qi[N_RNA]; // alt allele freqs in RNA reps
   real<lower=2> c; // concentration parameter of beta prior for qi
}

model {
   // Parameters:
   c ~ gamma(1.1, 0.0005); // concentration parameter for prior on qi
   for(i in 1:N_RNA)
      qi[i] ~ betaModeConc(q,c);

   // Likelihoods:
   for(i in 1:N_DNA)
      a[i] ~ binomial(a[i]+b[i],p);
   for(i in 1:N_RNA)
      k[i] ~ binomial(k[i]+m[i],qi[i]);
}

generated quantities {
   real<lower=0> theta; // effect size (odds ratio)
   theta=q/(1-q)/(p/(1-p));
}


