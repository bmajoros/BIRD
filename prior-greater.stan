functions {
   // This function can be used via: g~gammaModeSD(mode,sd);
   real gammaModeSD_lpdf(real parm,real m,real sd) {
      real r=(m+sqrt(m^2+4*sd^2))/(2*sd^2);
      real s=1+m*r;
      return gamma_lpdf(parm|r,s);
   }

   // This function can be used via: p~betaModeConc(mode,concentration);
   real betaModeConc_lpdf(real parm,real m,real c) {
      return beta_lpdf(parm|m*(c-2)+1, (1-m)*(c-2)+1);
   }  
}

data {
   real<lower=0,upper=1> v; // alt allele freq in VCF file
   int<lower=0> N_DNA;      // number of DNA replicates
   int<lower=0> a[N_DNA];   // DNA alt read counts
   int<lower=0> b[N_DNA];   // DNA ref read counts
   int<lower=0> N_RNA;      // number of RNA replicates
   int<lower=0> k[N_RNA];   // RNA alt read counts
   int<lower=0> m[N_RNA];   // RNA ref read counts
}

parameters {
   real<lower=0,upper=1> p; // alt allele freq in DNA library
   real<lower=0,upper=1> qi[N_RNA]; // alt allele freqs in RNA replicates
   real<lower=1> theta;
   real<lower=2> c1; // concentration parameter of beta prior for qi
   real<lower=2> c2; // concentration parameter of beta prior for p
   real<lower=0> c3; // std.dev. parameter of lognormal prior for theta
}

transformed parameters { // ORDER MATTERS!
   real<lower=0,upper=1> q; // alt allele freq in RNA (theoretical)
   q=theta*p/(1.0-p+theta*p);
}

model {
   // Parameters
   c1 ~ gamma(1.1, 0.005);
   c2 ~ gamma(1.1, 0.005);
   1/c3^2 ~ gammaModeSD(1,1); // c3^2 is precision
   target+=log(2)-3*log(c3); // Jacobian for gamma prior on precision
   log(theta)/c3 ~ normal(0,1);
   target+=-log(theta)-log(c3); // Jacobian for lognormal(0,1/c3)
   p ~ betaModeConc(v,c2);
   for(i in 1:N_RNA)
      qi[i] ~ betaModeConc(q,c1);
}

generated quantities {
   real LL=0;
   for(i in 1:N_DNA)
      LL+=binomial_lpmf(a[i]|a[i]+b[i],p);
   for(i in 1:N_RNA)
      LL+=binomial_lpmf(k[i]|k[i]+m[i],qi[i]);
}

