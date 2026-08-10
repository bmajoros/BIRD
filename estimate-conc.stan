functions {
   // This function can be used via: p~betaModeConc(mode,concentration);
   real betaModeConc_lpdf(real parm,real m,real c) {
      return beta_lpdf(parm|m*(c-2)+1, (1-m)*(c-2)+1);
   }}

data {
   int<lower=0> N_POOLS;          // number of pools
   int<lower=0> MAX_DNA;          // maximum # of DNA reps
   int<lower=0> N_DNA[N_POOLS];   // # DNA replicates
   real<lower=0,upper=1> pop_freq[N_POOLS]; // Population alt allele freqs
   int<lower=0> a[N_POOLS,MAX_DNA]; // DNA alt read counts
   int<lower=0> b[N_POOLS,MAX_DNA]; // DNA ref read counts
   //int D[N_POOLS];
   //for(j in 1:N_POOLS) {
   //   D[j]=to_int(pop_freq[j]*5);
   //   }
}

parameters {
   real<lower=2.001> pop_conc;        // Concentration of beta prior on p
   real<lower=0,upper=1> p[N_POOLS]; // alt allele freq in DNA library
}

model {
   // Parameters:
   //pop_conc ~ gamma(1.01, 0.0005); // concentration parameter for prior on p
   pop_conc ~ gamma(1, 0.25);
   //2~binomial(10,p[1]);
   //1~binomial(10,p[2]);
   //1~binomial(10,p[3]);
   for(j in 1:N_POOLS) {
      p[j] ~ betaModeConc(pop_freq[j],pop_conc);
      //D[j]~binomial(10,p[j]);
      for(i in 1:N_DNA[j]) {
      	 a[j,i] ~ binomial(a[j,i]+b[j,i],p[j]);
         //100*a[j,i] ~ binomial(100*a[j,i]+100*b[j,i],p[j]);
	 //print("counts: ",a[j,i],"   ",b[j,i]);
	 //print("pj: ",p[j]);
	 //print("conc: ",pop_conc);
	 }
   }
}



