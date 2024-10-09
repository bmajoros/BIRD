data {
   int dna_alt; // a;// alt DNA
   int dna_ref; // b; // ref DNA
   int rna_alt; // k; // alt RNA
   int rna_ref; // m; // ref RNA
   real beta1;
   real beta2;
}
parameters {
   real<lower=0.001,upper=0.999> p;
   real<lower=0.001,upper=1000> theta;
}
transformed parameters {
   real<lower=0,upper=1> q=theta*p/(1.0-p+theta*p);
}
model {
   theta ~ lognormal(0,1);
   p ~ beta(beta1,beta2);
   dna_alt ~ binomial(dna_alt+dna_ref,p);
   rna_alt ~ binomial(rna_alt+rna_ref,q);
}




