data {
   int dna_alt; // a;// alt DNA
   int dna_ref; // b; // ref DNA
   int rna_alt; // k; // alt RNA
   int rna_ref; // m; // ref RNA
   real mu;     // mean of prior on r_ref
   real sigma2; // variance of prior on r_ref
   real alpha;  // shape parm of gamma prior on lambda
   real beta;   // rate parm of gamma prior on lambda
}
parameters {
   real<lower=0.001,upper=1000> r_ref;
   real<lower=0.001,upper=1000> theta;
   real<lower=0.001,upper=1000> lambda_ref; // mean of Poisson for DNA count
   real<lower=0.001,upper=1000> lambda_alt; // mean of Poisson for DNA count
}
transformed parameters {
   real<lower=0.001,upper=1000> r_alt=theta*r_ref;
}
model {
   theta ~ lognormal(0,1);
   r_ref ~ lognormal(mu,sigma2);
   lambda_alt ~ gamma(alpha,beta);
   lambda_ref ~ gamma(alpha,beta);
   rna_alt ~ poisson(r_alt*lambda_alt);
   dna_alt ~ poisson(lambda_alt);
   rna_ref ~ poisson(r_ref*lambda_ref);
   dna_ref ~ poisson(lambda_ref);
}




