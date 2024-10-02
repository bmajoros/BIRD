functions {
   real PPNB_lpmf(int y,int x,real alpha,real beta,real ratio) {
      // This is a custom negative-binomial implementation that makes use
      // of a prior with parameters alpha and beta.  Ratio is RNA/DNA.
      real numer=lgamma(y+x+alpha)+y*log(ratio)+(x+alpha)*log(beta+1);
      real denom=lgamma(y+1)+lgamma(x+alpha)+(y+x+alpha)*log(beta+ratio+1);
      return numer-denom;
   }}
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
}
transformed parameters {
   real<lower=0.001,upper=1000> r_alt=theta*r_ref;
}
model {
   theta ~ lognormal(0,1);
   r_ref ~ lognormal(mu,sigma2);
   target += PPNB_lpmf(rna_alt|dna_alt,alpha,beta,r_alt);
   target += PPNB_lpmf(rna_ref|dna_ref,alpha,beta,r_ref);
}




