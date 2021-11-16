data{/////////////////////////////////////////////////////////////////////
    int<lower=1> S;    // number of sites or samples (nrow)
    int<lower=1> K[S];   // number of replicates per site (ncol)
    int<lower=0> N[S]; // detected or not at that site/replicate
    int z[S];   
}
parameters{/////////////////////////////////////////////////////////////////////
    real<lower=0,upper=1> psi;  //commonness parameter
    real<lower=0,upper=1> p11; //true positive detection rate
    real<lower=0,upper=1> p10; //false positive detection rate
}
transformed parameters{/////////////////////////////////////////////////////////////////////
}
model{/////////////////////////////////////////////////////////////////////
  real p[S];
  
    for (i in 1:S){
			z[i] ~ bernoulli(psi);
			p[i] = z[i]*p11 + (1-z[i])*p10;
			N[i] ~ binomial(K[i], p[i]);
	}; 
  
  //priors
  psi ~ beta(5,5); 
  p11 ~ beta(5,5); 
  p10 ~ beta(1,50);
}
generated quantities{
  real<lower=0,upper=1> Occupancy_prob[S];    //after inferring parameters above, now calculate occupancy probability for each observation. Equation from Lahoz-Monfort et al. 2015
  for (i in 1:S){
       Occupancy_prob[i]  = (psi*(p11^N[i])*(1-p11)^(K[i]-N[i])) 
       / ((psi*(p11^N[i])*(1-p11)^(K[i]-N[i])) 
          + (((1-psi)*(p10^N[i]))*((1-p10)^(K[i]-N[i])))
         );
      }
}

