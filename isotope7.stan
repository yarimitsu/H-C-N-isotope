# Isotope model for marine, coastal and freshwater

data {
  int<lower=1> N;   // number of observations
  int<lower=1> J;   // number of groups
  int<lower=1> K;   // number of sources
  int M;            // number of taxa  
  int tx [J];       // taxa of each species
  int grp[N];       // species or genus of each observation
  vector[N] dC_obs; // data: 13C of consumer tissues
  vector[N] dN_obs; // data: 15N of consumer tissues
  vector[N] dH_obs; // data: 2H of consumer tissues
  vector[J] dN_g;   // data: mean 15N of each group
  vector[J] dH_w;   // covariate: 2H of water
}

parameters {
  real dN_base;                        // 15N of primary consumers
  ordered[3] dC_ord;
  //real dC1;                            // 13C marine source
  //real dC2;                            // 13C coastal source
  //real <upper = -22> dC3;              // 13C freshwater source
  real dN[K];                          // 15N marine, coastal and freshwater sources
  real dH[K];                          // 2H marine, coastal and freshwater sources
  simplex[K] phi[J];                   // phi is defined as a unit simplex and thus sum(phi)=1
  real <upper = 2> Delta_C[M];         // trophic fractionation of C
  real Delta_N[M];                     // trophic fractionation of N
  real Delta_H;                        // trophic fractionation of H
  real<lower=0, upper = 0.6> omega[J]; // proportion of 2H due to ambient water dH_w
  real<lower=0> fblki;             // blki fractionation
  real<lower=0> sigma_src[3,K];    // sigma source parameters marine, coastal, freshwater for dC, dH, dN 
  real<lower=0> sigma_frc[2,J];    // sigma C_tot, N_tot
  real<lower=0> sigma_frcH;        // sigma Delta_H
  real<lower=0> sigma_omega[J];    // sigma omega
  real<lower=0> sigma_fblki;       // sigma blki fractionation

}
transformed parameters {
  real C_tot[J];
  real N_tot[J];
  real dC_exp[J];
  real dN_exp[J];
  real dH_exp[J];
  real sigma_C[J];
  real sigma_N[J];
  real sigma_H[J];
  real <lower = 1.8, upper = 5> tau[J];
  real dC1;                            // 13C marine source
  real dC2;                            // 13C coastal source
  real dC3;              // 13C freshwater source

  dC1 <- dC_ord[2];
  dC2 <- dC_ord[3];
  dC3 <- dC_ord[1];
  
  for (j in 1:J) {
    tau[j] <- 2 + (dN_g[j] - dN_base)/Delta_N[tx[j]];
    C_tot[j]  <- Delta_C[tx[j]] * (tau[j]-1);
    N_tot[j]  <- Delta_N[tx[j]] * (tau[j]-1); 
    dC_exp[j] <- (phi[j][1] * (dC1 + C_tot[j])) + 
                 (phi[j][2] * (dC2 + C_tot[j])) + 
                 (phi[j][3] * (dC3 + C_tot[j]));
    dN_exp[j] <- (phi[j][1] * (dN[1] + N_tot[j])) + 
                 (phi[j][2] * (dN[2] + N_tot[j])) + 
                 (phi[j][3] * (dN[3] + N_tot[j]));
   }


  dH_exp[1] <-   (omega[1] * dH_w[1]) + ((1 - omega[1]) * (
                 (phi[1][1] * (dH[1] + Delta_H)) + 
                 (phi[1][2] * (dH[2] + Delta_H)) + 
                 (phi[1][3] * (dH[3] + Delta_H)))) + fblki;
  
  for (j in 2:J) {
    dH_exp[j] <- (omega[j] * dH_w[j]) + ((1 - omega[j]) * (
                 (phi[j][1] * (dH[1] + Delta_H)) + 
                 (phi[j][2] * (dH[2] + Delta_H)) + 
                 (phi[j][3] * (dH[3] + Delta_H))));
    }
  for (j in 1:J) {
    sigma_C[j] <- sqrt((pow(phi[j][1],2) * (pow(sigma_src[1,1],2) + pow(sigma_frc[1,j],2))) + 
                       (pow(phi[j][2],2) * (pow(sigma_src[1,2],2) + pow(sigma_frc[1,j],2))) + 
                       (pow(phi[j][3],2) * (pow(sigma_src[1,3],2) + pow(sigma_frc[1,j],2))));
    sigma_N[j] <- sqrt((pow(phi[j][1],2) * (pow(sigma_src[2,1],2) + pow(sigma_frc[2,j],2))) + 
                       (pow(phi[j][2],2) * (pow(sigma_src[2,2],2) + pow(sigma_frc[2,j],2))) + 
                       (pow(phi[j][3],2) * (pow(sigma_src[2,3],2) + pow(sigma_frc[2,j],2))));
     }
  
    sigma_H[1] <- sqrt((pow(phi[1][1],2) * (pow(sigma_src[3,1],2) + pow(sigma_frcH,2) + pow(sigma_omega[1],2) + pow(sigma_fblki,2))) + 
                       (pow(phi[1][2],2) * (pow(sigma_src[3,2],2) + pow(sigma_frcH,2) + pow(sigma_omega[1],2) + pow(sigma_fblki,2))) + 
                       (pow(phi[1][3],2) * (pow(sigma_src[3,3],2) + pow(sigma_frcH,2) + pow(sigma_omega[1],2) + pow(sigma_fblki,2))));

  for (j in 2:J) {
    sigma_H[j] <- sqrt((pow(phi[j][1],2) * (pow(sigma_src[3,1],2) + pow(sigma_frcH,2) + pow(sigma_omega[j],2))) + 
                       (pow(phi[j][2],2) * (pow(sigma_src[3,2],2) + pow(sigma_frcH,2) + pow(sigma_omega[j],2))) + 
                       (pow(phi[j][3],2) * (pow(sigma_src[3,3],2) + pow(sigma_frcH,2) + pow(sigma_omega[j],2))));
     }
}

model {
  int jj;
  int ii;
 
  real eps_C;
  real eps_H;
  real eps_N;

  // Priors
  dN_base ~ normal(7.03, 1.38);   // copepod 15N, measured
  dC1 ~ normal(-24.2, 0.8);     // Offshore marine, <100m Bulk POM (Wu et al 1999)
  dC2 ~ normal(-19.1, 1.2);     // Coastal, measured
  dC3 ~ normal(-26.4, 1.47);    // Freshwater POM (Geary 1988 p 80)
  dN[1] ~ normal(3.6, 0.2);       // offshore marine SPOM (Wu et al 1997, p 298)
  dN[2] ~ normal(3.1, 0.6);       // coastal, measured POM
  dN[3] ~ normal(4.4, 3.9);       // freshwater, measured POM 
  dH[1] ~ normal(-7.4, 1.0);      // Marine, measured
  dH[2] ~ normal(-15, 20);        // coastal, measured
  dH[3] ~ normal(-113.0, 10.9);   // freshwater, measured
  Delta_C[M] ~ normal(0.4, 1.3);     // 13C fractionation per trophic level (Post 2002)
  //Delta_C[M] ~ uniform(0, 1);
  Delta_N[1] ~ normal(3.0, 0.9);  //bird.liver
  Delta_N[2] ~ normal(2.2, 0.7);  //bird.blood
  Delta_N[3] ~ normal(3.2, 1.9);  //fish
  Delta_N[4] ~ normal(2.3, 0.9);  //macrozoop
  Delta_N[5] ~ normal(2.3, 0.9);  //microzoop
  Delta_N[6] ~ normal(2.2,1.1);   //YOYfish
  Delta_H ~ normal(-163.7, 27);   // 2H fractionation between water and phytoplankton
  fblki ~ uniform(30,150);
  sigma_src[1][K] ~ uniform(0.0,10000);
  sigma_src[2][K] ~ uniform(0.0,10000);
  sigma_src[3][K] ~ uniform(0.0,10000);
  sigma_frcH ~ uniform(0.0,10000);
  sigma_omega[J] ~ uniform(0.0,10000);
  sigma_fblki ~ uniform(0.0,10000);
  omega[1] ~ normal(1-pow((1-0.23), tau[1]-1), 0.02);
  //omega[1] ~ uniform(0.18, 50);
  omega[2] ~ normal(0.23, 0.03);
  omega[3] ~ normal(0.33, 0.1);
  omega[4] ~ normal(0.23, 0.03);
  omega[5] ~ normal(0.23, 0.03);
  omega[6] ~ normal(0.33, 0.1);
  omega[7] ~ normal(0.33, 0.1);
  omega[8] ~ normal(1-pow((1-0.23), tau[8]-1), 0.02);
  omega[9] ~ normal(1-pow((1-0.23), tau[9]-1), 0.02);
  //omega[8] ~ uniform(0.18, 0.50);
  //omega[9] ~ uniform(0.18, 0.50);
  omega[10] ~ normal(0.23, 0.03);
  omega[11] ~ normal(0.33, 0.1);
  omega[12] ~ normal(0.33, 0.1);
  omega[13] ~ normal(0.23, 0.03);
  omega[14] ~ normal(0.23, 0.03);
  omega[15] ~ normal(0.33, 0.1);
  omega[16] ~ normal(0.33, 0.1);
  omega[17] ~ normal(0.33, 0.1);

  for (j in 1:J) {
    sigma_frc[1,j] ~ uniform(0.0,10000);
    sigma_frc[2,j] ~ uniform(0.0,10000);
  }

  // Likelihood
  for (i in 1:N) {
    jj <- grp[i];
    eps_C <- dC_obs[i] - dC_exp[jj];
    eps_N <- dN_obs[i] - dN_exp[jj];
    eps_H <- dH_obs[i] - dH_exp[jj];
    eps_C ~ normal(0.0, sigma_C[jj]);
    eps_N ~ normal(0.0, sigma_N[jj]);
    eps_H ~ normal(0.0, sigma_H[jj]);
  }
}
generated quantities {
  int jj;
  real dH_adj[N]; // d2H in consumer due to diet
  real srcC[J,K]; // fractionation corrected carbon sources
  real srcN[J,K]; // fractionation corrected nitrogen sources
  real srcH[J,K]; // fractionation corrected hydrogen sources

 for (i in 1:26) { 
    dH_adj[i] <- (dH_obs[i] - fblki - (omega[1] * dH_w[1]))/ (1 - omega[1]);
  }
 
  for (i in 27:N) {
    jj <- grp[i];
    dH_adj[i] <- (dH_obs[i] - (omega[jj] * dH_w[jj]))/(1 - omega[jj]);
  }
  
  for (j in 1:J) {
    for (k in 1:K) {
      srcC[j,1] <- dC1 + C_tot[j];// + normal_rng(0.0, sigma_C[j]);
      srcC[j,2] <- dC2 + C_tot[j];// + normal_rng(0.0, sigma_C[j]);
      srcC[j,3] <- dC3 + C_tot[j];// + normal_rng(0.0, sigma_C[j]);
      srcN[j,k] <- dN[k] + N_tot[j];// + normal_rng(0.0, sigma_N[j]);
      srcH[j,k] <- dH[k] + Delta_H;// + normal_rng(0.0, sigma_H[j]);
    }
 }
}
