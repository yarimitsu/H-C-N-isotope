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
  vector[N] dH_w_obs;   // covariate: 2H of consumer tissues
  
}
transformed data{
  int jj;
//  vector[J] n_g;
  vector[J] n_w;
//  vector[J] dN_g;   // data: mean 15N of each group
  vector[J] dH_w;   // covariate: 2H of water
  
  for(i in 1:J) {
//    n_g[i] <- 0.0;
    n_w[i] <- 0.0;
//    dN_g[i] <- 0.0;
    dH_w[i] <- 0.0;
  }
  
//  for (i in 1:N) {
//   jj <- grp[i];
//    n_g[jj] <- n_g[jj] + 1;
//    dN_g[jj] <- dN_g[jj] + dN_obs[i];
//  }

  for (i in 1:N) {
    jj <- grp[i];
    n_w[jj] <- n_w[jj] + 1;
    dH_w[jj] <- dH_w[jj] + dH_w_obs[i];
  }
  
  for (j in 1:J) {
//    dN_g[j] <- dN_g[j] / n_g[j];
    dH_w[j] <- dH_w[j] / n_w[j];
  }

  print("dH_w");
  print(dH_w);
//  print("dN_g");
//  print(dN_g);
}

parameters {
// real dN_base;                       // 15N of primary consumers
  real dN_g[J];                       // 15N of consumers by group for tau
  ordered[3] dC_ord;                   // carbon sources are ordered vector
  real dN[K];                          // 15N marine, coastal and freshwater sources
  real dH[K];                          // 2H marine, coastal and freshwater sources
  simplex[K] phi[J];                   // phi is defined as a unit simplex and thus are nonnegative and sum(phi)=1
  real phi_fw;                        // hyperparameter - fw source contribution to food web
  real<lower=0, upper = 1> Delta_C[M]; // trophic fractionation of C
  real Delta_N[M];                     // trophic fractionation of N
  real Delta_H;                        // trophic fractionation of H
  real<lower=0, upper = 0.7> omega[J]; // proportion of 2H due to ambient water dH_w
  real<lower=0> fblki;                 // blki fractionation
  real<lower=0> sigma_src[3,K];        // sigma source parameters marine, coastal, freshwater for dC, dH, dN 
  real<lower=0> sigma_frc[2,J];        // sigma C_tot, N_tot
  real<lower=0> sigma_frcH;            // sigma Delta_H
  real<lower=0> sigma_omega[J];        // sigma omega
  real<lower=0> sigma_fblki;           // sigma blki fractionation
}
transformed parameters {
  real <lower=1.8, upper = 5.5> tau[J];//trophic level by group
  real C_tot[J];
  real N_tot[J];
  real dC_exp[J];
  real dN_exp[J];
  real dH_exp[J];
  real sigma_C[J];
  real sigma_N[J];
  real sigma_H[J];
  real dC1;                            // 13C marine source
  real dC2;                            // 13C coastal source
  real dC3;                            // 13C freshwater source

  dC1 <- dC_ord[2];
  dC2 <- dC_ord[3];
  dC3 <- dC_ord[1];
  
  for (j in 1:J) {
   tau[j] <- phi[j][1]*(dN_g[j] - dN[1])/Delta_N[tx[j]] +
             phi[j][2]*(dN_g[j] - dN[2])/Delta_N[tx[j]] + 
             phi[j][3]*(dN_g[j] - dN[3])/Delta_N[tx[j]] + 1;
  }
  
//  for (i in 1:N) {
//  jj <- grp[i]
//   tau_ind[i] <- phi[jj][1]*(dN_obs[i] - dN[1])/Delta_N[tx[jj]] +
//             phi[jj][2]*(dN_obs[i] - dN[2])/Delta_N[tx[jj]] + 
//             phi[jj][3]*(dN_obs[i] - dN[3])/Delta_N[tx[jj]] + 1;
//  }
  
  
  for (j in 1:J) {
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
  int jjj;
  int ii;
 
  real eps_C;
  real eps_H;
  real eps_N;

// Priors
//  dN_base ~ normal(7.03, 1.38);   // copepod 15N, measured
  dN_g[1] ~ normal(15.91, 0.53);  // d15N BLKI
  dN_g[2] ~ normal(7.55, 1.22);   // d15N  Bulk Zoop
  dN_g[3] ~ normal(13.23, 0.41);  // d15N  Capelin
  dN_g[4] ~ normal(7.03, 01.38);  // d15N  Copepod
  dN_g[5] ~ normal(11.07, 1.21);  // d15N  Euphausia
  dN_g[6] ~ normal(14.75, 0.37);  // d15N  Eulachon
  dN_g[7] ~ normal(13.44, 0.61);  // d15N  herring
  dN_g[8] ~ normal(14.75, 0.36);  // d15N  KIMU
  dN_g[9] ~ normal(15.18, 0.35);  // d15N  MAMU
  dN_g[10] ~ normal(12.44, 0.47); // d15N  Neomysis
  dN_g[11] ~ normal(13.47, 0.75); // d15N  Pollock
  dN_g[12] ~ normal(12.29, 0.42); // d15N  sandlance
  dN_g[13] ~ normal(11.19, 0.53); // d15N  Themisto
  dN_g[14] ~ normal(9.92, 0.56);  // d15N  Thysanoessa
  dN_g[15] ~ normal(11.55, 0.28); // d15N  YOY capelin
  dN_g[16] ~ normal(11.53, 0.45); // d15N  YOY herring
  dN_g[17] ~ normal(11.67, 0.59); // d15N  YOY pollock
  
  dC1 ~ normal(-23.5, 0.8);       // Offshore marine, surface sediment 7 sites GOA, Walinski et al 2009
  dC2 ~ normal(-19.1, 1.2);       // Coastal, measured
  dC3 ~ normal(-26.4, 1.5);      // Freshwater POM (Geary 1988 p 80)
  //dN[1] ~ normal(3.6, 0.2);     // offshore marine SPOM (Wu et al 1997, p 298)
  dN[1] ~ normal(3.6, 0.7);       // offshore marine surface sediment GOA Walinski et al 2009
  dN[2] ~ normal(3.1, 0.6);       // coastal, measured POM
  dN[3] ~ normal(4.4, 3.9);       // freshwater, measured POM 
  //dH[1] ~ normal(-7.4, 1.0);    // Marine, measured
  dH[1] ~ normal(-7.4,10);        // Marine, measured *10SD
  dH[2] ~ normal(-15.3, 20);        // coastal, measured
  dH[3] ~ normal(-113.0, 10.9);   // freshwater, measured
  phi_fw ~ normal(phi[3], sigma_src[3][K]); //freshwater contribution to food web
  Delta_C[M] ~ normal(0.4, 1.3);  // 13C fractionation per trophic level (Post 2002)
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
  //omega[1] ~ normal(1-pow((1-0.23), tau[1]-1), 0.02);
  //omega[2] ~ normal(0.23, 0.03);
  //omega[3] ~ normal(0.33, 0.1);
  //omega[4] ~ normal(0.23, 0.03);
  //omega[5] ~ normal(0.23, 0.03);
  //omega[6] ~ normal(0.33, 0.1);
  //omega[7] ~ normal(0.33, 0.1);
  //omega[8] ~ normal(1-pow((1-0.23), tau[8]-1), 0.02);
  //omega[9] ~ normal(1-pow((1-0.23), tau[9]-1), 0.02);
  //omega[10] ~ normal(0.23, 0.03);
  //omega[11] ~ normal(0.33, 0.1);
  //omega[12] ~ normal(0.33, 0.1);
  //omega[13] ~ normal(0.23, 0.03);
  //omega[14] ~ normal(0.23, 0.03);
  //omega[15] ~ normal(0.33, 0.1);
  //omega[16] ~ normal(0.33, 0.1);
  //omega[17] ~ normal(0.33, 0.1);

  omega[1] ~ normal(1-pow((1-0.23), tau[1]-1), .05);
  omega[2] ~ normal(0.23, .03);
  omega[3] ~ normal(0.33, .1);
  omega[4] ~ normal(0.23, .03);
  omega[5] ~ normal(0.23, .03);
  omega[6] ~ normal(0.33, .1);
  omega[7] ~ normal(0.33, .1);
  omega[8] ~ normal(1-pow((1-0.23), tau[8]-1), .05);
  omega[9] ~ normal(1-pow((1-0.23), tau[9]-1), .05);
  omega[10] ~ normal(0.23, .03);
  omega[11] ~ normal(0.33, .1);
  omega[12] ~ normal(0.33, .1);
  omega[13] ~ normal(0.23, .03);
  omega[14] ~ normal(0.23, .03);
  omega[15] ~ normal(0.33, .1);
  omega[16] ~ normal(0.33, .1);
  omega[17] ~ normal(0.33, .1);

  for (j in 1:J) {
    sigma_frc[1,j] ~ uniform(0.0,10000);
    sigma_frc[2,j] ~ uniform(0.0,10000);
  }

  // Likelihood
  for (i in 1:N) {
    jjj <- grp[i];
    eps_C <- dC_obs[i] - dC_exp[jjj];
    eps_N <- dN_obs[i] - dN_exp[jjj];
    eps_H <- dH_obs[i] - dH_exp[jjj];
    eps_C ~ normal(0.0, sigma_C[jjj]);
    eps_N ~ normal(0.0, sigma_N[jjj]);
    eps_H ~ normal(0.0, sigma_H[jjj]);
  }
}
generated quantities {
  int jjj;
  real dH_adj[N]; // d2H in consumer due to diet
  real srcC[J,K]; // fractionation corrected carbon sources
  real srcN[J,K]; // fractionation corrected nitrogen sources
  real srcH[J,K]; // fractionation corrected hydrogen sources
  
  real eps_C;
  real eps_N;
  real eps_H;
  vector[N] dC_sim; // data: 13C of consumer tissues
  vector[N] dN_sim; // data: 15N of consumer tissues
  vector[N] dH_sim; // data: 2H of consumer tissues

  for (i in 1:N) {
    jjj <- grp[i];
    eps_C <- normal_rng(0.0, sigma_C[jjj]);
    eps_N <- normal_rng(0.0, sigma_N[jjj]);
    eps_H <- normal_rng(0.0, sigma_H[jjj]);
    dC_sim[i] <- dC_exp[jj] + eps_C;
    dN_sim[i] <- dN_exp[jj] + eps_N;
    dH_sim[i] <- dH_exp[jj] + eps_H;
  }
  
 for (i in 1:26) { 
    dH_adj[i] <- (dH_obs[i] - fblki - (omega[1] * dH_w[1]))/ (1 - omega[1]);
  }
 
  for (i in 27:N) {
    jjj <- grp[i];
    dH_adj[i] <- (dH_obs[i] - (omega[jjj] * dH_w[jjj]))/(1 - omega[jjj]);
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
