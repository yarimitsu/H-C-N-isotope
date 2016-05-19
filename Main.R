# Clean up R
rm(list=ls())

# R libraries
library(rstan)
library(reshape2)
library(ggplot2)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# =====================================================================#
# DATA
# =====================================================================#
dat <- read.csv("data/glacier_iso.csv", header = TRUE)
dH_obs <- as.vector(dat$d2H)
dC_obs <- as.vector(dat$d13C)
dN_obs <- as.vector(dat$d15N)
dH_w <- dat$d2H_w
#tau <- dat$TL
sp <- factor(dat$Species, ordered = TRUE)
taxa <- c(1, 5, 3, 5, 4, 3, 3, 2, 2, 4, 3, 3, 4, 4, 6, 6, 6) #taxa-tissue type for each species group
grp <- as.numeric(sp)
sp <- unique(grp)
tx <- as.numeric(taxa)
N <- nrow(dat) 
J <- length(unique(grp))
K <- 3 # number of sources
M <- length(unique(taxa))
#tau_g <- aggregate(tau, by = list(grp), FUN = mean)[2]
#tau_g <- as.numeric(tau_g[1:J, 1]) - 1
# dN_g <- aggregate(dN_obs, by = list(grp), FUN = mean)[2]
# dN_g <- round(as.numeric(dN_g[1:J, 1]),2)
# dC_g <- aggregate(dC_obs, by = list(grp), FUN = mean)[2]
# dC_g <- round(as.numeric(dC_g[1:J, 1]),2)
# dH_w_g <- aggregate(dH_w, by = list(grp), FUN = mean)[2]
# dH_w_g <- as.numeric(dH_w_g[1:J, 1])


#load("fit_CM.RData")
#fit <- fit_CM
# =====================================================================#
# STAN MODEL (MCMC)
# =====================================================================#
# Data

isotope_dat <- list(N = N, J = J, K = K, M = M, grp = grp, 
                    dC_obs = dC_obs, dN_obs = dN_obs, dH_obs = dH_obs, dH_w_obs = dH_w, 
                    tx = tx)


isotope_init <- list(list(dN_base = 7.03, dC1 = -23.5, dC2 = -19.0,  dC3 = -24.0,
                          dN = c(  3.6,   3.2,    3.5),
                          dH = c( -7.1, -15.0, -115.0),
                          phi = matrix(c(0.2,0.6,0.2), J, K),
                          Delta_C = rep(0.5, M), Delta_N = rep(3.0, M), Delta_H = -163.0,
                          omega = rep(0.23, J), fblki = 70,
                          sigma_src = matrix(rep(0.5, 3), 3, K),
                          sigma_frc = matrix(rep(1.0, 2), 2, J),
                          sigma_frcH = 1,
                          sigma_omega = rep(0.5, J), sigma_fblki = 1.0))


# STAN model
mod <- stan_model(file = 'isotope7.stan')
# Run MCMC
warmup <- 1e3
iter <- 1e5
thin <- 50
cat((iter-warmup)/thin, "samples will be saved\n")
fit <- sampling(object = mod, data = isotope_dat, init = isotope_init,
                warmup = warmup, iter = iter, thin = thin, chains = 1)

# SIMULATION
sim_init <- list(list(dN_base = 7.03, dC1 = -23.5, dC2 = -19.0,  dC3 = -24.0,
                      dN = c(  3.6,   3.2,    3.5),
                      dH = c( -7.1, -15.0, -115.0),
                      phi = matrix(c(0.2,0.6,0.2), J, K),
                      Delta_C = rep(0.5, M), Delta_N = rep(3.0, M), Delta_H = -163.0,
                      omega = rep(0.23, J), fblki = 70,
                      sigma_src = matrix(rep(0.5, 3), 3, K),
                      sigma_frc = matrix(rep(1.0, 2), 2, J),
                      sigma_frcH = 1,
                      sigma_omega = rep(0.5, J), sigma_fblki = 1.0))
sim <- sampling(object = mod, data = isotope_dat, init = sim_init,
                iter = 1, warmup = 0, chains = 1, algorithm = "Fixed_param")
sims <- extract(sim)
sim_dat <- list(N = N, J = J, K = K, M = M, grp = grp, 
                    dC_obs = as.numeric(sims$dC_sim), dN_obs = as.numeric(sims$dN_sim), dH_obs = as.numeric(sims$dH_sim),dH_w_obs = dH_w,
                    tx = tx)
fit_sim <- sampling(object = mod, data = sim_dat, init = isotope_init,
                 warmup = warmup, iter = iter, thin = thin, chains = 1)

# Plot model fit
source("Plot.R")

#print(fit)
#write.csv(summary(fit), "f.csv")
#save(fit, file = "fit.RData")

trace <- extract(fit)
tau_post_in <- data.frame(tau = trace$tau)
tau_post <- melt(data.frame(Sample = 1:nrow(tau_post_in), tau_post_in), id.vars = "Sample")
head(tau_post)
tau_postmean <- as.vector(by(tau_post[, "value"], tau_post$variable, mean))
tau_postsd <- as.vector(by(tau_post[, "value"], tau_post$variable, sd))
tau_postm <- as.matrix(cbind(Species = levels(dat$Species), mean = round(tau_postmean,1), sd = round(tau_postsd,2)))
tau_postm

phi_post_in <- data.frame(phi = trace$phi)
phi_post <- melt(data.frame(Sample = 1:nrow(phi_post_in), phi_post_in), id.vars = "Sample")
head(phi_post)
phi_postmedian <- as.vector(by(phi_post[, "value"], phi_post$variable, median))
phi_postsd <- as.vector(by(phi_post[, "value"], phi_post$variable, sd))
phi_postm <- as.matrix(cbind(Species = rep(levels(dat$Species)), ID = levels(phi_post$variable), median = round(phi_postmedian,2), sd = round(phi_postsd,2)))
phi_postm
