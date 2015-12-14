
#=============================================================================#
# Extract posteriors from the STAN MCMC
#=============================================================================#
trace <- extract(fit_CM)
names(trace)
post_in <- data.frame(dC1 = trace$dC1, dC2 = trace$dC2,dC3 = trace$dC3,dN = trace$dN, dH = trace$dH,
                      Delta_C = trace$Delta_C, Delta_N = trace$Delta_N, Delta_H = trace$Delta_H,frac_H_blki = trace$fblki,
                      omega = trace$omega, tau = trace$tau, phi = trace$phi, sigma_fblki = trace$sigma_fblki,
                      sigma_frc = trace$sigma_frc, sigma_frcH = trace$sigma_frcH,
                      sigma_src = trace$sigma_src,
                      lp = trace$lp__)
post <- melt(data.frame(Sample = 1:nrow(post_in), post_in), id.vars = "Sample")
head(post_in)
names(post_in)
post$Source <- NA
post$Species <- NA
for (i in 1:J)
{
    post$Source[post$variable %in%  paste0("phi.",i,".1")] <- "Offshore"
    post$Source[post$variable %in%  paste0("phi.",i,".2")] <- "Coastal"
    post$Source[post$variable %in%  paste0("phi.",i,".3")] <- "Freshwater"
    post$Species[post$variable %in% paste0("phi.",i,".1")] <- levels(factor(dat$Species, ordered = TRUE))[i]
    post$Species[post$variable %in% paste0("phi.",i,".2")] <- levels(factor(dat$Species, ordered = TRUE))[i]
    post$Species[post$variable %in% paste0("phi.",i,".3")] <- levels(factor(dat$Species, ordered = TRUE))[i]
    post$Species[post$variable %in% paste0("omega.",i)] <- levels(factor(dat$Species, ordered = TRUE))[i]
    post$Species[post$variable %in% paste0("tau.",i)] <- levels(factor(dat$Species, ordered = TRUE))[i]
    post$Species[post$variable %in% paste0("sigma_frc.1.",i)] <- levels(factor(dat$Species, ordered = TRUE))[i]
    post$Species[post$variable %in% paste0("sigma_frc.2.",i)] <- levels(factor(dat$Species, ordered = TRUE))[i]
}

# posteriors for dC_exp, dN_exp, and dH_exp
m <- as.matrix(fit_CM)
head(m)
gpost <- NULL
for (j in 1:J)
{
  block <- data.frame(levels(factor(dat$Species, ordered = TRUE))[j],
                      m[ , grep("dC_exp", colnames(m))[j]],
                      m[ , grep("dN_exp", colnames(m))[j]],
                      m[ , grep("dH_exp", colnames(m))[j]],
                      m[ , grep("sigma_C", colnames(m))[j]],
                      m[ , grep("sigma_N", colnames(m))[j]],
                      m[ , grep("sigma_H", colnames(m))[j]])
  gpost <- rbind(gpost, block)
}
gpost <- data.frame(gpost)
names(gpost) <- c("Species","dC","dN","dH", "sigma_C", "sigma_N", "sigma_H")
gpost$sp1 <- as.numeric(gpost$Species)

# posterior predictive distribution
library(plyr)
groups <- count(factor(grp))
group_sample <- matrix(groups[1:J, 2])

gpostpd_C <- data.frame(matrix(NA, nrow = max(group_sample), ncol = J))
names(gpostpd_C) <- c(levels(dat$Species))
gpostpd_N <- data.frame(matrix(NA, nrow = max(group_sample), ncol = J))
names(gpostpd_N) <- c(levels(dat$Species))
gpostpd_H <- data.frame(matrix(NA, nrow = max(group_sample), ncol = J))
names(gpostpd_H) <- c(levels(dat$Species))

for (i in 1:J) {
  gpostpd_C[1:group_sample[i], i] <- rnorm(group_sample[i], gpost$dC[gpost$sp1 == i], gpost$sigma_C[gpost$sp1 == i])
  gpostpd_N[1:group_sample[i], i] <- rnorm(group_sample[i], gpost$dN[gpost$sp1 == i], gpost$sigma_N[gpost$sp1 == i])
  gpostpd_H[1:group_sample[i], i] <- rnorm(group_sample[i], gpost$dH[gpost$sp1 == i], gpost$sigma_H[gpost$sp1 == i])
}

pdC <- melt(gpostpd_C, na.rm = T)
pdN <- melt(gpostpd_N, na.rm = T)
pdH <- melt(gpostpd_H, na.rm = T)

post_pred <- cbind(data.frame(Species = pdC[,1]), dC_ppd = pdC[,2], dN_ppd = pdN[,2], dH_ppd = pdH[,2])

# fractionation corrected source posterior
names(trace)
post2_in <- data.frame(srcC = trace$srcC, srcH = trace$srcH, srcN = trace$srcN)
post2 <- melt(data.frame(Sample = 1:nrow(post2_in), post2_in), id.vars = "Sample")
head(post2_in)
post2$Source <- NA
post2$Species <- NA
for (i in 1:J)
{
  post2$Source[post2$variable %in%  paste0("srcC.",i,".1")] <- "Offshore"
  post2$Source[post2$variable %in%  paste0("srcC.",i,".2")] <- "Coastal"
  post2$Source[post2$variable %in%  paste0("srcC.",i,".3")] <- "Freshwater"
  post2$Species[post2$variable %in% paste0("srcC.",i,".1")] <- levels(factor(dat$Species, ordered = TRUE))[i]
  post2$Species[post2$variable %in% paste0("srcC.",i,".2")] <- levels(factor(dat$Species, ordered = TRUE))[i]
  post2$Species[post2$variable %in% paste0("srcC.",i,".3")] <- levels(factor(dat$Species, ordered = TRUE))[i]
  post2$Source[post2$variable %in%  paste0("srcH.",i,".1")] <- "Offshore"
  post2$Source[post2$variable %in%  paste0("srcH.",i,".2")] <- "Coastal"
  post2$Source[post2$variable %in%  paste0("srcH.",i,".3")] <- "Freshwater"
  post2$Species[post2$variable %in% paste0("srcH.",i,".1")] <- levels(factor(dat$Species, ordered = TRUE))[i]
  post2$Species[post2$variable %in% paste0("srcH.",i,".2")] <- levels(factor(dat$Species, ordered = TRUE))[i]
  post2$Species[post2$variable %in% paste0("srcH.",i,".3")] <- levels(factor(dat$Species, ordered = TRUE))[i]
  post2$Source[post2$variable %in%  paste0("srcN.",i,".1")] <- "Offshore"
  post2$Source[post2$variable %in%  paste0("srcN.",i,".2")] <- "Coastal"
  post2$Source[post2$variable %in%  paste0("srcN.",i,".3")] <- "Freshwater"
  post2$Species[post2$variable %in% paste0("srcN.",i,".1")] <- levels(factor(dat$Species, ordered = TRUE))[i]
  post2$Species[post2$variable %in% paste0("srcN.",i,".2")] <- levels(factor(dat$Species, ordered = TRUE))[i]
  post2$Species[post2$variable %in% paste0("srcN.",i,".3")] <- levels(factor(dat$Species, ordered = TRUE))[i]
}

i <- grepl("srcC", post2$variable)
spost <- data.frame(cbind(Source = post2[i, 4], Species = post2[i, 5]))
spost$dC <- post2[i, 3]
i <- grepl("srcH", post2$variable)
spost$dH <- post2[i, 3]
i <- grepl("srcN", post2$variable)
spost$dN <- post2[i, 3]
spost$Source <- factor(spost$Source, levels = c("Offshore","Coastal","Freshwater"))
spost$Species <- factor(spost$Species, levels = c("Bulk Zoop", "Copepod", "E. pacifica", "Thysanoessa", "T. libellula", "Neomysis", "YOY Capelin", "YOY Herring", "YOY Pollock", "Eulachon", "Capelin", "Herring", "Sandlance", "Pollock", "KIMU", "MAMU", "BLKI"))

# Bayesian point estimate for dH_adj, which approximates the proportion of 2H_diet

post3_in <- data.frame(dH_adj = trace$dH_adj)
post3 <- melt(data.frame(Sample = 1:nrow(post3_in), post3_in), id.vars = "Sample")
head(post3)
dH_adj_postmean <- as.vector(by(post3[, "value"], post3$variable, mean))

#=============================================================================#
# MCMC trace plots
#=============================================================================#

i <- grepl("phi", post$variable)
d1 <- post[i,]
d1$Species <- factor(d1$Species, levels = c("Bulk Zoop", "Copepod", "E. pacifica", "Thysanoessa", "T. libellula", "Neomysis", "YOY Capelin", "YOY Herring", "YOY Pollock", "Eulachon", "Capelin", "Herring", "Sandlance", "Pollock", "KIMU", "MAMU", "BLKI"))
p <- ggplot(data = d1) + 
  facet_wrap(Species ~ Source, scales = "free_y", ncol = 6) + 
  geom_line(aes(Sample, value)) + 
  ylab(bquote(paste("Source Proportion", ~phi))) + xlab ("\nSample")
png("figs/mcmc_par_trace_phi.png", width = 10, height = 8, units = "in", res = 300)
plot(p)
dev.off()

post$Labs <- NA
post$Labs[post$variable %in%  paste0("dC.", 1)] <- "Offshore Carbon"
post$Labs[post$variable %in%  paste0("dC.", 2)] <- "Coastal Carbon"
post$Labs[post$variable %in%  paste0("dC.", 3)] <- "Freshwater Carbon"
post$Labs[post$variable %in%  paste0("dH.", 1)] <- "Offshore Hydrogen"
post$Labs[post$variable %in%  paste0("dH.", 2)] <- "Coastal Hydrogen"
post$Labs[post$variable %in%  paste0("dH.", 3)] <- "Freshwater Hydrogen"
post$Labs[post$variable %in%  paste0("dN.", 1)] <- "Offshore Nitrogen"
post$Labs[post$variable %in%  paste0("dN.", 2)] <- "Coastal Nitrogen"
post$Labs[post$variable %in%  paste0("dN.", 3)] <- "Freshwater Nitrogen"
post$Labs[post$variable %in%  "sigma_src.1.1"] <- "Sigma Offshore Carbon"
post$Labs[post$variable %in%  "sigma_src.1.2"] <- "Sigma Coastal Carbon"
post$Labs[post$variable %in%  "sigma_src.1.3"] <- "Sigma Freshwater Carbon"
post$Labs[post$variable %in%  "sigma_src.2.1"] <- "Sigma Offshore Nitrogen"
post$Labs[post$variable %in%  "sigma_src.2.2"] <- "Sigma Coastal Nitrogen"
post$Labs[post$variable %in%  "sigma_src.2.3"] <- "Sigma Freshwater Nitrogen"
post$Labs[post$variable %in%  "sigma_src.3.1"] <- "Sigma Offshore Hydrogen"
post$Labs[post$variable %in%  "sigma_src.3.2"] <- "Sigma Coastal Hydrogen"
post$Labs[post$variable %in%  "sigma_src.3.3"] <- "Sigma Freshwater Hydrogen"

for (i in 1:J) {
post$Labs[post$variable %in% paste0("tau.",i)] <- levels(factor(dat$Species, ordered = TRUE))[i]
}

for (i in 1:M) {
#post$Labs[post$variable %in%  paste0("Delta_C.",i)] <- levels(factor(dat$Taxa, ordered = TRUE))[i]
post$Labs[post$variable %in%  paste0("Delta_N.",i)] <- levels(factor(dat$Taxa, ordered = TRUE))[i]
}

post$Labs[post$variable %in%  "Delta_C"] <- "Delta C"

post$Labs[post$variable %in%  "Delta_H"] <- "Delta H"

post$Labs[post$variable %in%  "fblki"] <- "f BLKI"
# post$Labs <- factor(post$Labs, levels = c("Offshore Carbon","Coastal Carbon", "Freshwater Carbon", "Offshore Nitrogen", "Coastal Nitrogen", "Freshwater Nitrogen" , "Offshore Hydrogen", "Coastal Hydrogen", "Freshwater Hydrogen",
#                                          "Sigma Offshore Carbon","Sigma Coastal Carbon", "Sigma Freshwater Carbon", "Sigma Offshore Nitrogen", "Sigma Coastal Nitrogen", "Sigma Freshwater Nitrogen" , "Sigma Offshore Hydrogen", "Sigma Coastal Hydrogen", "Sigma Freshwater Hydrogen", "Delta C", "Delta N", "Delta H", "f BLKI"))

i <- grepl("dC", post$variable) | grepl("dN", post$variable) | grepl("dH", post$variable) | grepl("sigma_src", post$variable)
d1 <- post[i,]
head(d1)
p <- ggplot(data = d1) +
  facet_wrap(~Labs,  scales = "free_y", ncol = 3) +
  geom_line(aes(Sample, value))
png("figs/mcmc_par_trace_sources.png", width = 8, height = 8, units = "in", res = 300)
plot(p)
dev.off()

i <- grepl("Delta", post$variable) 
d1 <- post[i,]
head(d1)
p <- ggplot(data = d1) + facet_wrap(~variable,  scales = "free_y", ncol = 4) + geom_line(aes(Sample, value))
png("figs/mcmc_par_trace_fractionation.png", width = 8, height = 4, units = "in", res = 300)
plot(p)
dev.off()

i <- grepl("Delta_N", post$variable) 
d1 <- post[i,]
head(d1)
p <- ggplot(data = d1) +
  facet_wrap(~Labs,  scales = "free_y", ncol = 3) +
  geom_line(aes(Sample, value)) + ylab("\nDelta N")
png("figs/mcmc_par_trace_Nfractionation.png", width = 6, height = 4, units = "in", res = 300)
plot(p)
dev.off()


i <- grepl("sigma_frc.1.", post$variable)
d1 <- post[i,]
d1$Species <- factor(d1$Species, levels = c("Bulk Zoop", "Copepod", "E. pacifica", "Thysanoessa", "T. libellula", "Neomysis", "YOY Capelin", "YOY Herring", "YOY Pollock", "Eulachon", "Capelin", "Herring", "Sandlance", "Pollock", "KIMU", "MAMU", "BLKI"))
p <- ggplot(data = d1) +  geom_line(aes(x = Sample, y = value)) +
  facet_wrap(~Species, scales = "free_y", ncol = 4) + ylab(expression(paste(sigma, " Carbon Trophic Fractionation"))) + xlab("\nSample")
png("figs/mcmc_par_trace_sigma_C_frac.png", width = 8, height = 6, units = "in", res = 300)
plot(p)
dev.off()

i <- grepl("sigma_frc.2.", post$variable)
d1 <- post[i,]
d1$Species <- factor(d1$Species, levels = c("Bulk Zoop", "Copepod", "E. pacifica", "Thysanoessa", "T. libellula", "Neomysis", "YOY Capelin", "YOY Herring", "YOY Pollock", "Eulachon", "Capelin", "Herring", "Sandlance", "Pollock", "KIMU", "MAMU", "BLKI"))
p <- ggplot(data = d1) +  geom_line(aes(x = Sample, y = value)) +
  facet_wrap(~Species, scales = "free_y", ncol = 4) + ylab(bquote(paste("\n", sigma, " Nitrogen Trophic Fractionation"))) + xlab("\nSample")
png("figs/mcmc_par_trace_sigma_N_frac.png", width = 8, height = 6, units = "in", res = 300)
plot(p)
dev.off()

i <- grepl("omega", post$variable)
d1 <- post[i,]
d1$Species <- factor(d1$Species, levels = c("Bulk Zoop", "Copepod", "E. pacifica", "Thysanoessa", "T. libellula", "Neomysis", "YOY Capelin", "YOY Herring", "YOY Pollock", "Eulachon", "Capelin", "Herring", "Sandlance", "Pollock", "KIMU", "MAMU", "BLKI"))
p <- ggplot() + geom_line(data = d1, aes(x = Sample, y = value)) +  
  facet_wrap(~Species, scales = "free_y", ncol = 4) + ylab("Omega\n") + xlab ("\nSample")
png("figs/mcmc_par_trace_omega.png", width = 8, height = 6, units = "in", res = 300)
plot(p)
dev.off()

#=============================================================================#
# MCMC histograms
#=============================================================================#
# Find the range of each of the variables for plotting priors
rng <- apply(post_in, 2, range)
# Create a different data frame for plotting the priors
xx <- apply(rng, 2, function(x) seq(x[1], x[2], length.out = 1000))
dp <- xx
head(dp)
blki_omega <- 1-(1-0.23)^(mean(xx[, "tau.1"])-1)
kimu_omega <- 1-(1-0.23)^(mean(xx[, "tau.8"])-1)
mamu_omega <- 1-(1-0.23)^(mean(xx[, "tau.9"])-1)
pri <- c(function(x){ dnorm(x, -24.2, 0.8) },      # dC.1 marine
         function(x){ dnorm(x, -19.1, 1.2) },      # dC.2 coastal
         function(x){ dnorm(x, -26.4, 1.47) },     # dC.3 freshwater
         function(x){ dnorm(x, 3.6, 0.2) },        # dN.1 marine
         function(x){ dnorm(x, 3.1, 0.6) },        # dN.2 coastal
         function(x){ dnorm(x, 4.4, 3.9) },        # dN.3 freshwater
         function(x){ dnorm(x, -7.4, 1.0) },       # dH.1 marine
         function(x){ dnorm(x, -15.3, 3.6) },      # dH.2 coastal
         function(x){ dnorm(x, -113.0, 10.9) },    # dH.3 freshwater
         replicate(M, function(x){ dnorm(x, 0.4, 1.3) }),  # Delta_C
         #function(x){ dnorm(x, 0.4, 1.3) },      # Delta_C
         function(x){ dnorm(x, 3.0, 0.9) },      # Delta_N   
         function(x){ dnorm(x, 2.2, 0.7) },      # Delta_N
         function(x){ dnorm(x, 3.2, 1.9) },      # Delta_N
         function(x){ dnorm(x, 2.3, 0.9) },      # Delta_N
         function(x){ dnorm(x, 2.3, 0.9) },     # Delta_N
         function(x){ dnorm(x, 2.2, 1.1) },      # Delta_N
         function(x){ dnorm(x, -163.7, 27.0) },   # Delta_H
         function(x){ dunif(x, 30, 150) },         # fblki
         #function(x){ dnorm(x, blki_omega, 0.02) },# omega BLKI 
         function(x){ dunif(x, 0.18, 0.50) }, 
         function(x){ dnorm(x, 0.23, 0.03) },      # omega Bulk zoop
         function(x){ dnorm(x, 0.33, 0.1) },       # omega capelin
         function(x){ dnorm(x, 0.23, 0.03) },      # omega copepod
         function(x){ dnorm(x, 0.23, 0.03) },      # omega epacifica
         function(x){ dnorm(x, 0.33, 0.1) },       # omega eulachon
         function(x){ dnorm(x, 0.33, 0.1) },       # omega herring
         #function(x){ dnorm(x, kimu_omega, 0.02) },    # omega KIMU
         #function(x){ dnorm(x, mamu_omega, 0.02) },   # omega MAMU
         function(x){ dunif(x, 0.18, 0.50) }, 
         function(x){ dunif(x, 0.18, 0.50) }, 
         function(x){ dnorm(x, 0.23, 0.03) },  # omega neomysis
         function(x){ dnorm(x, 0.33, 0.1) },   # omega pollock
         function(x){ dnorm(x, 0.33, 0.1) },   # omega sandlance
         function(x){ dnorm(x, 0.23, 0.03) },  # omega T libellula
         function(x){ dnorm(x, 0.23, 0.03) },  # omega Thysanoessa
         function(x){ dnorm(x, 0.33, 0.1) },   # omega YOY capelin
         function(x){ dnorm(x, 0.33, 0.1) },   # omega YOY herring
         function(x){ dnorm(x, 0.33, 0.1) },   # omega YOY pollock
         replicate(J, function(x){NA}),                        # tau
         replicate(3*J, function(x){ NA }),                    # phis
         function(x){ dunif(x, 0.0, 10000)},                   # sigma_fblki
         replicate(2*J, function(x){ dunif(x, 0.0, 10000) }),  # sigma_frc
         function(x){ dunif(x, 0.0, 10000) },    # sigma_frcH
         replicate(3*K, function(x){ dunif(x, 0.0, 10000) }),  # sigma_src
         function(x){ NA }                                     # lp
)

for (ii in 1:ncol(xx)) {
    dp[,ii] <- pri[[ii]](xx[,ii])
}
xx <- melt(xx)
dp <- melt(dp)
dp$x <- xx$value
names(dp) <- c("Var1","variable","value1","x")
dp$Species <- NA
dp$Labs <- NA
for (i in 1:J)
{
  dp$Species[dp$variable %in% paste0("tau.",i)] <- levels(factor(dat$Species, ordered = TRUE))[i]
  dp$Species[dp$variable %in% paste0("phi.",i,".1")] <- levels(factor(dat$Species, ordered = TRUE))[i]
  dp$Species[dp$variable %in% paste0("phi.",i,".2")] <- levels(factor(dat$Species, ordered = TRUE))[i]
  dp$Species[dp$variable %in% paste0("phi.",i,".3")] <- levels(factor(dat$Species, ordered = TRUE))[i]
  dp$Species[dp$variable %in% paste0("omega.",i)] <- levels(factor(dat$Species, ordered = TRUE))[i]
  dp$Species[dp$variable %in% paste0("sigma_frc.1.",i)] <- levels(factor(dat$Species, ordered = TRUE))[i]
  dp$Species[dp$variable %in% paste0("sigma_frc.2.",i)] <- levels(factor(dat$Species, ordered = TRUE))[i]
  dp$Species[dp$variable %in%  paste0("sigma_frc.1", i)] <- levels(factor(dat$Species, ordered = TRUE))[i]
  dp$Species[dp$variable %in%  paste0("sigma_frc.2", i)] <- levels(factor(dat$Species, ordered = TRUE))[i]
}

for (i in 1:3){
  dp$Source[dp$variable %in%  paste0("sigma_src.", i, ".1")] <- "Offshore"
  dp$Source[dp$variable %in%  paste0("sigma_src.", i, ".2")] <- "Coastal"
  dp$Source[dp$variable %in%  paste0("sigma_src.", i, ".3")] <- "Freshwater"
  dp$Isotope[dp$variable %in%  paste0("sigma_src.1.", i)] <- "Carbon"
  dp$Isotope[dp$variable %in%  paste0("sigma_src.2.", i)] <- "Nitrogen"
  dp$Isotope[dp$variable %in%  paste0("sigma_src.3.", i)] <- "Hydrogen"
}

i <- grepl("phi", post$variable)
d1 <- post[i,]
d1$Species <- factor(d1$Species, levels = c("Bulk Zoop", "Copepod", "E. pacifica", "Thysanoessa", "T. libellula", "Neomysis", "YOY Capelin", "YOY Herring", "YOY Pollock", "Eulachon", "Capelin", "Herring", "Sandlance", "Pollock", "KIMU", "MAMU", "BLKI"))
d1$Source <- factor(d1$Source, levels = c("Offshore", "Coastal", "Freshwater"))
p <- ggplot(d1, aes(value)) +
    #facet_wrap(Species ~ Source, scales = "free", ncol = 3) +
    facet_grid(Species ~ Source, scales = "free") +
    geom_histogram(aes(y = ..density..)) +  xlab("\nSource Proportion") + ylab("Probability Density\n") +
    theme()
png("figs/mcmc_par_hist_phi.png", width = 8, height = 12, units = "in", res = 300)
plot(p)
dev.off()

i <- grepl("phi", post$variable)
d1 <- post[i,]
cpal <- c("#999999", "#E69F00", "#56B4E9")
d1$Source <- factor(d1$Source, levels = c("Offshore","Coastal","Freshwater"))
d1$Species <- factor(d1$Species, levels = c("Bulk Zoop", "Copepod", "E. pacifica", "Thysanoessa", "T. libellula", "Neomysis", "YOY Capelin", "YOY Herring", "YOY Pollock", "Eulachon", "Capelin", "Herring", "Sandlance", "Pollock", "KIMU", "MAMU", "BLKI"))
p <- ggplot(d1, aes(x = value)) + facet_wrap(~Species, scales = "free_y", ncol = 6) +
    geom_density(aes(fill = Source), alpha = 0.5) +  xlab("\nProportion of the diet") + ylab("Probability density\n") + 
    scale_fill_manual(values = cpal) + theme(legend.position = c(1,0), legend.justification = c(1,0))

png("figs/mcmc_par_density_phis.png", width = 10, height = 5, units = "in", res = 300)
plot(p)
dev.off()

i <- grepl("dC", post$variable) | grepl("dN", post$variable) | grepl("dH", post$variable) | grepl("sigma_src", post$variable)
d1 <- post[i,]
head(d1)

k <- grepl("dC", post$variable) | grepl("dN", post$variable) | grepl("dH", post$variable) | grepl("sigma_src", dp$variable)
dp2 <- dp[k,]
dp2$Labs <- NA
dp2$Labs[dp2$variable %in%  paste0("dC.", 1)] <- "Offshore Carbon"
dp2$Labs[dp2$variable %in%  paste0("dC.", 2)] <- "Coastal Carbon"
dp2$Labs[dp2$variable %in%  paste0("dC.", 3)] <- "Freshwater Carbon"
dp2$Labs[dp2$variable %in%  paste0("dH.", 1)] <- "Offshore Hydrogen"
dp2$Labs[dp2$variable %in%  paste0("dH.", 2)] <- "Coastal Hydrogen"
dp2$Labs[dp2$variable %in%  paste0("dH.", 3)] <- "Freshwater Hydrogen"
dp2$Labs[dp2$variable %in%  paste0("dN.", 1)] <- "Offshore Nitrogen"
dp2$Labs[dp2$variable %in%  paste0("dN.", 2)] <- "Coastal Nitrogen"
dp2$Labs[dp2$variable %in%  paste0("dN.", 3)] <- "Freshwater Nitrogen"
dp2$Labs[dp2$variable %in%  "sigma_src.1.1"] <- "Sigma Offshore Carbon"
dp2$Labs[dp2$variable %in%  "sigma_src.1.2"] <- "Sigma Coastal Carbon"
dp2$Labs[dp2$variable %in%  "sigma_src.1.3"] <- "Sigma Freshwater Carbon"
dp2$Labs[dp2$variable %in%  "sigma_src.2.1"] <- "Sigma Offshore Nitrogen"
dp2$Labs[dp2$variable %in%  "sigma_src.2.2"] <- "Sigma Coastal Nitrogen"
dp2$Labs[dp2$variable %in%  "sigma_src.2.3"] <- "Sigma Freshwater Nitrogen"
dp2$Labs[dp2$variable %in%  "sigma_src.3.1"] <- "Sigma Offshore Hydrogen"
dp2$Labs[dp2$variable %in%  "sigma_src.3.2"] <- "Sigma Coastal Hydrogen"
dp2$Labs[dp2$variable %in%  "sigma_src.3.3"] <- "Sigma Freshwater Hydrogen"
dp2$Labs <- factor(dp2$Labs, levels = c("Offshore Carbon","Coastal Carbon", "Freshwater Carbon", "Offshore Nitrogen", "Coastal Nitrogen", "Freshwater Nitrogen" , "Offshore Hydrogen", "Coastal Hydrogen", "Freshwater Hydrogen",
                                        "Sigma Offshore Carbon","Sigma Coastal Carbon", "Sigma Freshwater Carbon", "Sigma Offshore Nitrogen", "Sigma Coastal Nitrogen", "Sigma Freshwater Nitrogen" , "Sigma Offshore Hydrogen", "Sigma Coastal Hydrogen", "Sigma Freshwater Hydrogen"))

p <- ggplot(d1, aes(value)) +
  facet_wrap(~Labs, scales = "free", ncol = 3) +  
  geom_histogram(aes(y = ..density..))
# Add the priors
p <- p + geom_line(data = dp2, aes(x = x, y = value1), color = "red")  +  xlab(expression(paste("\nIsotope Ratio (\u2030)"))) + ylab("Probabilitiy Density\n")
png("figs/mcmc_par_hist_source.png", width = 8, height = 8, units = "in", res = 300)
plot(p)
dev.off()

i <- grepl("dC", post$variable) | grepl("dN", post$variable) | grepl("dH", post$variable)# | grepl("Delta", post$variable)
d1 <- post[i,]
p <- ggplot(d1, aes(value)) +
    facet_wrap(~Labs, scales = "free", ncol = 3) +
    geom_histogram(aes(y = ..density..)) + ylab("Probability Density\n") + xlab(expression(paste("\nIsotope Ratio (\u2030)")))
j <- grepl("dC", dp$variable) | grepl("dN", dp$variable) | grepl("dH", dp$variable)#| grepl("Delta", dp$variable)
dp1 <- dp[i,]
head(dp1)
dp1$Labs <- NA
dp1$Labs[dp1$variable %in%  paste0("dC.", 1)] <- "Offshore Carbon"
dp1$Labs[dp1$variable %in%  paste0("dC.", 2)] <- "Coastal Carbon"
dp1$Labs[dp1$variable %in%  paste0("dC.", 3)] <- "Freshwater Carbon"
dp1$Labs[dp1$variable %in%  paste0("dH.", 1)] <- "Offshore Hydrogen"
dp1$Labs[dp1$variable %in%  paste0("dH.", 2)] <- "Coastal Hydrogen"
dp1$Labs[dp1$variable %in%  paste0("dH.", 3)] <- "Freshwater Hydrogen"
dp1$Labs[dp1$variable %in%  paste0("dN.", 1)] <- "Offshore Nitrogen"
dp1$Labs[dp1$variable %in%  paste0("dN.", 2)] <- "Coastal Nitrogen"
dp1$Labs[dp1$variable %in%  paste0("dN.", 3)] <- "Freshwater Nitrogen"
#dp1$Labs[dp1$variable %in%  paste0("Delta_C")] <- "Delta C"
#dp1$Labs[dp1$variable %in%  paste0("Delta_N")] <- "Delta N"
#dp1$Labs[dp1$variable %in%  paste0("Delta_H")] <- "Delta H"

p <- p + geom_line(data = dp1, aes(x = x, y = value1), color = "red")
png("figs/mcmc_par_hist_source.png", width = 8, height = 6, units = "in", res = 300)
plot(p)
dev.off()

i <- grepl("Delta_C", post$variable)
d1 <- post[i,]
p <- ggplot(d1, aes(value)) +
  #facet_wrap(~Labs, scales = "free", ncol = 3) +
  geom_histogram(aes(y = ..density..)) + ylab("Probability Density\n") + xlab(expression(paste("\nDelta C (\u2030)")))
j <- grepl("Delta_C", dp$variable)
dp1 <- dp[i,]
head(dp1)
p <- p + geom_line(data = dp1, aes(x = x, y = value1), color = "red")
png("figs/mcmc_par_hist_DeltaC.png", width = 8, height = 6, units = "in", res = 300)
plot(p)
dev.off()

i <- grepl("Delta_N", post$variable)
d1 <- post[i,]
p <- ggplot(d1, aes(value)) +
  facet_wrap(~Labs, scales = "free", ncol = 3) +
  geom_histogram(aes(y = ..density..)) + ylab("Probability Density\n") + xlab(expression(paste("\nDelta N (\u2030)")))
j <- grepl("Delta_N", dp$variable)
dp1 <- dp[i,]
head(dp1)
for (i in 1:M){
  dp1$Labs[dp1$variable %in%  paste0("Delta_N.", i)] <- levels(factor(dat$Taxa, ordered = TRUE))[i]
}
p <- p + geom_line(data = dp1, aes(x = x, y = value1), color = "red")
png("figs/mcmc_par_hist_DeltaN.png", width = 6, height = 4, units = "in", res = 300)
plot(p)
dev.off()


i <- grepl("sigma_frc.1", post$variable)
d1 <- post[i,]
d1$Species <- factor(d1$Species, levels = c("Bulk Zoop", "Copepod", "E. pacifica", "Thysanoessa", "T. libellula", "Neomysis", "YOY Capelin", "YOY Herring", "YOY Pollock", "Eulachon", "Capelin", "Herring", "Sandlance", "Pollock", "KIMU", "MAMU", "BLKI"))
p <- ggplot(d1, aes(value)) +
  facet_wrap(~Species, scales = "free", ncol = 4) +
  geom_histogram(aes(y = ..density..)) + xlab(expression(paste(sigma, " Carbon Trophic Fractionation"))) + ylab("Probability Density\n")
# Add the priors
j <- grepl("sigma_frc.1", dp$variable)
p <- p + geom_line(data = dp[j,], aes(x = x, y = value1), color = "red")
png("figs/mcmc_par_hist_sigma_C_frac.png", width = 8, height = 6, units = "in", res = 300)
plot(p)
dev.off()

i <- grepl("sigma_frc.2", post$variable)
d1 <- post[i,]
d1$Species <- factor(d1$Species, levels = c("Bulk Zoop", "Copepod", "E. pacifica", "Thysanoessa", "T. libellula", "Neomysis", "YOY Capelin", "YOY Herring", "YOY Pollock", "Eulachon", "Capelin", "Herring", "Sandlance", "Pollock", "KIMU", "MAMU", "BLKI"))
p <- ggplot(d1, aes(value)) +
  facet_wrap(~Species, scales = "free", ncol = 4) +
  geom_histogram(aes(y = ..density..)) + xlab(expression(sigma~"Nitrogen Trophic Fractionation")) + ylab("Probability Density\n")
# Add the priors
j <- grepl("sigma_frc.2", dp$variable)
p <- p + geom_line(data = dp[j,], aes(x = x, y = value1), color = "red")
png("figs/mcmc_par_hist_sigma_N_frac.png", width = 8, height = 6, units = "in", res = 300)
plot(p)
dev.off()

i <- grepl("omega", post$variable)
d1 <- post[i,]
d1$Species <- factor(d1$Species, levels = c("Bulk Zoop", "Copepod", "E. pacifica", "Thysanoessa", "T. libellula", "Neomysis", "YOY Capelin", "YOY Herring", "YOY Pollock", "Eulachon", "Capelin", "Herring", "Sandlance", "Pollock", "KIMU", "MAMU", "BLKI"))
p <- ggplot(d1, aes(x = value)) +
  facet_wrap(~Species, scales = "free", ncol = 4) + xlab("\nOmega") + ylab("Probability Density\n") +
  geom_histogram(aes(y = ..density..))
j <- grepl("omega", dp$variable)
dp1 <- dp[j,]
dp1$Species <- factor(dp1$Species, levels = c("Bulk Zoop", "Copepod", "E. pacifica", "Thysanoessa", "T. libellula", "Neomysis", "YOY Capelin", "YOY Herring", "YOY Pollock", "Eulachon", "Capelin", "Herring", "Sandlance", "Pollock", "KIMU", "MAMU", "BLKI"))
p <- p + geom_line(data = dp1, aes(x = x, y = value1), color = "red")
png("figs/mcmc_par_hist_omega.png", width = 8, height = 6, units = "in", res = 300)
plot(p)
dev.off()

i <- grepl("tau", post$variable)
d1 <- post[i,]
d1$Species <- factor(d1$Species, levels = c("Bulk Zoop", "Copepod", "E. pacifica", "Thysanoessa", "T. libellula", "Neomysis", "YOY Capelin", "YOY Herring", "YOY Pollock", "Eulachon", "Capelin", "Herring", "Sandlance", "Pollock", "KIMU", "MAMU", "BLKI"))
p <- ggplot(d1, aes(x = value)) +
  facet_wrap(~Species, scales = "free", ncol = 4) + xlab("\nTrophic Level") + ylab("Probability Density\n") +
  geom_histogram(aes(y = ..density..))
png("figs/mcmc_par_hist_tau.png", width = 8, height = 6, units = "in", res = 300)
plot(p)
dev.off()

i <- grepl("frac_H_blki", post$variable) | grepl("sigma_fblki", post$variable)
d1 <- post[i,]
p <- ggplot(d1, aes(x = value)) +
  facet_wrap(~variable, scales = "free", ncol = 2) + ylab("Probability Density\n") + xlab(expression(paste(delta^{2}, "H (\u2030)"))) +
  geom_histogram(aes(y = ..density..))
j <- grepl("frac_H_blki", dp$variable)| grepl("sigma_fblki", dp$variable)
dp1 <- dp[j,]
p <- p + geom_line(data = dp1, aes(x = x, y = value1), color = "red")
png("figs/mcmc_par_hist_fblki.png", width = 6, height = 2, units = "in", res = 300)
plot(p)
dev.off()

# i <- grepl("DeltaC", post$variable) | grepl("frac_H_blki", post$variable)
# d1 <- post[i,]
# p <- ggplot(d1, aes(x = value)) +
#   facet_wrap(~variable, scales = "free", ncol = 4) + xlab(expression(paste("\nIsotope Ratio (\u2030)"))) + ylab("Probability Density\n") +
#   geom_histogram(aes(y = ..density..))
# j <- grepl("Delta", dp$variable) | grepl("frac_H_blki", dp$variable)
# dp1 <- dp[j,]
# p <- p + geom_line(data = dp1, aes(x = x, y = value1), color = "red")
# png("figs/mcmc_par_hist_fractionation.png", width = 8, height = 10, units = "in", res = 300)
# plot(p)
# dev.off()

#=============================================================================#
# Observed vs. expected
#=============================================================================#
labN <- expression(paste(delta^{15}, "N (\u2030)"))
labC <- expression(paste(delta^{13}, "C (\u2030)"))
labH <- expression(paste(delta^{2}, "H (\u2030)"))

gpost$Species <- factor(gpost$Species, levels = c("Bulk Zoop", "Copepod", "E. pacifica", "Thysanoessa", "T. libellula", "Neomysis", "YOY Capelin", "YOY Herring", "YOY Pollock", "Eulachon", "Capelin", "Herring", "Sandlance", "Pollock", "KIMU", "MAMU", "BLKI"))
dat$Taxa <- factor(dat$Taxa, levels = c("bird liver", "bird blood", "fish", "YOY fish", "macroz","microz"))
dat$Species <- factor(dat$Species, levels = c("Bulk Zoop", "Copepod", "E. pacifica", "Thysanoessa", "T. libellula", "Neomysis", "YOY Capelin", "YOY Herring", "YOY Pollock", "Eulachon", "Capelin", "Herring", "Sandlance", "Pollock", "KIMU", "MAMU", "BLKI"))


library(scales)
integer_breaks <- function(n = 5, ...) {
  breaker <- pretty_breaks(n, ...)
  function(x) {
    breaks <- breaker(x)
    breaks[breaks == floor(breaks)]
  }
}

# Posterior
p <- ggplot() +
    stat_density2d(data = gpost, aes(x = dC, y = dH, color = factor(Species))) +
    geom_point(data = dat, aes(x = d13C, y = d2H, color = factor(Species), shape = factor(Taxa))) + 
    xlab(labC) + ylab(labH) +
    guides(color = guide_legend(title = "Species"), shape = guide_legend(title = "Taxa"))
png("figs/fit_dC_dH.png", width = 8, height = 6, units = "in", res = 300)
plot(p)
dev.off()

p <- ggplot() +  facet_wrap(~Species, scales = "free") +
  geom_point(data = dat, aes(x = d13C, y = d2H, shape = Taxa)) +
  stat_density2d(data = gpost, aes(x = dC, y = dH, fill = ..level.., alpha = ..level..), geom = "polygon") +
  xlab(labC) + ylab(labH) + guides(fill = FALSE, alpha = FALSE)
png("figs/fit_dC_dH_species.png", width = 10, height = 6, units = "in", res = 300)
plot(p)
dev.off()

p <- ggplot() +  facet_wrap(~Species, scales = "free", ncol = 6) +
  geom_point(data = dat, aes(x = d13C, y = d2H, shape = Taxa)) +
  geom_point(data = post_pred, aes(x = dC_ppd, y = dH_ppd), color = "#CC0000", alpha = 0.3) +
  stat_density2d(data = gpost, aes(x = dC, y = dH, fill = ..level.., alpha = ..level..), geom = "polygon") +
  xlab(labC) + ylab(labH) + guides(color = guide_legend(title = "Species")) + guides(fill = FALSE, alpha = FALSE) + theme(legend.position = c(1,0), legend.justification = c(1,0), legend.text = element_text(size = 14), legend.title = element_blank()) +
  scale_x_continuous(breaks = integer_breaks()) 
png("figs/fit_dC_dH_species_ppd.png", width = 11, height = 6, units = "in", res = 300)
plot(p)
dev.off()

p <- ggplot() +
    geom_density2d(data = gpost, aes(x = dC, y = dN, color = factor(Species))) +
    geom_point(data = dat, aes(x = d13C, y = d15N, color = factor(Species), shape = Taxa)) +
    xlab(labC) + ylab(labN) + guides(color = guide_legend(title = "Species"))
png("figs/fit_dC_dN.png", width = 8, height = 6, units = "in", res = 300)
plot(p)
dev.off()

p <- ggplot() +  facet_wrap(~Species, scales = "free", ncol = 6) +
  geom_point(data = dat, aes(x = d13C, y = d15N, shape = Taxa)) +
  geom_point(data = post_pred, aes(x = dC_ppd, y = dN_ppd), color = "#CC0000", alpha = 0.3) +
  stat_density2d(data = gpost, aes(x = dC, y = dN, fill = ..level.., alpha = ..level..), geom = "polygon") +
  xlab(labC) + ylab(labN) + guides(fill = FALSE, alpha = FALSE) +  theme(legend.position = c(1,0), legend.justification = c(1,0), legend.text = element_text(size = 14), legend.title = element_blank()) +
  scale_x_continuous(breaks = integer_breaks()) 
png("figs/fit_dC_dN_species_ppd.png", width = 10, height = 6, units = "in", res = 300)
plot(p)
dev.off()

p <- ggplot() +
  geom_point(data = dat, aes(x = d13C, y = d15N, shape = Taxa)) +
  stat_density2d(data = gpost, aes(x = dC, y = dN, fill = ..level.., alpha = ..level..), geom = "polygon") +
  xlab(labC) + ylab(labN)  + facet_wrap(~Species, scales = "free") + guides(fill = FALSE, alpha = FALSE)
png("figs/fit_dC_dN_species.png", width = 10, height = 6, units = "in", res = 300)
plot(p)
dev.off()

p <- ggplot() +  facet_wrap(~Species, scales = "free") +
  geom_point(data = dat, aes(x = d15N, y = d2H, shape = Taxa)) +
  geom_point(data = post_pred, aes(x = dN_ppd, y = dH_ppd), color = "#CC0000", alpha = 0.3) +
  stat_density2d(data = gpost, aes(x = dN, y = dH, fill = ..level.., alpha = ..level..), geom = "polygon") +
  xlab(labN) + ylab(labH) + guides(fill = FALSE, alpha = FALSE)
png("figs/fit_dN_dH_species_ppd.png", width = 10, height = 6, units = "in", res = 300)
plot(p)
dev.off()

p <- ggplot() +
    geom_density2d(data = gpost, aes(x = dN, y = dH, color = factor(Species))) +
    geom_point(data = dat, aes(x = d15N, y = d2H, color = factor(Species), shape = Taxa)) +
    xlab(labN) + ylab(labH) + guides(color = guide_legend(title = "Species"))
png("figs/fit_dN_dH.png", width = 8, height = 6, units = "in", res = 300)
plot(p)
dev.off()

p <- ggplot() +
  geom_point(data = dat, aes(x = d15N, y = d2H, shape = Taxa)) +
  stat_density2d(data = gpost, aes(x = dN, y = dH, fill = ..level.., alpha = ..level..), geom = "polygon") +
  #geom_density2d(data = gpost, aes(x = dN, y = dH, color = factor(Species))) +
  xlab(labN) + ylab(labH) + facet_wrap(~Species, scales = "free") + guides(fill = FALSE, alpha = FALSE)
png("figs/fit_dN_dH_species.png", width = 10, height = 6, units = "in", res = 300)
plot(p)
dev.off()

#=============================================================================#
# Isotope mixing space plots
#=============================================================================#

p <- ggplot() + 
  stat_density2d(data = spost, aes(x = dC, y = dH, fill = Source, color = Source, geom = "polygon")) + 
  scale_color_manual(values = cpal) + guides(shape = FALSE) + theme(legend.position = c(1,0),  legend.justification = c(1,0.2)) +
  geom_point(data = dat, aes(x = d13C, y = dH_adj_postmean, shape = Taxa)) +  facet_wrap(~ Species, scales = "free", ncol = 6) + xlab(labC) + ylab(expression(paste(delta^{2},H[diet]~'(\u2030)')))  +
  scale_x_continuous(breaks = integer_breaks())
  
png("figs/fit_dC_dH_src.png", width = 10, height = 5, units = "in", res = 300)
plot(p)
dev.off()  

p <- ggplot() + 
  stat_density2d(data = spost, aes(x = dN, y = dH, fill = Source, color = Source, geom = "polygon")) +
  scale_color_manual(values = cpal) +
  geom_point(data = dat, aes(x = d15N, y = dH_adj_postmean, shape = Taxa)) +
  xlab(labN) + ylab(expression(paste(delta^{2},H[diet]~'(\u2030)'))) + facet_wrap(~Species, scales = "free") 
png("figs/fit_dN_dH_src.png", width = 9.5, height = 6, units = "in", res = 300)
plot(p)
dev.off() 

p <- ggplot() + 
  stat_density2d(data = spost, aes(x = dC, y = dN, fill = Source, color = Source, geom = "polygon")) +
  scale_color_manual(values = cpal) +
  geom_point(data = dat, aes(x = d13C, y = d15N, shape = Taxa)) +
  xlab(labC) + ylab(labN) + facet_wrap(~Species, ncol = 4, scales = "free") 
png("figs/fit_dC_dN_src.png", width =9.5, height = 6, units = "in", res = 300)
plot(p)
dev.off()

#########################################################################################################
# is there a trend in d2H_adjust relative to dN? If yes, then omega_j is not correct
# no need for any posteriors in this plot
#########################################################################################################
p <- ggplot() + 
  geom_point(data = dat, aes(x = d15N, y = dH_adj_postmean, color = factor(Species), shape = Taxa)) +
  xlab(labN) + ylab(expression(paste(delta^{2},H[diet]~'(\u2030)'))) + guides(color = guide_legend(title = "Species"))
png("figs/fit_dHdiet_dN.png", width = 8, height = 6, units = "in", res = 300)
plot(p)
dev.off()

#p <- ggplot() +
#    geom_density2d(data = gpost, aes(x = dC, y = dH)) +
#    geom_point(data = dat, aes(x = d13C, y = d2H, shape = factor(year), color = factor(fjord))) + 
#    xlab(labC) + ylab(labH) +
#    facet_wrap(~Species) + 
#    guides(color = guide_legend(title = "Fjord"),
#           shape = guide_legend(title = "Year"))

#png("figs/fit_dC_dH_v2.png", width = 8, height = 6, units = "in", res = 300)
#plot(p)
#dev.off()


#p <- ggplot() +
#    geom_density2d(data = gpost, aes(x = dN, y = dH)) +
#    geom_point(data = dat, aes(x = d15N, y = d2H, shape = factor(year), color = factor(fjord))) + 
#    xlab(labN) + ylab(labH) +
#    facet_wrap(~Species) + 
#    guides(color = guide_legend(title = "Fjord"),
#           shape = guide_legend(title = "Year"))

#png("figs/fit_dN_dH_v2.png", width = 8, height = 6, units = "in", res = 300)
#plot(p)
#dev.off()

#=============================================================================#

