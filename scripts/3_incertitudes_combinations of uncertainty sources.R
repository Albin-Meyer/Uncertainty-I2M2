rm(list=ls())

#-----------------------------------------------------------------------------------------------------#

# This script used the bootstrapped values of the parameters used in the calculation of the I2M2,
# (cf. script n°2) and/or Monte-Carlo simulations of the original metric values (cf. script n°1),
# before using these values in the bioindication process in order to simulate the uncertainty related to
# each of these uncertainty sources.

# Note: this script takes at least 15 hours for all calculations.
# This time could definitely be shortened by doing some more matrix calculations

#-----------------------------------------------------------------------------------------------------#

#### Libraries ####

library(foreach)      # For parallelization
library(parallel)     # For parallelization

#-----------------------------------------------------------------------------------------------------#

#### Loading data ####

## Output directory

out.dir <- "output/"

## Additional functions

source("scripts/function_quick_eqr_20190401.R", encoding = 'UTF-8')
source("scripts/function_ecol_state_20190326.R", encoding = 'UTF-8')
i2m2 <- function(y, de) {mean(apply(t(de), MAR = 1, function(x) weighted.mean(y, w = x)))}

## Data used for the development of the I2M2 index

desc <- read.csv2("data/Base Invertébrés (16-02-2016)/desc_details_new.csv", check.names = F, row.names = 1) # Characteristics of the site sampling events (SSE)
met <- read.csv2("data/metrics_20190315.csv", check.names = F, row.names = 1) # Raw metric values

# Simplified typology of the streams

typo_I2M2 <- read.csv2("data/typo_I2M2_incertitudes_20191129.csv")

# Day, for file names when exporting

day <- gsub("-","",Sys.Date())

# Loading generated random data or model parameters

load(paste0(out.dir,"index parameters_original and bootstrapped values_20241104.RData"))
load(paste0(out.dir,"sigmet_20241029.RData"))

#### Run ####

#-----------------------------------------------------------------------------------------------------#

#### Initialization step ####

# Needed later in the script

# Some basic vectors

# n.run = 10000
n.run = 100
nm.met <- c("Shannon","ASPT","Polyvoltinism","Ovoviviparity","Richness") # Metric names
typo <- unique(met$typo)

# Seed for random processes
# Of note, I did not use a given seed when I did the analyses presented in the associated paper

set.seed(123)

## Parameters for parallelization

# Number of logical cores

n.cores <- parallel::detectCores() - 2

# Local cluster

my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK") # PSOCK for Windows (the other one for Linux)
print(my.cluster)

# Cluster is registered for %dopar%

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered() # Should be TRUE
foreach::getDoParWorkers() # Should be equal to the number of cores

#-----------------------------------------------------------------------------------------------------#

#### EQR calculations ####

# Used later for uncertainty analyses when they do not involve inter-operator variability

eqr <- data.frame(matrix(ncol = 5, nrow = nrow(met)))
for (i in seq(eqr)) for (j in 1:nrow(eqr)) {
  eqr[j,i] <- (met[j,i] - real.worst[i]) / (real.best[as.character(desc$typo_I2M2_V2[j]),i] - real.worst[i])
}
eqr[eqr < 0] <- 0
eqr[eqr > 1] <- 1
colnames(eqr) <- nm.met

#-----------------------------------------------------------------------------------------------------#

#### Descriptive statistics ####

# Numbers of LIRRs per stream type

N.LIRR <- c()
for (i in seq(typo)) {
  met.typo <- met[which(met$typo == typo[i] & met$RR == "LIRR T"),]
  N.LIRR[i] <- nrow(met.typo)
}
range(N.LIRR)
summary(N.LIRR)

# Global numbers of IRRs and LIRRs

sum(met$RR == "IRR T")
sum(met$RR == "LIRR T")

#-----------------------------------------------------------------------------------------------------#

##### Single uncertainty sources #####

#-----------------------------------------------------------------------------------------------------#

#### WGHT ####

## Calculations of the uncertainties linked only to WGHT values

ind <- foreach(run = seq(n.run), .combine = cbind) %dopar% {
  ghost.de <- real.de
  for (i in seq(ghost.de)) {
    ghost.de[,i] <- de[[i]][sample(1:10000, 1),] # Random sampling among the bootstrapped values
  }
  ind.run <- apply(eqr, MAR = 1, i2m2, de = ghost.de) # Calculation of the I2M2 value
  ind.run
}

eco.ind <- foreach(i = 1:nrow(ind), .combine = rbind) %dopar% {
  eco.ind.tmp <- ecol.state(ind[i,], desc[i,], typo_I2M2) # Assignment of the ecological class to each SSE
  eco.ind.tmp
}

save(eco.ind, file = paste0(out.dir,"uncertainty_WGHT_eco_",day,".RData"))

#-----------------------------------------------------------------------------------------------------#

#### BEST ####

## Calculations of the uncertainties linked only to BEST values

ind.b <- foreach(run = seq(n.run), .combine = cbind) %dopar% {
  ghost.best <- real.best
  for (i in 1:nrow(ghost.best)) {
    ghost.best[i,] <- typo.best[[rownames(ghost.best)[i]]][sample(1:10000, 1),]
  }
  # EQR calculations
  ghost.eqr <- quick.eqr(met, desc, nm.met = nm.met, real.worst = real.worst, real.best = ghost.best)
  # I2M2 calculations
  ind.b.run <- apply(ghost.eqr, MAR = 1, i2m2, de = real.de)
  ind.b.run
}

eco.b <- foreach(i = 1:nrow(ind.b), .combine = rbind) %dopar% {
  eco.ind.tmp <- ecol.state(ind.b[i,], desc[i,], typo_I2M2) # Assignment of the ecological class to each SSE
  eco.ind.tmp
}

save(eco.b, file = paste0(out.dir,"uncertainty_BEST_eco_",day,".RData"))

# Analysis of the results per stream type

# N.LIRR = N.OPC = N.INCER = INCER = N.GM = GM = c()
# for (i in seq(typo)) {
#   # Number of LIRRs per type
#   met.typo <- met[which(met$typo == typo[i] & met$RR == "LIRR T"),]
#   N.LIRR[i] <- nrow(met.typo)
#   # Other analyses
#   tmp <- results[which(desc$typo_I2M2_V2 == typo[i]),]
#   N.OPC[i] <- nrow(tmp) # Number of SSEs per stream
#   N.INCER[i] <- sum(tmp$incer30 == 2)
#   INCER[i] <- sum(tmp$incer30 == 2)/nrow(tmp)
#   N.GM[i] <- sum(tmp$incer30 == 2 & tmp$mn.eco > 2 & tmp$mn.eco < 3)
#   GM[i] <- sum(tmp$incer30 == 2 & tmp$mn.eco > 2 & tmp$mn.eco < 3)/nrow(tmp)
# }
# 
# final <- data.frame(typo, N.OPC, N.LIRR, N.INCER, INCER, N.GM, GM)
# write.csv2(final, file = paste0("output/uncertainty_BEST_per stream type_",day,".csv"))

#-----------------------------------------------------------------------------------------------------#

#### WORST ####

## Calculations of the uncertainties linked only to WORST values

ind.w <- foreach(run = seq(n.run), .combine = cbind) %dopar% {
  ghost.worst <- unlist(worst[sample(1:10000, 1),])
  # EQR calculations
  ghost.eqr <- quick.eqr(met, desc, nm.met = nm.met, real.worst = ghost.worst, real.best = real.best)
  # I2M2 calculations
  ind.w.run <- apply(ghost.eqr, MAR = 1, i2m2, de = real.de)
  ind.w.run
}

eco.w <- foreach(i = 1:nrow(ind.w), .combine = rbind) %dopar% {
  eco.ind.tmp <- ecol.state(ind.w[i,], desc[i,], typo_I2M2) # Assignment of the ecological class to each SSE
  eco.ind.tmp
}

save(eco.w, file = paste0(out.dir,"uncertainty_WORST_eco_",day,".RData"))

#-----------------------------------------------------------------------------------------------------#

#### OBS ####

## Calculations of the uncertainties linked only to the inter-operator variability of metric values

ind.m <- foreach(run = seq(n.run), .combine = cbind) %dopar% {
  # METRIQUES
  ghost.met <- met
  for (i in seq(sigmet)) {
    ghost.met[,i] <- sapply(met[,i], function(x) rnorm(n = 1, mean = x, sd = sigmet[i]))
  }
  # EQR calculations
  ghost.eqr <- quick.eqr(ghost.met, desc, nm.met = nm.met, real.worst = real.worst, real.best = real.best)
  # I2M2 calculations
  ind.m.run <- apply(ghost.eqr, MAR = 1, i2m2, de = real.de)
  ind.m.run
}

eco.m <- foreach(i = 1:nrow(ind.m), .combine = rbind) %dopar% {
  eco.tmp <- ecol.state(ind.m[i,], desc[i,], typo_I2M2) # Assignment of the ecological class to each SSE
  eco.tmp
}

save(eco.m, file = paste0(out.dir,"uncertainty_OBS_eco_",day,".RData"))

#-----------------------------------------------------------------------------------------------------#

##### Dual uncertainty sources #####

#-----------------------------------------------------------------------------------------------------#

#### WGHT + BEST ####

## Calculations of the uncertainties linked to the combination
## of uncertainties linked to WGHT and BEST values

ind.bde <- foreach(run = seq(n.run), .combine = cbind) %dopar% {
  # BEST
  ghost.best <- real.best
  for (i in 1:nrow(ghost.best)) {
    ghost.best[i,] <- typo.best[[rownames(ghost.best)[i]]][sample(1:10000, 1),]
  }
  # WGHT
  ghost.de <- real.de
  for (i in seq(ghost.de)) {
    ghost.de[,i] <- de[[i]][sample(1:10000, 1),]
  }
  # EQR calculations
  ghost.eqr <- quick.eqr(met, desc, nm.met = nm.met, real.worst = real.worst, real.best = ghost.best)
  # I2M2 calculations
  ind.bde.run <- apply(ghost.eqr, MAR = 1, i2m2, de = ghost.de)
  ind.bde.run
}

eco.bde <- foreach(i = 1:nrow(ind.bde), .combine = rbind) %dopar% {
  eco.tmp <- ecol.state(ind.bde[i,], desc[i,], typo_I2M2) # Assignment of the ecological class to each SSE
  eco.tmp
}

save(eco.bde, file = paste0(out.dir,"uncertainty_BWg_eco_",day,".RData"))

#-----------------------------------------------------------------------------------------------------#

#### WGHT + WORST ####

## Calculations of the uncertainties linked to the combination
## of uncertainties linked to WGHT and WORST values

ind.wde <- foreach(run = seq(n.run), .combine = cbind) %dopar% {
  # WORST
  ghost.worst <- unlist(worst[sample(1:10000, 1),])
  # WGHT
  ghost.de <- real.de
  for (i in seq(ghost.de)) {
    ghost.de[,i] <- de[[i]][sample(1:10000, 1),]
  }
  # EQR calculations
  ghost.eqr <- quick.eqr(met, desc, nm.met = nm.met, real.worst = ghost.worst, real.best = real.best)
  # I2M2 calculations
  ind.wde.run <- apply(ghost.eqr, MAR = 1, i2m2, de = ghost.de)
  ind.wde.run
}

eco.wde <- foreach(i = 1:nrow(ind.wde), .combine = rbind) %dopar% {
  eco.tmp <- ecol.state(ind.wde[i,], desc[i,], typo_I2M2) # Assignment of the ecological class to each SSE
  eco.tmp
}

save(eco.wde, file = paste0(out.dir,"uncertainty_WoWg_eco_",day,".RData"))

#-----------------------------------------------------------------------------------------------------#

#### WGHT + OBS ####

## Calculations of the uncertainties linked to the combination
## of uncertainties linked to WGHT values and inter-operator variability

ind.mde <- foreach(run = seq(n.run), .combine = cbind) %dopar% {
  # OBS (inter-operator variability)
  ghost.met <- met
  for (i in seq(sigmet)) {
    ghost.met[,i] <- sapply(met[,i], function(x) rnorm(n = 1, mean = x, sd = sigmet[i]))
  }
  # WGHT
  ghost.de <- real.de
  for (i in seq(ghost.de)) {
    ghost.de[,i] <- de[[i]][sample(1:10000, 1),]
  }
  # EQR calculations
  ghost.eqr <- quick.eqr(ghost.met, desc, nm.met = nm.met, real.worst = real.worst, real.best = real.best)
  # I2M2 calculations
  ind.mde.run <- apply(ghost.eqr, MAR = 1, i2m2, de = ghost.de)
  ind.mde.run
}

eco.mde <- foreach(i = 1:nrow(ind.mde), .combine = rbind) %dopar% {
  eco.tmp <- ecol.state(ind.mde[i,], desc[i,], typo_I2M2) # Assignment of the ecological class to each SSE
  eco.tmp
}

save(eco.mde, file = paste0(out.dir,"uncertainty_OBSWg_eco_",day,".RData"))

#-----------------------------------------------------------------------------------------------------#

#### BEST + WORST ####

## Calculations of the uncertainties linked to the combination
## of uncertainties linked to BEST and WORST values

ind.bw <- foreach(run = seq(n.run), .combine = cbind) %dopar% {
  # BEST
  ghost.best <- real.best
  for (i in 1:nrow(ghost.best)) {
    ghost.best[i,] <- typo.best[[rownames(ghost.best)[i]]][sample(1:10000, 1),]
  }
  # WORST
  ghost.worst <- unlist(worst[sample(1:10000, 1),])
  # EQR calculations
  ghost.eqr <- quick.eqr(met, desc, nm.met = nm.met, real.worst = ghost.worst, real.best = ghost.best)
  # I2M2 calculations
  ind.bw.run <- apply(ghost.eqr, MAR = 1, i2m2, de = real.de)
  ind.bw.run
}

eco.bw <- foreach(i = 1:nrow(ind.bw), .combine = rbind) %dopar% {
  eco.tmp <- ecol.state(ind.bw[i,], desc[i,], typo_I2M2) # Assignment of the ecological class to each SSE
  eco.tmp
}

save(eco.bw, file = paste0(out.dir,"uncertainty_BWo_eco_",day,".RData"))

#-----------------------------------------------------------------------------------------------------#

#### BEST + OBS ####

## Calculations of the uncertainties linked to the combination
## of uncertainties linked to BEST values and inter-operator variability

ind.mb <- foreach(run = seq(n.run), .combine = cbind) %dopar% {
  # OBS (inter-operator variability)
  ghost.met <- met
  for (i in seq(sigmet)) {
    ghost.met[,i] <- sapply(met[,i], function(x) rnorm(n = 1, mean = x, sd = sigmet[i]))
  }
  # BEST
  ghost.best <- real.best
  for (i in 1:nrow(ghost.best)) {
    ghost.best[i,] <- typo.best[[rownames(ghost.best)[i]]][sample(1:10000, 1),]
  }
  # EQR calculations
  ghost.eqr <- quick.eqr(ghost.met, desc, nm.met = nm.met, real.worst = real.worst, real.best = ghost.best)
  # I2M2 calculations
  ind.mb.run <- apply(ghost.eqr, MAR = 1, i2m2, de = real.de)
  ind.mb.run
}

eco.mb <- foreach(i = 1:nrow(ind.mb), .combine = rbind) %dopar% {
  eco.tmp <- ecol.state(ind.mb[i,], desc[i,], typo_I2M2) # Assignment of the ecological class to each SSE
  eco.tmp
}

save(eco.mb, file = paste0(out.dir,"uncertainty_OBSB_eco_",day,".RData"))

#-----------------------------------------------------------------------------------------------------#

#### WORST + OBS ####

## Calculations of the uncertainties linked to the combination
## of uncertainties linked to WORST values and inter-operator variability

ind.mw <- foreach(run = seq(n.run), .combine = cbind) %dopar% {
  # OBS (inter-operator variability)
  ghost.met <- met
  for (i in seq(sigmet)) {
    ghost.met[,i] <- sapply(met[,i], function(x) rnorm(n = 1, mean = x, sd = sigmet[i]))
  }
  # WORST
  ghost.worst <- unlist(worst[sample(1:10000, 1),])
  # EQR calculations
  ghost.eqr <- quick.eqr(ghost.met, desc, nm.met = nm.met, real.worst = ghost.worst, real.best = real.best)
  # I2M2 calculations
  ind.mw.run <- apply(ghost.eqr, MAR = 1, i2m2, de = real.de)
  ind.mw.run
}

eco.mw <- foreach(i = 1:nrow(ind.mw), .combine = rbind) %dopar% {
  eco.tmp <- ecol.state(ind.mw[i,], desc[i,], typo_I2M2) # Assignment of the ecological class to each SSE
  eco.tmp
}

save(eco.mw, file = paste0(out.dir,"uncertainty_OBSWo_eco_",day,".RData"))

#-----------------------------------------------------------------------------------------------------#

##### Triple uncertainty sources #####

#-----------------------------------------------------------------------------------------------------#

#### WGHT + BEST + WORST ####

## Calculations of the uncertainties linked to the combination
## of uncertainties linked to WGHT, BEST and WORST values

ind.bwde <- foreach(run = seq(n.run), .combine = cbind) %dopar% {
  # BEST
  ghost.best <- real.best
  for (i in 1:nrow(ghost.best)) {
    ghost.best[i,] <- typo.best[[rownames(ghost.best)[i]]][sample(1:10000, 1),]
  }
  # WORST
  ghost.worst <- unlist(worst[sample(1:10000, 1),])
  # WGHT
  ghost.de <- real.de
  for (i in seq(ghost.de)) {
    ghost.de[,i] <- de[[i]][sample(1:10000, 1),]
  }
  # EQR calculations
  ghost.eqr <- quick.eqr(met, desc, nm.met = nm.met, real.worst = ghost.worst, real.best = ghost.best)
  # I2M2 calculations
  ind.bwde.run <- apply(ghost.eqr, MAR = 1, i2m2, de = ghost.de)
  ind.bwde.run
}

eco.bwde <- foreach(i = 1:nrow(ind.bwde), .combine = rbind) %dopar% {
  eco.tmp <- ecol.state(ind.bwde[i,], desc[i,], typo_I2M2) # Assignment of the ecological class to each SSE
  eco.tmp
}

save(eco.bwde, file = paste0(out.dir,"uncertainty_BWoWg_eco_",day,".RData"))

#-----------------------------------------------------------------------------------------------------#

#### WGHT + BEST + OBS ####

## Calculations of the uncertainties linked to the combination
## of uncertainties linked to WGHT and BEST values and to inter-operator variability

ind.mbde <- foreach(run = seq(n.run), .combine = cbind) %dopar% {
  # OBS (inter-operator variability)
  ghost.met <- met
  for (i in seq(sigmet)) {
    ghost.met[,i] <- sapply(met[,i], function(x) rnorm(n = 1, mean = x, sd = sigmet[i]))
  }
  # BEST
  ghost.best <- real.best
  for (i in 1:nrow(ghost.best)) {
    ghost.best[i,] <- typo.best[[rownames(ghost.best)[i]]][sample(1:10000, 1),]
  }
  # WGHT
  ghost.de <- real.de
  for (i in seq(ghost.de)) {
    ghost.de[,i] <- de[[i]][sample(1:10000, 1),]
  }
  # EQR calculations
  ghost.eqr <- quick.eqr(ghost.met, desc, nm.met = nm.met, real.worst = real.worst, real.best = ghost.best)
  # I2M2 calculations
  ind.mbde.run <- apply(ghost.eqr, MAR = 1, i2m2, de = ghost.de)
  ind.mbde.run
}

eco.mbde <- foreach(i = 1:nrow(ind.mbde), .combine = rbind) %dopar% {
  eco.tmp <- ecol.state(ind.mbde[i,], desc[i,], typo_I2M2) # Assignment of the ecological class to each SSE
  eco.tmp
}

save(eco.mbde, file = paste0(out.dir,"uncertainty_OBSBWg_eco_",day,".RData"))

#-----------------------------------------------------------------------------------------------------#

#### WGHT + WORST + OBS ####

## Calculations of the uncertainties linked to the combination
## of uncertainties linked to WGHT and WORST values and to inter-operator variability

ind.mwde <- foreach(run = seq(n.run), .combine = cbind) %dopar% {
  # OBS (inter-operator variability)
  ghost.met <- met
  for (i in seq(sigmet)) {
    ghost.met[,i] <- sapply(met[,i], function(x) rnorm(n = 1, mean = x, sd = sigmet[i]))
  }
  # WORST
  ghost.worst <- unlist(worst[sample(1:10000, 1),])
  # WGHT
  ghost.de <- real.de
  for (i in seq(ghost.de)) {
    ghost.de[,i] <- de[[i]][sample(1:10000, 1),]
  }
  # EQR calculations
  ghost.eqr <- quick.eqr(ghost.met, desc, nm.met = nm.met, real.worst = ghost.worst, real.best = real.best)
  # I2M2 calculations
  ind.mwde.run <- apply(ghost.eqr, MAR = 1, i2m2, de = ghost.de)
  ind.mwde.run
}

eco.mwde <- foreach(i = 1:nrow(ind.mwde), .combine = rbind) %dopar% {
  eco.tmp <- ecol.state(ind.mwde[i,], desc[i,], typo_I2M2) # Assignment of the ecological class to each SSE
  eco.tmp
}

save(eco.mwde, file = paste0(out.dir,"uncertainty_OBSWoWg_eco_",day,".RData"))

#-----------------------------------------------------------------------------------------------------#

#### BEST + WORST + OBS ####

## Calculations of the uncertainties linked to the combination
## of uncertainties linked to BEST and WORST values and to inter-operator variability

ind.mbw <- foreach(run = seq(n.run), .combine = cbind) %dopar% {
  # OBS (inter-operator variability)
  ghost.met <- met
  for (i in seq(sigmet)) {
    ghost.met[,i] <- sapply(met[,i], function(x) rnorm(n = 1, mean = x, sd = sigmet[i]))
  }
  # BEST
  ghost.best <- real.best
  for (i in 1:nrow(ghost.best)) {
    ghost.best[i,] <- typo.best[[rownames(ghost.best)[i]]][sample(1:10000, 1),]
  }
  # WORST
  ghost.worst <- unlist(worst[sample(1:10000, 1),])
  # EQR calculations
  ghost.eqr <- quick.eqr(ghost.met, desc, nm.met = nm.met, real.worst = ghost.worst, real.best = ghost.best)
  # I2M2 calculations
  ind.mbw.run <- apply(ghost.eqr, MAR = 1, i2m2, de = real.de)
  ind.mbw.run
}

eco.mbw <- foreach(i = 1:nrow(ind.mbw), .combine = rbind) %dopar% {
  eco.tmp <- ecol.state(ind.mbw[i,], desc[i,], typo_I2M2) # Assignment of the ecological class to each SSE
  eco.tmp
}

save(eco.mbw, file = paste0(out.dir,"uncertainty_OBSBWo_eco_",day,".RData"))

#-----------------------------------------------------------------------------------------------------#

##### All four sources of uncertainty #####

#-----------------------------------------------------------------------------------------------------#

#### WGHT + BEST + WORST + OBS #### 

## Calculations of the uncertainties linked to the combination
## of uncertainties linked to WGHT, BEST and WORST values and to inter-operator variability

ind.mbwde <- foreach(run = seq(n.run), .combine = cbind) %dopar% {
  # OBS (inter-operator variability)
  ghost.met <- met
  for (i in seq(sigmet)) {
    ghost.met[,i] <- sapply(met[,i], function(x) rnorm(n = 1, mean = x, sd = sigmet[i]))
  }
  # BEST
  ghost.best <- real.best
  for (i in 1:nrow(ghost.best)) {
    ghost.best[i,] <- typo.best[[rownames(ghost.best)[i]]][sample(1:10000, 1),]
  }
  # WORST
  ghost.worst <- unlist(worst[sample(1:10000, 1),])
  # WGHT
  ghost.de <- real.de
  for (i in seq(ghost.de)) {
    ghost.de[,i] <- de[[i]][sample(1:10000, 1),]
  }
  # EQR calculations
  ghost.eqr <- quick.eqr(ghost.met, desc, nm.met = nm.met, real.worst = ghost.worst, real.best = ghost.best)
  # I2M2 calculations
  ind.mbwde.run <- apply(ghost.eqr, MAR = 1, i2m2, de = ghost.de)
  ind.mbwde.run
}

eco.mbwde <- foreach(i = 1:nrow(ind.mbwde), .combine = rbind) %dopar% {
  eco.tmp <- ecol.state(ind.mbwde[i,], desc[i,], typo_I2M2) # Assignment of the ecological class to each SSE
  eco.tmp
}

save(eco.mbwde, file = paste0(out.dir,"uncertainty_OBSWoWgB_eco_",day,".RData"))

# Stats, for Figure 7

sd.eqr <- apply(ind.mbwde, MAR = 1, sd) # Standard deviation of all EQR values, per SSE, for Figure 7
save(sd.eqr, file = paste0(out.dir,"sd_OBSWoWgB_",day,".RData"))

#-----------------------------------------------------------------------------------------------------#



