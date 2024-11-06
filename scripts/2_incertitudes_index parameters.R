rm(list=ls())

#-----------------------------------------------------------------------------------------------------#

# This script allows for the calculations of parameters used in the calculation of the I2M2,
# either as they were calculated during the development of the I2M2 (ie. Mondy et al. 2012),
# or calculated from randomized dataset (random sampling with replacement of the development dataset)

#-----------------------------------------------------------------------------------------------------#

#### Libraries ####

library(foreach)      # For parallelization
library(parallel)     # For parallelization

#-----------------------------------------------------------------------------------------------------#

#### Loading data ####

## Output directory

out.dir <- "output/"

## Additional function

source("scripts/function_calcul_DE_20181109.R") # For calculations of WGHT values

## Data used for the development of the I2M2 index

desc <- read.csv2("data/Base Invertébrés (16-02-2016)/desc_details_new.csv", check.names = F, row.names = 1) # Characteristics of the site sampling events (SSE)
met <- read.csv2("data/metrics_20190315.csv", check.names = F, row.names = 1) # Raw metric values

# Simplified typology of the streams

typo_I2M2 <- read.csv2("data/typo_I2M2_incertitudes_20191129.csv")

# Day, for file names when exporting

day <- gsub("-","",Sys.Date())

#### Run ####

#-----------------------------------------------------------------------------------------------------#

#### Initialization step ####

# Needed later in the script

# Some basic vectors

n.run = 10000
nm.pre <- colnames(desc)[32:48] # Names of the pressure categories
nm.met <- c("Shannon","ASPT","Polyvoltinism","Ovoviviparity","Richness") # Metric names
typo <- unique(met$typo)
typ <- c(T,T,F,F,T) # "Increasing" metrics (T) or "Decreasing" (F)

# Some empty tables and objects

worst <- data.frame(matrix(nrow = n.run, ncol = 5))
real.worst = c()
real.best <- data.frame(matrix(nrow = length(typo), ncol = 5))
colnames(worst) = colnames(real.best) = nm.met

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

#### BEST values ####

## As calculated when the I2M2 was developed ('real' values)
# Calculated from the LIRRs

for (i in seq(typo)) {
  met.typo <- met[which(met$typo == typo[i] & met$RR == "LIRR T"),]
  for (t in seq(typ)) {
    if (typ[t] == T) { # "Increasing" metrics
      real.best[i,t] <- quantile(met.typo[,t], probs = seq(0,1,0.05), type = 7)[20]
    } else {           # "Decreasing" metrics
      real.best[i,t] <- quantile(met.typo[,t], probs = seq(0,1,0.05), type = 7)[2]
    }
  }
}
rownames(real.best) <- typo

## Bootstrapped values, for simulating uncertainty

typo.best <- list()
N.LIRR = c()
for (i in seq(typo)) {
  print(i)
  met.typo <- met[which(met$typo == typo[i] & met$RR == "LIRR T"),]
  N.LIRR[i] <- nrow(met.typo)
  
  best <- foreach(run = seq(n.run)) %dopar% {
    id <- sample(1:nrow(met.typo), nrow(met.typo), replace = T) # Random sampling of LIRRs
    bst <- met.typo[id,]
    best.run <- vector(length = length(typ))
    for (t in seq(typ)) {
      if (typ[t] == T) { # "Increasing" metrics
        best.run[t] <- quantile(bst[,t], probs = seq(0,1,0.05), type = 7)[20]
      } else {           # "Decreasing" metrics
        best.run[t] <- quantile(bst[,t], probs = seq(0,1,0.05), type = 7)[2]
      }
    }
    best.run
  }
  best <- do.call(rbind, best)
  colnames(best) <- c("Shannon","ASPT","Polyvoltinism","Ovoviviparity","Richness")
  typo.best[[i]] <- best
}
names(typo.best) <- typo

#-----------------------------------------------------------------------------------------------------#

#### WORST values ####

## As calculated when the I2M2 was developed ('real' values)
# Calculated from the IRRs

for (t in seq(typ)) {
  if (typ[t] == T) { # "Increasing" metrics
    real.worst[t] <- quantile(met[which(met$RR == "IRR T"),t], probs = seq(0,1,0.05), type = 7)[2]
  } else {           # "Decreasing" metrics
    real.worst[t] <- quantile(met[which(met$RR == "IRR T"),t], probs = seq(0,1,0.05), type = 7)[20]
  }
}
names(real.worst) <- colnames(worst)

## Bootstrapped values, for simulating uncertainty

irr <- met[which(met$RR == "IRR T"),]

worst <- foreach(run = seq(n.run)) %dopar% {
  id <- sample(1:nrow(irr), nrow(irr), replace = T) # Random sampling of IRRs
  bst <- irr[id,]
  worst.run <- vector(length = length(typ))
  for (t in seq(typ)) {
    if (typ[t] == T) { # "Increasing" metrics
      worst.run[t] <- quantile(bst[,t], probs = seq(0,1,0.05), type = 7)[2]
    } else {           # "Decreasing" metrics
      worst.run[t] <- quantile(bst[,t], probs = seq(0,1,0.05), type = 7)[20]
    }
  }
  worst.run
}
worst <- do.call(rbind, worst)
colnames(worst) <- c("Shannon","ASPT","Polyvoltinism","Ovoviviparity","Richness")

#-----------------------------------------------------------------------------------------------------#

#### EQR calculations ####

# Used in the next section for the calculations of WGHT values

eqr <- data.frame(matrix(ncol = 5, nrow = nrow(met)))
for (i in seq(eqr)) for (j in 1:nrow(eqr)) {
  eqr[j,i] <- (met[j,i] - real.worst[i]) / (real.best[as.character(desc$typo_I2M2_V2[j]),i] - real.worst[i])
}
eqr[eqr < 0] <- 0
eqr[eqr > 1] <- 1
colnames(eqr) <- nm.met

#-----------------------------------------------------------------------------------------------------#

#### WGHT values ####

# Weight values (WGHT) = discrimination efficiency values ("de")

## As calculated when the I2M2 was developed ('real' values)

real.de <- matrix(nrow = 5, ncol = length(nm.pre))
for (p in seq(nm.pre)) {
  tmp <- eqr[which(met$RR == "LIRR T" | desc[,nm.pre[p]] %in% c("Intermediate","Poor","Bad")),]
  RR.tmp <- met$RR[which(met$RR == "LIRR T" | desc[,nm.pre[p]] %in% c("Intermediate","Poor","Bad"))]
  RR.tmp <- ifelse(RR.tmp %in% "LIRR T", "LIRR", "IRR")
  levels(RR.tmp) <- c("IRR","LIRR")
  real.de[,p] <- calcul.de(tmp, RR.tmp)[,1]
}
real.de <- data.frame(real.de)
dimnames(real.de) <- list(nm.met, nm.pre)

## Bootstrapped values, for simulating uncertainty

de <- list()

for (p in seq(nm.pre)) {
  print(p)
  print(nm.pre[p])
  tmp <- eqr[which(met$RR == "LIRR T" | desc[,nm.pre[p]] %in% c("Intermediate","Poor","Bad")),]
  RR.tmp <- met$RR[which(met$RR == "LIRR T" | desc[,nm.pre[p]] %in% c("Intermediate","Poor","Bad"))]
  tmp$RR <- RR.tmp
  tmp$RR <- ifelse(tmp$RR %in% "LIRR T", "LIRR", "IRR")
  
  de.p <- foreach(run = seq(n.run)) %dopar% {
    id <- sample(1:nrow(tmp), nrow(tmp), replace = T) # Random sampling
    bst <- tmp[id,]
    calcul.de(bst[,1:5], bst$RR)[,1]
  }
  de[[p]] <- do.call(rbind, de.p)
}

#-----------------------------------------------------------------------------------------------------#

#### Output saving ####

## Save the original and bootstrapped values of all parameters in the output directory

save(real.best, typo.best, real.worst, worst, real.de, de, file = paste0(out.dir,"index parameters_original and bootstrapped values_",day,".RData"))

#-----------------------------------------------------------------------------------------------------#
