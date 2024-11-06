rm(list=ls())

#-----------------------------------------------------------------------------------------------------#

# This final script does some final data preparation before generating all remaining figures.

#-----------------------------------------------------------------------------------------------------#

#### Librairies ####

library(ggplot2)
library(ggupset)
library(nls2)
library(conover.test)

#-----------------------------------------------------------------------------------------------------#

#### Loading data ####

## Output directory

out.dir <- "output/"

## Additional functions

source("scripts/function_quick_eqr_20190401.R", encoding = 'UTF-8')
source("scripts/function_ecol_state_20190326.R", encoding = 'UTF-8')
i2m2 <- function(y, de) {mean(apply(t(de), MAR = 1, function(x) weighted.mean(y, w = x)))}

## Chargement des tableaux

desc <- read.csv2("data/Base Invertébrés (16-02-2016)/desc_details_new.csv", check.names = F, row.names = 1)
met <- read.csv2("data/metrics_20190315.csv", check.names = F, row.names = 1)

# Simplified typology of the streams

typo_I2M2 <- read.csv2("data/typo_I2M2_incertitudes_20191129.csv")
typo_I2M2$lim.cat <- factor(typo_I2M2$PB) # Using Poor-Bad boundary as an indicator of intercal type for streams
levels(typo_I2M2$lim.cat) <- c("cat5","cat1","cat4","cat3","cat2") # Simplified intercal categories (types)

# Boundaries of the ecological quality classes

limits <- data.frame(lim.PB = c(0.118,0.148,0.153,0.155,0.166), lim.MP = c(0.236,0.295,0.306,0.31,0.332),
                      lim.GM = c(0.354,0.443,0.46,0.464,0.498), lim.HG = c(0.605,0.665,0.665,0.676,0.665))
rownames(limits) <- c("cat5","cat1","cat4","cat3","cat2")

# Loading generated random data or model parameters

load(paste0(out.dir,"index parameters_original and bootstrapped values_20241104.RData"))
load(paste0(out.dir,"sigmet_20241029.RData"))

# Loading standard deviations of EQR values, for the combination of all uncertainty sources

load(paste0(out.dir,"sd_OBSWoWgB_20241105.RData"))

# Day, for file names when exporting

day <- gsub("-","",Sys.Date())

#### Run ####

#-----------------------------------------------------------------------------------------------------#

#### Initialization step ####

## Needed later in the script

# Some basic vectors

n.run <- 10000
nm.met <- c("Shannon","ASPT","Polyvoltinism","Ovoviviparity","Richness") # Metric names
typo <- unique(met$typo)
nm.ce <- c("Bad","Poor","Moderate","Good","High") # Names of ecological quality classes

# Numbers of LIRR sites and SSEs per stream type

N.LIRR = N.OPC = c()
for (i in seq(typo)) {
  met.typo <- met[which(met$typo == typo[i] & met$RR == "LIRR T"),]
  desc.typo <- desc[which(met$typo == typo[i] & met$RR == "LIRR T"),]
  N.LIRR[i] <- nrow(met.typo)
  N.OPC[i] <- length(unique(desc.typo$cd_site))
}
names(N.LIRR) = names(N.OPC) = typo

#-------------------------------------------------------------------------------------------------------#

#### Original quality classes ####

## Original I2M2 values (EQR)

ghost.eqr <- quick.eqr(met, desc, nm.met = nm.met, real.worst = real.worst, real.best = real.best)
eqr2 <- apply(ghost.eqr, MAR = 1, i2m2, de = real.de)

## Assignment of the ecological class to each SSE

cee.opc <- c()
for (i in seq(eqr2)) {
  cee.opc[i] <- ecol.state(eqr2[i], desc[i,], typo_I2M2)
}
cee.opc <- factor(cee.opc)
levels(cee.opc) <- 5:1 # From 5 = Bad, to 1 = High

#-------------------------------------------------------------------------------------------------------#

#### Mismatch Frequency (MF) ####

## Relative frequency of quality class mismatches

# Chargement des données de classes d'état d'écologique

sel <- list.files(out.dir, pattern = "uncertainty_")

# Names of uncertainty sources

src.incer <- gsub("uncertainty_", "", sel)
src.incer <- gsub("_eco_[0-9]*.RData$", "", src.incer)
src.incer

# Calculation loop

eco.count <- matrix(nrow = nrow(desc), ncol = length(sel), data = NA)
for (i in seq(sel)) { # i = 12 for the combination of all sources of uncertainty
  
  # Data loading
  
  eco.data <- load(paste0(out.dir, sel[i]))
  eco <- get(eco.data)
  rm(list = eco.data)
  
  # Percentage of a given SSE to be in each ecological class
  
  zeros <- rep(0, nrow(eco))
  eco.count.tmp <- data.frame(V1 = zeros, V2 = zeros, V3 = zeros, V4 = zeros, V5 = zeros)
  for (j in 1:5) { # Here, the j index correspond to the "number" of the quality class
    eco.count.tmp[,j] <- apply(eco, MAR = 1, function(x) sum(x[!is.na(x)] == j)/ncol(eco))
  }
  
  # Mismatch Frequency (MF)
  
  p.cee = c()
  for (j in 1:nrow(eco.count.tmp)) {
    p.cee[j] <- eco.count.tmp[j,cee.opc[j]]
  }
  eco.count.tmp$mf <- 1 - p.cee
  
  # Output MF in main table
  
  eco.count[,i] <- eco.count.tmp$mf
  
}
eco.count <- as.data.frame(eco.count)
colnames(eco.count) <- src.incer

#-----------------------------------------------------------------------------------------------------#

#### Fig.3 ####

# Average MF per combination of uncertainty sources

mean.mf <- apply(eco.count, MAR = 2, mean)
mf.long <- data.frame(mf = mean.mf[1:15])

combination <- c("BEST", "BEST-WGHT", "WORST-BEST", "WORST-BEST-WGHT", "OBS", "OBS-BEST", "OBS-BEST-WGHT", "OBS-WORST-BEST",
                 "OBS-WGHT", "OBS-WORST", "OBS-WORST-WGHT", "OBS-WORST-BEST-WGHT", "WGHT", "WORST", "WORST-WGHT")

mf.long$combination <- combination
mf.long <- mf.long[c(5,1,14,13,6,10,9,3,2,15,8,7,11,4,12),]
mf.long$combination <- strsplit(mf.long$combination, "-")

# Upset-plot

svg(paste0(out.dir,"Fig3_Average MF_",day,".svg"))
ggplot(mf.long) +
  aes(x = combination, y = mf, label = sprintf("%.3f", mf)) +
  geom_col() +
  scale_x_upset(order_by = "degree") +
  ylim(0,0.15) +
  geom_text(nudge_y = 0.005) +
  labs(x = "", y = "Average MF") +
  theme_light() +
  theme(panel.grid.major.x = element_blank(), axis.title.y = element_text(size = 15))
dev.off()

#-----------------------------------------------------------------------------------------------------#

#### Fig.4 ####

## Average MF as a function of the number of LIRRs per stream type

# Dataframe

mf.typo.lirr <- data.frame(mf.mean = tapply(eco.count$BEST, factor(desc$typo_I2M2_V2), mean), # BEST only
                           N.LIRR = N.LIRR[match(as.character(1:66), names(N.LIRR))],
                           N.OPC = N.OPC[match(as.character(1:66), names(N.OPC))])

# Non-linear regression

reg.mf <- nls(mf.mean ~ a*exp(N.LIRR*b) + c, data = mf.typo.lirr, # BEST only
               start = list(a = 1, b = -0.01, c = 0.1))
a = coef(reg.mf)[1]
b = coef(reg.mf)[2]
c = coef(reg.mf)[3]
x <- seq(6,135,0.01)
y <- a*exp(x*b) + c

mf.model <- lm(log(mf.mean) ~ log(a*exp(N.LIRR*b) + c), data = mf.typo.lirr[mf.typo.lirr$mf.mean > 0,]) # With 10k repetitions there should not be any zero, but for testing with less repetitions there can be
stats.mod1 <- summary(mf.model)
mf.lm <- lm(mf.mean ~ N.LIRR, data = mf.typo.lirr) # Linear model
summary(mf.lm)

reg.mf2 <- nls(mf.mean ~ a*exp(N.OPC*b) + c, data = mf.typo.lirr,
               start = list(a = 1, b = -0.1, c = 0.1))
a = coef(reg.mf2)[1]
b = coef(reg.mf2)[2]
c = coef(reg.mf2)[3]
x2 <- seq(1,33,0.01)
y2 <- a*exp(x2*b) + c

mf.model.2 <- lm(log(mf.mean) ~ log(a*exp(N.OPC*b) + c), data = mf.typo.lirr[mf.typo.lirr$mf.mean > 0,]) # With 10k repetitions there should not be any zero, but for testing with less repetitions there can be
stats.mod2 <- summary(mf.model.2)

# Figure

y.range = c(-0.01,0.21)
cex.lab = 1.3
svg(paste0(out.dir,"Fig4_MF per type_",day,".svg"), height = 7, width = 12)
par(mfrow = c(1,2))
# Fig.4A
plot(mf.mean ~ N.LIRR, data = mf.typo.lirr, xlim = c(-1,151), ylim = y.range,
     xlab = "Number of SSEs considered as LIRRs", ylab = "Average MF", cex.lab = cex.lab)
lines(y ~ x, col = "red", lwd = 3, lty = 2)
text(151, range(mf.typo.lirr$mf.mean)[2]*0.95,
     col = "red", cex = 1.3, pos = 2,
     label = paste0("R² = ", round(stats.mod1$r.squared,4),"\np = ", round(stats.mod1$coefficients[2,4],4)))
# Fig.4B
plot(mf.mean ~ N.OPC, data = mf.typo.lirr, ylim = y.range,
     xlab = "Number of sites considered as LIRRs", ylab = "Average MF", cex.lab = cex.lab)
lines(y2 ~ x2, col = "red", lwd = 3, lty = 2)
text(range(mf.typo.lirr$N.OPC)[2], range(mf.typo.lirr$mf.mean)[2]*0.95,
     col = "red", cex = 1.3, pos = 2,
     label = paste0("R² = ", round(stats.mod2$r.squared,4),"\np = ", round(stats.mod2$coefficients[2,4],4)))
dev.off()

# Some basic statistics

range(mf.typo.lirr$N.LIRR)
range(mf.typo.lirr$mf.mean[mf.typo.lirr$N.LIRR < 50])
range(mf.typo.lirr$mf.mean[mf.typo.lirr$N.LIRR >= 50])

#-----------------------------------------------------------------------------------------------------#

#### Fig.5 ####

## MF per ecological class

cee <- factor(cee.opc)
levels(cee) <- nm.ce
summary(cee)
levels(cee) <- paste0(nm.ce,"\n(N = ", summary(cee),")") # Include number of SSEs in each class
# cee <- factor(cee, levels = levels(cee)[5:1])

# Figure

svg(paste0(out.dir,"Fig5_mf per class_",day,".svg"))
par(mfrow = c(1,1))
boxplot(eco.count$OBSWoWgB ~ cee, xlab = "Ecological status", ylab = "Mismatch frequency (MF)", cex.lab = 1.2,
        col = c("red","orange","yellow","green","blue"), range = 0)
dev.off()

# Kruskal-Wallis test

levels(cee) <- nm.ce # Without the number of SSEs in each class
kruskal.test(eco.count$OBSWoWgB ~ cee)
conover.test(eco.count$OBSWoWgB, g = cee, method = "holm") # Post-hoc multiple comparison tests

#-----------------------------------------------------------------------------------------------------#

#### Figures 6 & S1 ####

## MF as a function of EQR values

# Loading results obtained with all uncertainty sources combined together

i <- which(grepl("OBSWoWgB", sel)) # Check the index number of the combination with all uncertainty sources
eco.data <- load(paste0(out.dir, sel[i]))
eco <- get(eco.data)
rm(list = eco.data)

# Adding intercal categories to EQR table

data.eqr <- data.frame(sd.eqr, eqr2)
data.eqr$lim.cat <- typo_I2M2$lim.cat[match(desc$typo_nationale, typo_I2M2$typo_nationale)]

##### Fig.S1 ####

## Based on Kelly et al. 2009 (https://doi.org/10.1007/s10750-009-9872-z)

# Polynomial regression

poly.reg <- nls2(sd.eqr ~ 0.015 + a * eqr2 + b * eqr2^c, data = data.eqr,
                 start = list(a = 0.1, b = -0.1, c = 3))

cf <- coef(poly.reg)
x <- seq(0,1,0.001)
y <- 0.015 + cf[1]*x + cf[2]*x^cf[3]

# Figure

svg(paste0(out.dir,"FigS1_sd per eqr_",day,".svg"))
plot(sd.eqr ~ eqr2, xlab = "I2M2 index (EQR)", ylab = "Standard deviation", cex.lab = 1.2)
lines(y ~ x, col = "red", lwd = 3, lty = 2)
dev.off()

##### Fig.6 ####

# Mismatch Frequency (MF)

zeros <- rep(0, nrow(eco))
eco.count.fig6 <- data.frame(V1 = zeros, V2 = zeros, V3 = zeros, V4 = zeros, V5 = zeros)
for (j in 1:5) {
  eco.count.fig6[,j] <- apply(eco, MAR = 1, function(x) sum(x == j)/ncol(eco))
}

# Data formating

data.fig.eco <- unlist(eco.count.fig6)
data.fig.eco <- data.frame(prob = data.fig.eco, eqr = rep(eqr2, 5), ce = rep(nm.ce, each = 10074), lim.cat = rep(data.eqr$lim.cat, 5))
x <- seq(0,1,0.001)
class.boundaries <- limits[c(2,5,4,3,1),]

# Figure loop

svg(paste0(out.dir,"Fig6_proba per class and SSE_",day,".svg"), width = 12, height = 13.5)
par(mfrow = c(3,2))
for (i in 1:5) {
  
# Calculation of probabilities
# Based on Kelly et al. 2009 (https://doi.org/10.1007/s10750-009-9872-z)
  
  p.bad <- pnorm(x, class.boundaries$lim.PB[i], 0.015 + cf[1]*class.boundaries$lim.PB[i] + cf[2]*class.boundaries$lim.PB[i]^cf[3], lower.tail = F)
  p.poor <- pnorm(x, class.boundaries$lim.MP[i], 0.015 + cf[1]*class.boundaries$lim.MP[i] + cf[2]*class.boundaries$lim.MP[i]^cf[3], lower.tail = F)
  p.moderate <- pnorm(x, class.boundaries$lim.GM[i], 0.015 + cf[1]*class.boundaries$lim.GM[i] + cf[2]*class.boundaries$lim.GM[i]^cf[3], lower.tail = F)
  p.good <- pnorm(x, class.boundaries$lim.HG[i], 0.015 + cf[1]*class.boundaries$lim.HG[i] + cf[2]*class.boundaries$lim.HG[i]^cf[3], lower.tail = F)
  
  pp.bad <- p.bad
  pp.poor <- p.poor - p.bad
  pp.moderate <- p.moderate - p.poor
  pp.good <- p.good - p.moderate
  pp.high <- 1 - p.good
  
  data.fig.eco.kelly <- data.frame(prob = c(pp.bad, pp.poor, pp.moderate, pp.good, pp.high), eqr = rep(x,5), ce = rep(nm.ce, each = length(x)))
  
  # Figure
  
  plot(prob ~ eqr, data = subset(data.fig.eco, lim.cat == paste0("cat",i)), col = c("red","green","blue","yellow","orange")[factor(data.fig.eco$ce[data.fig.eco$lim.cat == paste0("cat",i)])],
       xlab = "I2M2 index (EQR)", ylab = "Confidence of class", cex.lab = 1.5, xlim = c(0,1))
  lines(prob ~ eqr, data = data.fig.eco.kelly, col = "black", lwd = 2)
  abline(v = class.boundaries$lim.PB[i], lwd = 3, lty = 2)
  abline(v = class.boundaries$lim.MP[i], lwd = 3, lty = 2)
  abline(v = class.boundaries$lim.GM[i], lwd = 3, lty = 2)
  abline(v = class.boundaries$lim.HG[i], lwd = 3, lty = 2)
}
dev.off()

summary(data.eqr$lim.cat)

#-----------------------------------------------------------------------------------------------------#


