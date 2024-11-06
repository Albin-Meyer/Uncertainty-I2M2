rm(list=ls())

#-----------------------------------------------------------------------------------------------------#

# This script takes the results from a preliminary inter-operator variability study
# and uses them to simulate 10,000 values for each site sampling event and each
# metric constitutive of the I2M2 index

# This script is also used to generate Table 1 and Figure S2

#-----------------------------------------------------------------------------------------------------#

#### Libraries ####

library(fitdistrplus)

#-----------------------------------------------------------------------------------------------------#

#### Loading data ####

## Output directory

out.dir <- "output/"

## Data

# Inter-operator results

i2m2_2009 <- read.csv2("data/var_2009_2015-06-09/I2M2.csv") # I2M2 values for 2009
eqr_2009 <- read.csv2("data/var_2009_2015-06-09/EQR.csv") # EQR metric values from 2009
metriques_2009 <- read.csv2("data/var_2009_2015-06-09/Metriques.csv") # Raw metric values from 2009
i2m2_2010 <- read.csv2("data/var_2010_2015-06-09/I2M2.csv") # I2M2 values for 2010
eqr_2010 <- read.csv2("data/var_2010_2015-06-09/EQR.csv") # EQR metric values from 2010
metriques_2010 <- read.csv2("data/var_2010_2015-06-09/Metriques.csv") # Raw metric values from 2010

# Metric values for the dataset used for the development of the I2M2

met <- read.csv2("data/metrics_20190315.csv", check.names = F, row.names = 1)

# Information for sampling events

infos_2009 <- read.table("data/opeconts_2009.txt", h=T)
infos_2010 <- read.table("data/opeconts_2010.txt", h=T)

# Day, for file names when exporting

day <- gsub("-","",Sys.Date())

# Stream names

riv <- c("ALLIER","ANCEDUNORD","AUSTREBERTHE","BEUVE","CHER","CLEURIE","COLAGNE","DESGES","DIVES","DOUX","EAULNE","ESCOURCE","GLAND",
         "HERAULT","JOYEUSE","LERGUE","LEYSSE","LOT","OISE","REMEL","ROUVRE","SAULDRE","SAYE","SEGRE")

#### Run  ####

#-----------------------------------------------------------------------------------------------------#

#### Tab.1 ####

# Average of the number of operators per site sampling event

opc09 <- summary(factor(infos_2009$riviere))
summary(opc09)
sd(opc09)

opc10 <- summary(factor(infos_2010$riviere))
summary(opc10)
sd(opc10)

all.opc <- merge(data.frame(opc09), data.frame(opc10), by = 0, all = T)
rownames(all.opc) <- all.opc$Row.names
all.opc <- all.opc[,-1]
all.opc[is.na(all.opc)] <- 0

# write.csv2(all.opc, file = "output/Tab1_all_opc_20191112.csv") # Table 1 of the article

#-----------------------------------------------------------------------------------------------------#

#### Data preparation ####

# Merging results from 2009 and 2010

res_2009 <- metriques_2009
res_2009 <- merge(res_2009, infos_2009, by.x="X", by.y="nb_opecont")

res_2010 <- metriques_2010
res_2010 <- merge(res_2010, infos_2010, by.x="X", by.y="nb_opecont")

res <- rbind(res_2009, res_2010)
res$annee <- as.factor(res$annee) # Year

#-----------------------------------------------------------------------------------------------------#

#### Inter-operator variability ####

# Calculations of inter-operator deviations from the average (for all metrics and the I2M2)

results <- list()
keeps <- c("Shannon..B1B2.","ASPT..B2B3.","Polyvoltinism..B1B2B3.","Ovoviviparity..B1B2B3.","Richness..B1B2B3.")
cv <- matrix(nrow = length(unique(factor(res_2009$no_riviere))) + length(unique(factor(res_2010$no_riviere))) , ncol = length(keeps)) # Coefficients of variation

row.n = 0
for (an in c(2009,2010)) {
  
  # Correction of the stream names vector, due to some difference between both years
  
  if(an == 2009) riv <- c("ALLIER","AUSTREBERTHE","BEUVE","CHER","CLEURIE","COLAGNE","DESGES","DIVES","DOUX","EAULNE","ESCOURCE","GLAND","HERAULT","JOYEUSE","LERGUE","LEYSSE","LOT","OISE","REMEL","ROUVRE","SAULDRE","SAYE","SEGRE")
  if(an == 2010) riv <- c("ALLIER","ANCEDUNORD","AUSTREBERTHE","BEUVE","CHER","CLEURIE","COLAGNE","DESGES","DIVES","DOUX","EAULNE","ESCOURCE","GLAND","HERAULT","JOYEUSE","LERGUE","LEYSSE","LOT","OISE","REMEL","ROUVRE","SAULDRE","SAYE","SEGRE")
  
  for(i in 1:length(riv)) {
    
    # Calculations of the average value
    
    res_riv <- subset(res, subset = (riviere %in% riv[[i]] & annee %in% an))
    tmp <- res_riv[,keeps]
    tmp.mn <- colMeans(tmp)
    
    row.n = row.n + 1
    cv[row.n,] <- apply(tmp, MAR = 2, sd)/tmp.mn # Coefficients of variation, for Figure S2
    
    # Calculations of the deviations
    
    results[[i]] <- data.frame(t(t(tmp) - tmp.mn))
    
    # Finalization step
    
    if (i == length(riv)) {
      out <- do.call(rbind, results)
      assign(paste0("ecarts_", an), out)
      results <- list()
    }
  }
}

# Comparisons of the distributions of deviations observed in 2009 and 2010

for (i in seq(ecarts_2009)) {
  print(colnames(ecarts_2009)[i])
  print(var.test(ecarts_2009[,i], ecarts_2010[,i]))
  print(wilcox.test(ecarts_2009[,i], ecarts_2010[,i])) # Mann-Whitney tests
  boxplot(list(ecarts_2009[,i], ecarts_2010[,i]))
}

# Variances were sometimes significantly different, but distributions were always similar (rank-based tests)

## Merging results obtained both years

ecarts <- rbind(ecarts_2009, ecarts_2010)

# Distributions

sigmet <- c()
pdf(paste0(out.dir,"fit_ecarts_met_",day,".pdf"), width = 10)
for (i in seq(ecarts)) {
  distri <- fitdist(ecarts[,i], "norm", fix.arg = list(mean = 0))
  sigmet[i] <- distri$estimate
  plot(distri)
}
dev.off()

names(sigmet) <- colnames(ecarts)
save(sigmet, file = paste0(out.dir,"sigmet_",day,".RData"))

#-----------------------------------------------------------------------------------------------------#

#### Figure.S2 ####

## Coefficients of variation (CV)

cv <- as.data.frame(cv)

## CV of the values of the index I2M2

i2m2 <- data.frame(I2M2 = c(i2m2_2009$I2M2, i2m2_2010$I2M2), cd.opc = substr(c(i2m2_2009$X, i2m2_2010$X),1,6))
cv$I2M2 <- tapply(i2m2$I2M2, i2m2$cd.opc, sd)/tapply(i2m2$I2M2, i2m2$cd.opc, mean)

# Boxplots

colnames(cv) <- c("Shannon", "ASPT", "Polyvoltinism", "Ovoviviparity", "Richness", "I2M2")
colMeans(cv)
svg(paste0(out.dir,"Fig_S2_boxplot_cv_",day,".svg"), width = 8.5, height = 8.5)
boxplot(cv, col = c("white","white","white","white","white","grey"), ylab = "CV")
dev.off()

#-----------------------------------------------------------------------------------------------------#

