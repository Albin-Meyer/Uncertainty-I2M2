
#-------------------------------------------------------------------------------------------------------#

# Date of creation: 21/03/2019

#-------------------------------------------------------------------------------------------------------#

## Simplified function for calculations of EQR values

quick.eqr <- function(met, desc, nm.met, real.worst, real.best, ...) {
  
  # Empty output table
  
  eqr <- data.frame(matrix(ncol = 5, nrow = nrow(met)))
  
  # Calculations
  
  for (i in seq(eqr)) {
    eqr[,i] <- (met[,i] - real.worst[i]) / (real.best[as.character(desc$typo_I2M2_V2),i] - real.worst[i])
  }
  
  # Correction of values < 0 or > 1
  
  eqr[eqr < 0] <- 0
  eqr[eqr > 1] <- 1
  
  # Export
  
  colnames(eqr) <- nm.met
  eqr
  
}

#-------------------------------------------------------------------------------------------------------#
