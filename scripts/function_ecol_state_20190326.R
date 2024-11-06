
#-------------------------------------------------------------------------------------------------------#

# Day of creation: 22/03/2019

#-------------------------------------------------------------------------------------------------------#

## Function used for the assignment of the ecological class based on the value of the index (I2M2)

# Ecological classes are the following, from best to worst: High, Good, Moderate, Poor, Bad

ecol.state <- function(i2m2, desc, typo_I2M2, ...) {
  
  # Data
  
  I2M2 <- i2m2
  limites <- typo_I2M2[typo_I2M2$typo_nationale == as.character(desc$typo_nationale), c("HG","GM","MP","PB")] # Boundary values (HG High-Good, etc...)
  limites <- rev(unlist(c(Inf,limites,-Inf)))
  
  # Ecological class assignment
  
  # out <- cut(I2M2, limites, include.lowest = T)
  out <- cut(I2M2, limites, right = F)
  levels(out) <- 1:5
  out <- as.numeric(out)
  
  # Export
  
  out
  
}

#-------------------------------------------------------------------------------------------------------#
