
#-------------------------------------------------------------------------------------------------------#

# Function used for the calculations of discrimination efficiency values of multiple metrics
# Both DE0.25 (using the first quartile of LIRRs) and DE0.75 (using the third quartile of LIRRs)

calcul.de <- function(metrics, RR, type = 7, robust = NULL, ...) {
  
  N.IRR <- length(RR[which(RR == "IRR")]) # Number of IRRs
  Q.LIRR = DE0.25 = DE0.75 = NULL # Empty vectors
  if (ncol(metrics) == 1) {
    Q.LIRR <- quantile(metrics[which(RR == "LIRR"),1], type = type)
    DE0.25 <- length(which(metrics[which(RR == "IRR"),1] < Q.LIRR[2]))/N.IRR
    DE0.75 <- length(which(metrics[which(RR == "IRR"),1] > Q.LIRR[4]))/N.IRR
  } else {
    for (i in 1:ncol(metrics)) {
      Q.LIRR <- quantile(metrics[which(RR == "LIRR"),i], type = type)
      DE0.25[i] <- length(which(metrics[which(RR == "IRR"),i] < Q.LIRR[2]))/N.IRR
      DE0.75[i] <- length(which(metrics[which(RR == "IRR"),i] > Q.LIRR[4]))/N.IRR
    }
  }
  DE <- as.data.frame(cbind(DE0.25,DE0.75))
  rownames(DE) <- colnames(metrics)
  DE
  
}

#-------------------------------------------------------------------------------------------------------#
