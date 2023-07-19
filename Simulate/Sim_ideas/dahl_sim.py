library(BEDMatrix)

simfxn <- function(G, sig2hom, sig2het, tauhet, bin=FALSE, prev=.2, fixmean=TRUE, alpha1=NA, sd.alpha=.1 ){

  S  <- ncol(G)

  epsilon <- sqrt(1-sig2hom) * rnorm(nrow(G))

  y <- as.numeric(
   X %*% alpha +
   G %*% rnorm( S, sd=sqrt( sig2hom/S ) ) +
   epsilon
  )


