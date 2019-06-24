# count_xties and count_yties are booleans that indicate whether ties should count as 0.5 correctly ordered 
# and 0.5 incorrectly ordered pairs.  if count_*ties = 0, ties are ignored. 

naiveCI <- function(x, y, count_xties=0, count_yties=0, alternative=c("two.sided", "greater", "less"),
                     alpha=0.05, interval=c("confidence", "prediction")){
  alternative <- match.arg(alternative)
  interval <- match.arg(interval)
  
  myCompleteCases <- complete.cases(x,y)
  x <- x[myCompleteCases]
  y <- y[myCompleteCases]
  
  xmat <- matrix(rep(x, length(x)), ncol=length(x))
  ymat <- matrix(rep(y, length(y)), ncol=length(y))
 
  # Deltamat_ij = x_j - x_i 
  xdeltamat <- xmat - t(xmat)
  ydeltamat <- ymat - t(ymat)
  
  # RCI: sign(xdeltamat) * (abs(xdeltamat) > threshold)
  cimat <- sign(xdeltamat) * sign(ydeltamat)
  
  # Extracting the important vectors from the cimat calculation
  N <- length(x)
  Cvec <- rowSums(cimat == 1)
  Dvec <- rowSums(cimat == -1)

  tievec <- numeric(length(Cvec))
  # Subtract 1 for the diagonal, for which delta trivially equals 0. 
  if (count_xties == 1 & count_yties == 0){
    tievec <- rowSums(xdeltamat == 0) - 1
  } else if (count_xties == 0 & count_yties == 1){
    tievec <- rowSums(ydeltamat == 0) - 1
  } else if (count_xties == 1 & count_yties == 1){
    tievec <- rowSums(cimat == 0) - 1
  }
  Cvec <- Cvec + 0.5*tievec
  Dvec <- Dvec + 0.5*tievec

  # Computing statistical significance using the Pencina approximation:
  C <- sum(Cvec)
  CC <- sum(Cvec * (Cvec-1))
  D <- sum(Dvec)
  DD <- sum(Dvec * (Dvec-1))
  CD <- sum(Cvec*Dvec)

  cindex <- C/(C+D)
  varp <- 4*((D^2 * CC - 2*C*D*CD + C^2 * DD) / (C + D)^4) * N * (N-1) / (N-2)
  #browser()
  
  if (varp >= 0) {
    sterr <- sqrt(varp / (N-1))
    if (interval == "confidence"){
      ci <- qnorm(p=alpha/2, lower.tail=FALSE) * sterr
      p <- pnorm((cindex - 0.5)/sterr)
    } else{
      p <- pnorm((cindex - 0.5) / sterr)
      ci <- qnorm(p=alpha/2, lower.tail=FALSE) * sterr * sqrt(2)
    }
  } else {
    return(list(cindex=cindex, p.value=1, sterr=NA, lower=0, upper=0, 
           relevant.pairs.no=(C+D)/2, 
           cimat=cimat))
  } 
  
  return(list(cindex=cindex, 
              p.value=switch(alternative, less=p, greater=1-p, two.sided=2*min(p, 1-p)), 
              sterr=sterr, 
              lower=max(cindex - ci, 0), 
              upper=min(cindex + ci, 1), 
              relevant.pairs.no=(C+D)/2, 
              cimat=cimat))
}
