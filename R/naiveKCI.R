

naiveKCI <- function(x, y, 
                      mykernel = function(x) (1/(1+exp(-27.5512 * (abs(x) - 0.0800)))), 
                      compute.p = c(TRUE, FALSE), 
                      alternative=c("two.sided", "greater", "less"),
                      p.method=c(),
                      alpha=0.05, 
                      interval=c("confidence", "prediction"), 
                      returnAll=c(FALSE, TRUE), 
                     altkmat=0){
  
  #replace with isbool compute.p <- match.arg(compute.p)
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
  if (altkmat){
    kmat <- altkmat
  } else {
    kmat <- mykernel(abs(xdeltamat)) * mykernel(abs(ydeltamat))
  }
  kcimat <- sign(xdeltamat) * sign(ydeltamat) * kmat
  cindex <- (1 + sum(kcimat)/sum(abs(kcimat)))/2

  #### ####
  # Extracting the important vectors from the cimat calculation
  N <- length(x)
  Cvec <- rowSums(abs(kcimat) + kcimat)/2
  Dvec <- -rowSums(kcimat - abs(kcimat))/2
  
  
  # Computing statistical significance using the Pencina approximation:
  C <- sum(Cvec)
  CC <- sum(Cvec * (Cvec-1))
  D <- sum(Dvec)
  DD <- sum(Dvec * (Dvec-1))
  CD <- sum(Cvec*Dvec)
  
  #cindex <- C/(C+D)  ### This should match the above cindex. 
  varp <- 4*((D^2 * CC - 2*C*D*CD + C^2 * DD) / (C + D)^4) # * N * (N-1) / (N-2)
  
  
  if (is.finite(varp) && varp >= 0) {
    sterr <- sqrt(varp / (N-1))
    if (interval == "confidence"){
      ci <- qnorm(p=alpha/2, lower.tail=FALSE) * sterr
      p <- pnorm((cindex - 0.5)/sterr)
    } else{
      p <- pnorm((cindex - 0.5) / sterr)
      ci <- qnorm(p=alpha/2, lower.tail=FALSE) * sterr * sqrt(2)
    }
  } else {
    return(list(cindex=cindex, 
                p.value=1, 
                sterr=NA, 
                lower=NA, 
                upper=NA, 
                relevant.pairs.no=(C+D)/2))
  } 
  
  return(list(cindex=cindex, 
              p.value=switch(alternative, less=p, greater=1-p, two.sided=2*min(p, 1-p)), 
              sterr=sterr, 
              lower=max(cindex - ci, 0), 
              upper=min(cindex + ci, 1), 
              relevant.pairs.no=sum(abs(kcimat))/2))
}
