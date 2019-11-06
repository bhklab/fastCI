# count_xties and count_yties are booleans that indicate whether ties should count as 0.5 correctly ordered 
# and 0.5 incorrectly ordered pairs.  if count_*ties = 0, ties are ignored. 

naiveCI <- function(x, y, 
                    tie.method.x=c("ignore", "half"), 
                    tie.method.y=c("ignore", "half"), 
                    compute.p = c(TRUE, FALSE), 
                    alternative=c("two.sided", "greater", "less"),
                    p.method=c(),
                    alpha=0.05, 
                    interval=c("confidence", "prediction"), 
                    returnAll=c(FALSE, TRUE)){
  
  tie.method.x <- match.arg(tie.method.x)
  tie.method.y <- match.arg(tie.method.y)
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
  cimat <- sign(xdeltamat) * sign(ydeltamat)
  
  # Extracting the important vectors from the cimat calculation
  N <- length(x)
  Cvec <- rowSums(cimat == 1)
  Dvec <- rowSums(cimat == -1)

  tievec <- numeric(length(Cvec))
  # Subtract 1 for the diagonal, for which delta trivially equals 0. 
  if (tie.method.x == "half" & tie.method.y == "ignore"){
    tievec <- rowSums(xdeltamat == 0) - 1
  } else if (tie.method.x == "ignore" & tie.method.y == "half"){
    tievec <- rowSums(ydeltamat == 0) - 1
  } else if (tie.method.x == "half" & tie.method.y == "half"){
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
              relevant.pairs.no=(C+D)/2))
}
