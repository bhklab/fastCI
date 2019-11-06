# count_xties and count_yties are booleans that indicate whether ties should count as 0.5 correctly ordered 
# and 0.5 incorrectly ordered pairs.  if count_*ties = 0, ties are ignored. 

# This function has not yet been tested.

naiveRCI <- function(x, y, 
                     delta_x = 0.2, 
                     delta_y = 0.2, 
                     valid.logic = c("or", "and"),
                     tie.method.x = c("ignore", "half"), 
                     tie.method.y = c("ignore", "half"),
                     compute.p = c(TRUE, FALSE), 
                     alternative = c("two.sided", "greater", "less"),
                     p.method = c(), 
                     alpha=0.05, 
                     interval=c("confidence", "prediction")){
  
  alternative <- match.arg(alternative)
  interval <- match.arg(interval)
  valid.logic <- match.arg(valid.logic)
  tie.method.x = match.arg(tie.method.x)
  tie.method.y = match.arg(tie.method.y)
  
  myCompleteCases <- complete.cases(x,y)
  x <- x[myCompleteCases]
  y <- y[myCompleteCases]
  
  xmat <- matrix(rep(x, length(x)), ncol=length(x))
  ymat <- matrix(rep(y, length(y)), ncol=length(y))
 
  # Deltamat_ij = x_j - x_i 
  xdeltamat <- xmat - t(xmat)
  ydeltamat <- ymat - t(ymat)
  
  # RCI: sign(xdeltamat) * (abs(xdeltamat) > threshold)
  if (valid.logic == "and"){
    rcimat <- sign(xdeltamat) * sign(ydeltamat) * ((abs(xdeltamat) > delta_x) & (abs(ydeltamat) > delta_y))
  } else {
    rcimat <- sign(xdeltamat) * sign(ydeltamat) * ((abs(xdeltamat) > delta_x) | (abs(ydeltamat) > delta_y))
  }
  
  # Extracting the important vectors from the rcimat calculation
  N <- length(x)
  Cvec <- rowSums(rcimat == 1)
  Dvec <- rowSums(rcimat == -1)

  tievec <- numeric(length(Cvec))
  if (tie.method.x == "half" & tie.method.y == "ignore"){
    tievec <- rowSums(xdeltamat * (abs(xdeltamat) >= delta_x) == 0) - (diag(xdeltamat) == 0)
  } else if (tie.method.x == "ignore" & tie.method.y == "half"){
    tievec <- rowSums(ydeltamat * (abs(ydeltamat) >= delta_y) == 0) - (diag(ydeltamat) == 0)
  } else if (tie.method.x == "half" & tie.method.y == "half"){
    tievec <- rowSums(rcimat == 0)  - (diag(rcimat) == 0)
  }
  Cvec <- Cvec + 0.5*tievec
  Dvec <- Dvec + 0.5*tievec

  # Computing statistical significance using the Pencina approximation:
  C <- sum(Cvec)
  CC <- sum(Cvec * (Cvec-1))
  D <- sum(Dvec)
  DD <- sum(Dvec * (Dvec-1))
  CD <- sum(Cvec*Dvec)

  rcindex <- C/(C+D)
  varp <- 4*((D^2 * CC - 2*C*D*CD + C^2 * DD) / (C + D)^4) * N * (N-1) / (N-2)

  if (is.finite(varp) && varp >= 0) {
    sterr <- sqrt(varp / (N-1))
    if (interval == "confidence"){
      ci <- qnorm(p=alpha/2, lower.tail=FALSE) * sterr
      p <- pnorm((rcindex - 0.5)/sterr)
    } else{
      p <- pnorm((rcindex - 0.5) / sterr)
      ci <- qnorm(p=alpha/2, lower.tail=FALSE) * sterr * sqrt(2)
    }
  } else {
    return(list(rcindex=rcindex, 
                p.value=1, 
                sterr=NA, 
                lower=NA, 
                upper=NA, 
                relevant.pairs.no=(C+D)/2))
  } 
  
  return(list(rcindex=rcindex, 
              p.value=switch(alternative, less=p, greater=1-p, two.sided=2*min(p, 1-p)), 
              sterr=sterr, 
              lower=max(rcindex - ci, 0), 
              upper=min(rcindex + ci, 1), 
              relevant.pairs.no=(C+D)/2))
}
