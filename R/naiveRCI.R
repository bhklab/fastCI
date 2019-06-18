# count_xties and count_yties are booleans that indicate whether ties should count as 0.5 correctly ordered 
# and 0.5 incorrectly ordered pairs.  if count_*ties = 0, ties are ignored. 

# This function has not yet been tested.

naiveRCI <- function(x, y, delta_x = 0.2, delta_y = 0.2, count_xties=0, count_yties=0, logic.operator=c("and", "or"), 
                     alternative=c("two.sided", "greater", "less"),
                     alpha=0.05, interval=c("confidence", "prediction")){
  alternative <- match.arg(alternative)
  interval <- match.arg(interval)
  logic.operator <- match.arg(logic.operator)
  
  myCompleteCases <- complete.cases(x,y)
  x <- x[myCompleteCases]
  y <- y[myCompleteCases]
  
  xmat <- matrix(rep(x, length(x)), ncol=length(x))
  ymat <- matrix(rep(y, length(y)), ncol=length(y))
 
  # Deltamat_ij = x_j - x_i 
  xdeltamat <- xmat - t(xmat)
  ydeltamat <- ymat - t(ymat)
  
  # RCI: sign(xdeltamat) * (abs(xdeltamat) > threshold)
  if (logic.operator == "and"){
    rcimat <- sign(xdeltamat) * sign(ydeltamat) * ((abs(xdeltamat) > delta_x) & (abs(ydeltamat) > delta_y))
  } else {
    rcimat <- sign(xdeltamat) * sign(ydeltamat) * ((abs(xdeltamat) > delta_x) | (abs(ydeltamat) > delta_y))
  }
  
  # Extracting the important vectors from the rcimat calculation
  N <- length(x)
  Cvec <- rowSums(rcimat == 1)
  Dvec <- rowSums(rcimat == -1)

  tievec <- numeric(length(Cvec))
  if (count_xties == 1 & count_yties == 0){
    tievec <- rowSums(xdeltamat * (abs(xdeltamat) >= delta_x) == 0) - (diag(xdeltamat) == 0)
  } else if (count_xties == 0 & count_yties == 1){
    tievec <- rowSums(ydeltamat * (abs(ydeltamat) >= delta_y) == 0) - (diag(ydeltamat) == 0)
  } else if (count_xties == 1 & count_yties == 1){
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

  if (varp >= 0) {
    sterr <- sqrt(varp / (N-1))
    if (interval == "confidence"){
      ci <- qnorm(p=alpha/2, lower.tail=FALSE) * sterr
      p <- pnorm((rcindex - 0.5)/sterr)
    } else{
      p <- pnorm((rcindex - 0.5) / sterr)
      ci <- qnorm(p=alpha/2, lower.tail=FALSE) * sterr * sqrt(2)
    }
  } else {
    return(list(rcindex=rcindex, p.value=1, sterr=NA, lower=0, upper=0, 
           relevant.pairs.no=(C+D)/2))
  } 
  
  return(list(rcindex=rcindex, 
              p.value=switch(alternative, less=p, greater=1-p, two.sided=2*min(p, 1-p)), 
              sterr=sterr, 
              lower=max(rcindex - ci, 0), 
              upper=min(rcindex + ci, 1), 
              relevant.pairs.no=(C+D)/2, 
              rcimat = rcimat))  # delete this line
}
