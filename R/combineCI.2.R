## Combines CIs and p-values using approximation in case that it is necessary

combineCI.2 <-
function(x, nullTable){
  
  x <- x[,complete.cases(t(x))]
  
  total.N <- sum(x[2,])
  if(total.N >= 100){
    vec1 <- x[,1]
    total.out <- choose(vec1[2], 2)
    conc.out <- vec1[1] * total.out
    pars.out <- computeExpectedApproximation(vec1[2])
    x <- x[,-1, drop=FALSE]
    while(NCOL(x)){
      vec1 <- x[,1]
      pars.new <- computeExpectedApproximation(vec1[2])
      pars.out[1] <- pars.out[1] + pars.new[1]
      pars.out[2] <- sqrt(pars.out[2]^2 + pars.new[2]^2)
      conc.out <- conc.out + vec1[1]*choose(vec1[2], 2)
      total.out <- total.out + choose(vec1[2], 2)
      x <- x[,-1, drop=FALSE]
    }
    
    if(conc.out < pars.out[1]){
      prob <- pnorm(conc.out + 1,mean = pars.out[1], sd = pars.out[2])
      p.out <- 2*prob
    } else {
      prob <- pnorm(conc.out - 1,mean = pars.out[1], sd = pars.out[2], lower.tail = FALSE)
      p.out <- 2*prob
    }
    CI.out <- conc.out/total.out
    
  } else {
    if(missing(nullTable)){
      nullTable <- makeTableUpToN(max(x[2,]))
    }
    vec1 <- x[,1]
    x <- x[,-1, drop=FALSE]
    nullOut <- nullTable[[vec1[2]]]
    conc.out <- vec1[1]*choose(vec1[2], 2)
    total.out <- choose(vec1[2], 2)
    while(NCOL(x)){
      vec1 <- x[,1]
      nullOut <- convolve(nullOut, nullTable[[vec1[2]]], type='o')
      total.out <- total.out + choose(vec1[2],2)
      conc.out <- conc.out + vec1[1]*choose(vec1[2], 2)
      x <- x[,-1, drop=FALSE]
    }
    CI.out <- conc.out/total.out
    p.out <- getCIPvals(cumsum(nullOut), total.out - conc.out)
  }
  
  return(c("CI" = CI.out[[1]], "p" = p.out[[1]]))
  
}
