combineCI <- function(CI1, CI2, N1, N2, multvect1, multvect2){
  null1 <- nullCIDist(N1, multvect1, cumulative = FALSE)
  null2 <- nullCIDist(N2, multvect2, cumulative = FALSE)
  
  if(N1 < N2){
    smaller <- null1
    bigger <- null2
  } else {
    smaller <- null2
    bigger <- null1
  }
  # smaller <- pad(smaller, length(bigger))
  res.null <- convolve(bigger, smaller, type = "o")
  
  res.CI <- (choose(N1,2)*CI1 + choose(N2,2)*CI2)/(choose(N1,2) + choose(N2,2))
  res.pairs <- (choose(N1,2)*CI1 + choose(N2,2)*CI2)
  
  res.p <- getCIPvals(cumsum(res.null), res.pairs)
  
  return(c(res.CI, res.p))
}


# combineHelperExact <- function(x1, x2){
#   return(combineCI(x1[1], x2[1], x1[2], x2[2], 1, 1))
# }

# combineHelperContinous 

computeExpectedApproximation <- function(N){
  
  mean <- choose(N, 2)/2 
  
  var <- sqrt((2*N^3 + 3*N^2 - 5*N)/72)
  return(c("mean" = mean, "var"=var))
}




combineCI.2 <- function(x){
  
  x <- x[,complete.cases(t(x))]
  
  total.N <- sum(x[2,])
  if(total.N >= 100){
    vec1 <- x[,1]
    total.out <- choose(vec1[2], 2)
    conc.out <- vec1[1] * total.out
    pars.out <- computeExpectedApproximation(total.out)
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
    nullTable <- makeTableUpToN(max(x[2,]))
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



### useful for N tissue studies!


