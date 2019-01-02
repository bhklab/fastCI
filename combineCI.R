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
  
  browser()
  res.p <- getCIPvals(cumsum(res.null), res.pairs)
  
  return(c(res.CI, res.p))
}


### useful for N tissue studies!


