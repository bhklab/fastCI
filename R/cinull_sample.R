cinull.sample <- function(mycdf, n, range=c(0,1)){
  # Use for small values of mycdf
  if (length(mycdf) > 20000){
    stop("Use cinull.sample.reject instead of cinull.sample")
  }

  myprobs <- runif(n)
  mysample <- unlist(lapply(myprobs, function(x) min(which(mycdf > x))))
  mysample <- (mysample-1)/(length(mycdf) - 1)
  return(mysample)
}
