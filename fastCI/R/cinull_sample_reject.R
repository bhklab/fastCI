# Rejection sampling for comparison with permutation nulls:
cinull.sample.reject <- function(mypdf, range=c(0,1), n){
  # mypdf is a probability distribution function sampled uniformly on some interval (i.e. a histogram)
  mymean <- sum(mypdf * seq(range[1], range[2], 1/(length(mypdf)-1)))
  mysd <- sqrt(sum(mypdf * ((seq(range[1], range[2], 1/(length(mypdf)-1)) - 0.5)^2)))
  myratio <- max(mypdf/dnorm(x=seq(0,1,1/(length(mypdf)-1)), mean=mymean, sd=1.2*mysd))

  mysamples <- numeric(n)
  ii = 1
  while(ii <= n){
    x <- rnorm(1, mean=mymean, sd=mysd*1.2)
    u <- runif(1)

    if (x > range[1] & x < range[2] & u <= mypdf[round(x * (length(mypdf)-1))+1] / (dnorm(x, mean=mymean, sd=1.2*mysd) * myratio)){
      mysamples[ii] <- x
      ii <- ii + 1
    }
  }
  return(mysamples)
}
