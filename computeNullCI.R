# Make CI null distribution table
# Derives from https://cs.uwaterloo.ca/journals/JIS/VOL4/MARGOLIUS/inversions.pdf

library(polynom)
library(purrr)

#' Function nullCIDist: Computes the distribution of CI under the assumption that all permutations are equally likely.
#' Parameters: 
#'     n:        positive integer, length of the sequence, i.e. the number of elements
#'     makeplot  boolean, whether to output a plot of the cdf (default=0)
#'     outdir    string, path to output file for the values as .RData
#' Returns either a probability density distribution or a cumulative distribution (if cumulative=1)
#' over the number of discordant pairs k, with a vector of length choose(n,2)+1 corresponding to 0 
#' through choose(n,2) discordant pairs.  It is necessarily symmetric, and assumes that every
#' permutation is equally likely.  The function currently does not handle ties or mCI, and both
#' makeplot and outdir don't actually do anything (fix this). 
nullCIDist <- function(n, multvect=c(1), makeplot=0, outdir="", cumulative=1, return_polylist=0, force_sym=0){
  # Initialize to zeros; by symmetry, it's sufficient to count up to
  # choose(n,2)/2, but this is for thoroughness.  My dist is the probability
  # distribution on the number of inversions [0, choose(n,2)] on n elements.  
  
  ## Example:
  ##   mynull100 <- nullCIDist(n=100)
  if (sum(multvect) < n) {
    multvect <- c(multvect, rep(1, n - sum(multvect)))
  }
  
  # Given the multiplicities in v
  if (sum(multvect == 1) > 0){
    polylist <- list(getSimplePolyProduct(elts=sum(multvect == 1), range=sum(multvect == 1), norm=1))
  } else {
    polylist <- list(c(1))
  }

  multcount <- sum(multvect == 1)
  mvect <- multvect[multvect > 1]
  for (ii in mvect){
    polylist <- c(polylist, getSmartMultiplicity(multcount + ii, ii, norm=1))
    multcount <- multcount + ii
  }
  
  #mydist <- as.numeric(reduce(polylist, mult.p))
  mydist <- as.numeric(mult.plist(polylist))
  # There is a symmetry problem with mydist
  if (force_sym == 1){
    mydist[length(mydist):(length(mydist)/2 + 1)] = mydist[1:(length(mydist)/2)]
  }
  
  if (cumulative == 1){
    mydist <- cumsum(mydist)
  }
  
  if (return_polylist == 1){
    return(polylist) }  #for debugging
  else{
    return(as.numeric(mydist)) }
}


getCIPvals <- function(nullCIDist, discpairs, alternative=c("two.sided", "greater", "less")){
  if (nullCIDist[length(nullCIDist)] < 0.99) {
    return("Error: nullCIDist is not a CDF.  Run nullCIDist with argument cumulative=1")
  } else {
    alternative <- match.arg(alternative)
    ix_twosided <- min(discpairs, length(nullCIDist) - discpairs)
    
    switch (alternative, 
            "two.sided" = ifelse(discpairs == length(nullCIDist)/2, 
                                 1, 
                                 2*nullCIDist[ix_twosided+1]),
            "greater" = 1 - nullCIDist[discpairs+1],
            "less" = nullCIDist[discpairs+1]
    )
  }
}


# library(polynom) does allow operations on polynomials, but it is very slow
# compared to the CDF-trick for convolutions of c(1,1,1,...,1) implemented
# for the tie-less case in nullCIDist.
# getMultiplicityPoly computes the D_M(x) = [alpha alpha_{n}}]_{x} polynomial
# Note that alpha = sum_{i =1:n} alpha_{i} - i.e. it includes this additional
# multiplicity.  
# This function turns out to be basically unnecessary.
# Haha just kidding, 
getMultiplicityPoly <- function(elements, multiplicity){
  num <- getSimplePolyProduct(elements, multiplicity)
  denom <- getSimplePolyProduct(multiplicity, multiplicity)
  
  return(as.numeric(as.polynomial(num)/denom))
}


getSmartMultiplicity <- function(elements, multiplicity, norm=1) {
  if (multiplicity == 1) {
    return(getSimplePolyProduct(elements, elements, norm=norm))
  }
  
  numlist <- getSimplePolyList(elements, multiplicity, norm=norm)
  for (ii in multiplicity:2){
    x <- rep(1/ii, ii)
    # This is very slow because of the polynum modulo operation and is very ghetto
    # The operation doesn't check that there exists an element jj for which x is a factor.  
    # Fix this.  Moving the new element to the end 
    for (jj in 1:length(numlist)){
      if (sum(abs(as.numeric(as.polynomial(numlist[[jj]]) %% x))) < 10^-14){
        numlist <- c(numlist, list(as.numeric(divide.p(numlist[[jj]], x))))
        numlist[jj] <- c()
        break
      }
    }
  }

  if (sum(unlist(lapply(numlist, function(x) length(x) - 1))) != (elements - multiplicity) * multiplicity){
    
  }
  return(numlist) #(reduce(numlist,mult.p))
}


# This function uses the convolutional trick to multiply k polynomials of the
# form 1/m(1+x+x^2+...x^m) from m=n-k+1 to n. It is normalized so its coefficients sum
# to 1, i.e. it is a probability density function.  This is faster for computation. 
getSimplePolyProduct <- function(elts, range, norm=0){
  if (range > elts){
    return("Error: Range (or number of polynomials) is greater than the input number of elements")}
  if (range < 1 | range != round(range)){
    return("Error: Range must be a positive integer")}

  #initialize - initial length is sum of (elts-1):(elts-range) + 1
  #This initialization is a bit gross, but it's faster to start at the bottom of the range
  retpoly <- c(rep(1, (elts-range+1)), rep(0, range/2 * (2*elts-range-1) + 1 - (elts-range+1))) #* 1/(elts-range+1)
  count <- (elts - range + 1)
  if (range > 1){
    for (k in (elts-range+2):elts){
      retpoly[1:(count+k)] = (cumsum(c(retpoly[1:count], rep(0,k))) - cumsum(c(rep(0,k), retpoly[1:count])))# * 1/k
      count <- count + k-1
      # Normalize stepwise, floor out guys at 10^-200?
    }
  }
  if (norm != 0){
    retpoly <- retpoly / sum(retpoly)
  }
  return(retpoly[1:count])
}


getSimplePolyList <- function(elts, range, norm=0){
  if (range > elts) {
    return("Error: Range (or number of polynomials) is greater than the input number of elements")}
  if (range < 1 | range != round(range)) {
    return("Error: Range must be a positive integer")}
  a <- list()

  for (ii in (elts - range + 1):elts){
    k <- ifelse(norm != 0, 1/ii, 1)
    a <- c(a, list(rep(k,ii)))
  }
  return(a)
}


pad <- function(x, len, location = c("end", "start", "sym")){
  location = match.arg(location)
  if(length(x) > len){
    stop("Cannot do negative padding.")
  }
  if (location == "end"){
    return(c(x, numeric(len - length(x))))
  }
  if (location == "start"){
    return(c(numeric(len - length(x)), x))
  }
}


mult.p <- function(p1, p2, outOrder){
  if(missing(outOrder)){
    outOrder <- length(p1) + length(p2) - 1
  } 
  p1 <- pad(p1, outOrder)
  p2 <- pad(p2, outOrder)
  
  p1.fft <- fft(p1)
  p2.fft <- fft(p2)
  p3.fft <- p1.fft*p2.fft
  
  return(fft(1/length(p3.fft) * p3.fft, inverse = TRUE))
}


mult.plist <- function(plist, outOrder){
  # Fix the float overflow problem here with SmartMultipying
  if (missing(outOrder)) {
    outOrder <- sum(unlist(lapply(plist, length))) - length(plist) + 1
  }
  if (length(plist) == 1){
    return(plist[[1]])
  }
  
  product <- fft(pad(plist[[1]], outOrder))
  for (ii in 2:length(plist)){
    product <- product * fft(pad(plist[[ii]], outOrder))
  }
  
  return(fft(1/length(product) * product, inverse = TRUE))
}


divide.p <- function(p1, p2){
  outlen <- length(p1) - length(p2) + 1
  p1 <- pad(p1, length(p1)+1)
  p2 <- pad(p2, length(p1))
  
  p1.fft <- fft(p1)
  p2.fft <- fft(p2)
  p3.fft <- p1.fft / p2.fft
  
  return(fft(1/length(p3.fft) * p3.fft, inverse = TRUE)[1:outlen])
}


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