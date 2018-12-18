# Make CI null distribution table
# Derives from https://cs.uwaterloo.ca/journals/JIS/VOL4/MARGOLIUS/inversions.pdf

library(polynom)

#' Function makeNullCIDist: Computes the distribution of CI under the assumption that all permutations are equally likely.
#' Parameters: 
#'     n:        positive integer, length of the sequence, i.e. the number of elements
#'     makeplot  boolean, whether to output a plot of the cdf (default=0)
#'     outdir    string, path to output file for the values as .RData
#' Returns either a probability density distribution or a cumulative distribution (if cumulative=1)
#' over the number of discordant pairs k, with a vector of length choose(n,2)+1 corresponding to 0 
#' through choose(n,2) discordant pairs.  It is necessarily symmetric, and assumes that every
#' permutation is equally likely.  The function currently does not handle ties or mCI, and both
#' makeplot and outdir don't actually do anything (fix this). 
makeNullCIDist <- function(n, multvect, makeplot=0, outdir="", cumulative=1){
  # Initialize to zeros; by symmetry, it's sufficient to count up to
  # choose(n,2)/2, but this is for thoroughness.  My dist is the probability
  # distribution on the number of inversions [0, choose(n,2)] on n elements.  
  
  # Example:
  #   mynull100 <- makeNullCIDist(n=100)
  mydist <- rep(0, 2+choose(n,2))
  mydist[1] <- 1
  if (n > 1){
    for (ii in 2:n){
      t <- 1/ii * mydist[1:(choose(ii-1,2)+1)]
      mydist[1:(choose(ii,2)+2)] <- cumsum(c(t, rep(0,ii))) - cumsum(c(rep(0,ii), t))
    }
  }
  
  # Given the multiplicities in v
  multvect <- multvect[multvect > 1]
  multcount <- sum(unlist(lapply(multvect, function(x) choose(x,2))))
  for (ii in multvect){
    mydist <- divide.p(mydist, getSimplePolyProduct(ii,ii))
  }

  if (cumulative == 1){
    cumsum(mydist[1:(choose(n,2)+1)])
  } else {
    mydist[1:(choose(n,2) + 1 - multcount)]
  }
}


getCIPvals <- function(nullCIDist, n, discpairs, alternative=c("two.sided", "greater", "less")){
  if (length(nullCIDist) != choose(n, 2)+1){
    return("Error: length of input nullCIDist does not match n (should be choose(n,2)+1)")
  } else if ((discpairs < 0) | discpairs > choose(n,2)){
    return("Error: number of discordant pairs out of bounds [0, choose(n,2)]")
  } else if (nullCIDist[length(nullCIDist)] < 0.99) {
    return("Error: nullCIDist is not a CDF.  Run makeNullCIDist with argument cumulative=1")
  } else {
    
    alternative <- match.arg(alternative)
    ix_twosided <- min(discpairs, choose(n,2) - discpairs)
    
    switch (alternative, 
            "two.sided" = ifelse(discpairs == choose(n,2)/2, 
                                 1, 
                                 2*nullCIDist[ix_twosided+1]),
            "greater" = 1 - nullCIDist[discpairs+1],
            "less" = nullCIDist[discpairs+1]
    )
  }
}

# library(polynom) does allow operations on polynomials, but it is very slow
# compared to the CDF-trick for convolutions of c(1,1,1,...,1) implemented
# for the tie-less case in makeNullCIDist.
# getMultiplicityPoly computes the D_M(x) = [alpha alpha_{n}}]_{x} polynomial
# Note that alpha = sum_{i =1:n} alpha_{i} - i.e. it includes this additional
# multiplicity.  
# This function turns out to be basically unnecessary.
getMultiplicityPoly <- function(elements, multiplicity){
  num <- getSimplePolyProduct(elements, multiplicity)
  denom <- getSimplePolyProduct(multiplicity, multiplicity)
  
  return(as.numeric(as.polynomial(num)/denom))
}

# This function uses the convolutional trick to multiply k polynomials of the
# form 1/m(1+x+x^2+...x^m) from m=n-k+1 to n. It is normalized so its coefficients sum
# to 1, i.e. it is a probability density function.  This is faster for computation. 
getSimplePolyProduct <- function(elts, range){
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
    }
  }
  return(retpoly[1:count])
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

mult.p <-function(p1, p2, outOrder){
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

divide.p <- function(p1, p2){
  outlen <- length(p1) - length(p2) + 1
  p1 <- pad(p1, length(p1)+1)
  p2 <- pad(p2, length(p1))
  
  p1.fft <- fft(p1)
  p2.fft <- fft(p2)
  p3.fft <- p1.fft / p2.fft
  
  return(fft(1/length(p3.fft) * p3.fft, inverse = TRUE)[1:outlen])
}