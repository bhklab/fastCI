
pad <- function(x, len, location = c("end", "start", "sym")){
  
  location = match.arg(location)
  
  if(length(x) > len){
    stop("Cannot do negative padding.")
  }
  
  if(location == "end"){
    return(c(x, numeric(len - length(x))))
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



testNullCI <- function(n){
  
  p <- 1
  
  for (i in seq(2, n)){
    p <- 1/i*mult.p(p, rep(1,i), choose(i,2)+1)
  }
  
  return(Re(p))
}


testNullCI.poly <- function(n){
  
  p <- polynomial(1)
  
  for (i in seq(2, n)){
    p <- 1/i* p * rep(1,i)
  }
  
  return(Re(p))
}
