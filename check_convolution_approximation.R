source("~/Code/fastCI/computeNullCI.R")

n1 <- 10
n2 <- 15


naive_approximate_normal_from_density <- function(x){
  
  mean <- sum((seq_along(x) - 1)*x)
  
  sd = sqrt(sum((seq_along(x) - 1)^2*x) - sum((seq_along(x) - 1)*x)^2)
  
  return(c("mean" = mean, "sd"=sd))
  
}

computeExpectedApproximation <- function(N){
  
  mean <- choose(N, 2)/2 
  
  var <- sqrt((2*N^3 + 3*N^2 - 5*N)/72)
  return(c("mean" = mean, "var"=var))
}



plotConvolutionApprox <- function(n1, n2){
  
  plot.new()
  
  test1 <- nullCIDist(n1,1,cumulative = 0, force_sym = 1)
  test2 <- nullCIDist(n2,1, cumulative = 0, force_sym = 1)
  tt1 <- computeExpectedApproximation(n1)
  tt2 <- computeExpectedApproximation(n2)
  
  test <- convolve(test1, test2, type='o')
  
  tt <- tt1 + tt2
  tt[2] <- sqrt(tt1[2]^2 + tt2[2]^2)
  
  # lines(0:(length(test)-1), test, type='l', col='red')
  
  plot(0:(length(test)-1), dnorm(0:(length(test)-1), tt[1], tt[2]), col='blue', type='l')
  lines(0:(length(test)-1), test, type='l', col='red')
  
  
  legend("topleft", legend = c("True", "Normal Approx Convolution"), fill = c("red", "blue"))
}



plotConvolutionApprox3 <- function(n1, n2, n3){
  
  plot.new()
  
  test1 <- nullCIDist(n1,1,cumulative = 0, force_sym = 1)
  test2 <- nullCIDist(n2,1, cumulative = 0, force_sym = 1)
  test3 <- nullCIDist(n3,1, cumulative = 0, force_sym = 1)
  
  tt1 <- computeExpectedApproximation(n1)
  tt2 <- computeExpectedApproximation(n2)
  tt3 <- computeExpectedApproximation(n3)
  
  test <- convolve(convolve(test1, test2, type='o'), test2, type='o')
  
  tt <- tt1 + tt2
  tt[2] <- sqrt(tt1[2]^2 + tt2[2]^2)
  
  tt[1] <- tt[1] + tt3[1]
  tt[2] <- sqrt(tt[2]^2 + tt3[2]^2)
  
  
  # lines(0:(length(test)-1), test, type='l', col='red')
  
  plot(0:(length(test)-1), dnorm(0:(length(test)-1), tt[1], tt[2]), col='blue', type='l')
  lines(0:(length(test)-1), test, type='l', col='red')
  
  
  legend("topleft", legend = c("True", "Normal Approx Convolution"), fill = c("red", "blue"))
}


#### going to investigate whether KL divergence grows as we add more combine steps, compared to refitting a normal approximation at each step.


each.n <- 10

num.convs <- 10

library(entropy)

KL.analytic <- numeric(num.convs)
KL.refit <- numeric(num.convs)



test1 <- nullCIDist(each.n, 1, cumulative = 0, force_sym = 1)
tt1 <- computeExpectedApproximation(each.n)
tt2 <- computeExpectedApproximation(each.n)

for(i in seq_len(num.convs)){
  
  test1 <- convolve(test1, nullCIDist(each.n, 1, cumulative = 0, force_sym = 1), type = 'o')
  tt1[1] <- tt1[1] + tt2[1]
  tt1[2] <- sqrt(tt1[2]^2 + tt2[2]^2)
  
  
  dist <- dnorm(seq_along(test1)-1, tt1[1], tt1[2])
  
  pp <- naive_approximate_normal_from_density(test1)
  pp.dist <- dnorm(seq_along(test1)-1, pp[1], pp[2])
  
  KL.analytic[i] <- KL.plugin(test1, dist)
  KL.refit[i] <- KL.plugin(test1, pp.dist)
}

plot(KL.analytic, col="blue")
points(KL.refit, col="red")


plot(KL.refit - KL.analytic)


