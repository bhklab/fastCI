source("~/Code/fastCI/computeNullCI.R")


naive_approximate_normal_from_density <- function(x){
  
  max.x <- which.max(x)
  mean <- mean(which(x == x[max.x]))
  
  var = abs(which.min(abs(x - x[mean]/2)) - mean)/(sqrt(2*log(2)))
  
  return(c("mean" = mean - 1, "var"=var))
  
}

computeExpectedApproximation <- function(N){
  
  mean <- choose(N, 2)/2 
  
  var <- sqrt((2*N^3 + 3*N^2 - 5*N)/72)
  return(c("mean" = mean, "var"=var))
}


plotNormalApproxN <- function(N){
  
  plot.new()
  
  test <- nullCIDist(N,1,cumulative = 0, force_sym = 1)
  
  tt <- naive_approximate_normal_from_density(test)
  
  tt2 <- computeExpectedApproximation(N)
  
  plot(0:(length(test)-1), test, type='l', col='red')
  
  lines(0:(length(test)), dnorm(0:length(test), tt[1], tt[2]), col='blue')
  
  
  lines(0:(length(test)), dnorm(0:length(test), tt2[1], tt2[2]), col='green')
  
  
  legend("topleft", legend = c("True", "Fitted Approximation", "Expected Approximation"),fill = c("red", "blue", "green"))
}

library(entropy)

n.min <- 5
n.max <- 150

ns <- seq(n.min, n.max)

KL.fitted <- numeric(n.max - n.min + 1)
KL.analytic <- numeric(n.max - n.min + 1)

for(i in seq_len(n.max - n.min + 1)){
  
  exact.null <- nullCIDist(ns[i],1,cumulative = 0, force_sym = 1)
  
  pars <- naive_approximate_normal_from_density(exact.null)
  pars2 <- computeExpectedApproximation(ns[i])
  
  normal.fitted <- dnorm(0:(length(exact.null)-1), pars[1], pars[2])
  
  normal.analytic <- dnorm(0:(length(exact.null)-1), pars2[1], pars2[2])
  
  KL.fitted[i] <- KL.empirical(exact.null, normal.fitted)
  KL.analytic[i] <- KL.empirical(exact.null, normal.analytic)
  
}

plot(x = ns, y = KL.fitted, col="blue")
points(x = ns, y = KL.analytic, col = "green")
legend("topright", legend = c("Fitted Approximation", "Expected Approximation"),fill = c("blue", "green"))


plot(x = ns, y = -log10(KL.analytic), col = "green")
points(x = ns, y = -log10(KL.fitted), col="blue")
legend("topleft", legend = c("Fitted Approximation", "Expected Approximation"),fill = c("blue", "green"))


## The analytic formula clearly fits much better
## To check for whether this is a conservative approximation at alpha = 0.05, we will look at differencs in the tails as a function of N


tail.diffs <- numeric(n.max - n.min + 1)
cvs <- numeric(n.max - n.min + 1)

for(i in seq_len(n.max - n.min + 1)){
  
  exact.null <- nullCIDist(ns[i],1,cumulative = 1, force_sym = 1)
  
  critical.val <- which.min(abs(exact.null - 0.025)) # alpha = 0.05 two sided
  cvs[i] <- critical.val
  
  pars <- computeExpectedApproximation(ns[i])
  
  
  # exact.null is 0 based, so we subtract 1 to get actual number of inversions
  tail.diffs[i] <- pnorm(critical.val - 1, mean = pars[1], sd = pars[2]) - exact.null[critical.val]
  
}

plot(ns, tail.diffs)

## Above method underestimates the probability of the tails in this regime consistently, so we will add an inversion into our normal estimate (or alternatively subtract an inversion if using other tail). 

## Testing above for extra inversion added:

tail.diffs <- numeric(n.max - n.min + 1)
cvs <- numeric(n.max - n.min + 1)

for(i in seq_len(n.max - n.min + 1)){
  
  exact.null <- nullCIDist(ns[i],1,cumulative = 1, force_sym = 1)
  
  critical.val <- which.min(abs(exact.null - 0.025)) # alpha = 0.05 two sided
  cvs[i] <- critical.val
  
  pars <- computeExpectedApproximation(ns[i])
  
  
  # exact.null is 0 based, so we subtract 1 to get actual number of inversions
  tail.diffs[i] <- pnorm(critical.val, mean = pars[1], sd = pars[2]) - exact.null[critical.val]
  
}

plot(ns, tail.diffs)


## What is the effect of adding the extra inversion depending on the alpha level? sweeping through all critical values for n = 100

tail.diffs <- numeric(choose(100,2)/2)

null.dist <- nullCIDist(100)
pars <- computeExpectedApproximation(100)

for(i in 1:(choose(100,2)/2)){
  
  tail.diffs[i] <- pnorm(i , mean = pars[1], sd = pars[2]) - null.dist[i]
  
  
}

plot(1:(choose(100,2)/2), tail.diffs)



