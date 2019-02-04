# Petr is not allowed in this file
# 
# This script includes various simulations and trials to test the properties of CI.
# Author: Ian Smith
# 31 December 2018

library(ggplot2)
source("computeNullCI.R")
source("fastCI.R")

outdir <- "/Users/iansmith/Work/bhk/analysis/ci_analysis"

# Compare binomial distribution with vanilla CI for intuition
nlist = 50 #seq(10, 150, 10)
cisd <- numeric(length(nlist))
binomsd <- numeric(length(nlist))

for (ii in 1:length(nlist)){
  n <- nlist[ii]
  cinull <- nullCIDist(n, cumulative=0)
  binomnull <- dbinom(x=0:choose(n,2), size=choose(n,2), prob=0.5)
  plot(x=seq(0,1, 1/choose(n,2)), y=cinull, type="l", col="red", lwd=2, ylim=c(0, 1.1*max(binomnull)), xlab="Number of Inversions", ylab="Density", 
     main=paste("Probability Distribution on Inversions for Vanilla CI vs Binomial with C(N,2) independent trials, N =", n))
  lines(x=seq(0,1, 1/choose(n,2)), y=binomnull, type="l", col="blue", lwd=2)
  grid()
  legend(x="topright", c("Vanilla CI", "Binomial"), col=c("red","blue"), lwd=c(2,2))
  if (is.element(n, c(20, 50, 100))){
    dev.print(pdf, paste(outdir, paste("binom_nullcomp_dist_vs_ci_n=", n, ".pdf", sep=""), sep="/"))
  }

  cisd[ii] <- sqrt(sum(cinull * ((seq(0,1,1/choose(n,2)) - 0.5)^2)))
  binomsd[ii] <- 0.5/sqrt(choose(n,2))
}

plot(x=nlist, cisd, col="red", type="b", lwd=2, ylim=c(0, 1.1*max(cisd)), 
     xlab="Length of vectors", ylab="Standard Deviation", main="Distribution of Standard Dev for CI null compared to a binomial model vs vector length n")
lines(x=nlist, binomsd, col="blue", type="b", lwd=2)
legend(x="topright", c("Vanilla CI", "Binomial"), col=c("red","blue"), lwd=c(2,2))
grid()
dev.print(pdf, paste(outdir, "binom_nullcomp_sd_vs_ci_x=n.pdf", sep="/"), width=10, height=8)


# Validate Vanilla CI correctly computes the null distribution
n = 50
mynull <- nullCIDist(n=n)
pdfnull <- nullCIDist(n=n, cumulative=0)
perms <- lapply(1:1e5, function(x) fastCI(c(1:n), sample(1:n,n,replace=FALSE))$cindex)
mysample <- cinull.sample.reject(pdfnull, n=1e4)
discsample <- cinull.sample(mynull, n=1e4)

permhist <- lapply(seq(0,1,1/(length(mynull)-1)), function(x) sum(unlist(perms) < x))
plot(seq(0,1,1/(length(mynull)-1)), unlist(permhist), type="l", col="blue", lwd=2, xlim=c(0.3, 0.7),
     main=paste("CDF on Vanilla CI for n=", n), xlab="CI", ylab="Permutation count")
lines(seq(0,1,1/(length(mynull)-1)), length(perms)*mynull, type="l", col="red", lwd=2)
legend(x="bottomright", c("Permutation null", "Analytical null"), col=c("blue", "red"), lwd=c(2,2))
grid()
dev.print(pdf, paste(outdir, paste("nullcomp_perm_vs_analyticalnull_n=", n, ".pdf", sep=""), sep="/"), width=10, height=8)

plot(seq(0,1,1/(length(mynull)-1)), (unlist(permhist) - length(perms)*mynull)/length(perms))


myp <- numeric(1000)
discp <- numeric(1000)
for (ii in 1:1000){
  mysample <- cinull.sample.reject(pdfnull, n=1e4)
  discsample <- cinull.sample(mynull, n=1e4)
  myp[ii] <- ks.test(unlist(perms), mysample)$p.value
  discp[ii] <- ks.test(unlist(perms), discsample)$p.value
  if (ii %% 100 == 0){print(ii)}
}
hist(myp, breaks=seq(0,1,0.05), col=rgb(1,0,0,0.5), main="Distribution of KS-test p-values for permutation null vs sampling methods")
hist(discp, breaks=seq(0,1,0.05), col=rgb(0,0,1,0.5), add=T)
legend(x="topright", c("Continuous Rejection Sampling", "Discrete sampling"), col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), lty=c(1,1), lwd=c(10,10))
grid()
dev.print(pdf, paste(outdir, "nullcomp_perm_vs_sampling_kstest_n=50.pdf", sep="/"), width=10, height=8)

ggplot2()



# Problems: asymmetries in the calculation, especially on the right side of the distribution:
# Floating point errors?
mynull <- nullCIDist(n=10, cumulative=0)
plot(mynull[8000:8200])
points(mynull[3175:2976], col="blue")
# Because the distribution is symmetric, this should be a flat line at y = 0:
mynull <- nullCIDist(n=20, cumulative=0)
plot(mynull[1:(length(mynull)/2)] - mynull[length(mynull):(length(mynull)/2+1)])
plot((mynull[1:(length(mynull)/2)] - mynull[length(mynull):(length(mynull)/2+1)])/mynull[1:(length(mynull)/2)])

