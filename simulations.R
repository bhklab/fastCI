# Petr is not allowed in this file
# 
# This script includes various simulations and trials to test the properties of CI.
# Author: Ian Smith
# 31 December 2018

outdir <- "/Users/iansmith/Work/bhk/analysis/ci_analysis"

# Compare binomial distribution with vanilla CI for intuition
nlist = seq(10, 150, 10)
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



# Problems: asymmetries in the calculation, especially on the right side of the distribution:
# Floating point errors?
mynull <- nullCIDist(n=10, cumulative=0)
plot(mynull[8000:8200])
points(mynull[3175:2976], col="blue")
# Because the distribution is symmetric, this should be a flat line at y = 0:
mynull <- nullCIDist(n=150, cumulative=0)
plot(mynull[1:(length(mynull)/2)] - mynull[length(mynull):(length(mynull)/2+1)])
plot((mynull[1:(length(mynull)/2)] - mynull[length(mynull):(length(mynull)/2+1)])/mynull[1:(length(mynull)/2)])

