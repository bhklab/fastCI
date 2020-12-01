# Compare the exact analytical null to a permutation null
# As CI is non-parametric, the distribution from which samples are drawn is irrelevant
# up to the tie structure.  

n <- 100
iter <- 1e6
myfunc <- rnorm
myargs <- list(mean=0, sd=1, n=n)
outdir <- "~/Work/bhk/analysis/mci_manuscript/nulldist/"
label <- "ci_rnorm"
dsname <- "Norm(0,1)"

ci_vals <- t(sapply(seq(iter), FUN=function(x) {
  x <- do.call(myfunc, myargs)
  y <- do.call(myfunc, myargs)
  paired.concordance.index.new(predictions = x, observations = y, delta.pred=0, delta.obs=0, outx=TRUE, logic.operator="and", 
                               CPP=TRUE, p_method="Asymptotic", conf_int_method="Asymptotic")[c("cindex", "p.value", "relevant.pairs.no")]
}))


ci_vals <- as.data.frame(matrix(sapply(ci_vals, unlist), ncol=3, dimnames=list(rownames(ci_vals), colnames(ci_vals))))


ci_null <- nullCIDist(n=100, cumulative=0)
ci_null_cum <- nullCIDist(n=100, cumulative=1)

# Multiplying the density plot by choose(n,2) is necessary because the density plot is discrete on 
# seq(0, choose(n,2)) concordant pairs.  The sum over ci_null is 1, but because the "bin" sizes are 
# 1/choose(n,2), it's necessary to rescale it to match the density plot of the permutation null
pdf(paste(outdir, sprintf("AnalyticalNullComparisontoPermutation_%s_n=%d_%dk.pdf", label, n, round(iter/1000)), sep="/"), width=10, height=8)
plot(seq(0,choose(n,2))/choose(n,2), ci_null*choose(n,2), type="l", lwd=2.5, col=rgb(1,0,0, alpha=0.6), 
     xlab="Concordance Index", ylab="Probability Density", main="Comparison of Analytical Null with Permutation Null", xlim=c(0.3, 0.7))
lines(density(ci_vals$cindex, bw=0.001), col=rgb(0,0.2,0.8, alpha=0.6), lwd=2.5)
legend(x="topright", legend=c("Analytical Null", "Permutation Null N=1e6"), col=c(rgb(1,0,0, alpha=0.6), rgb(0,0.2,0.8, alpha=0.6)), lwd=c(2.5,2.5))
dev.off()

ci_pvals <- sapply(ci_vals$cindex * choose(n,2), function(x) getCIPvals(ci_null_cum, x, alternative="two.sided"))
ci_rank <- rank(-abs(ci_vals$cindex - 0.5))/length(ci_vals$cindex)
maxy <- 1.05 * max(max(-log10(ci_vals$p.value)), -log10(ci_rank), -log10(ci_pvals))

pdf(paste(outdir, sprintf("QQAnalyticalNull_vs_Approximation_%s_n=%d_%dk.pdf", label, n, round(iter/1000)), sep="/"), width=10, height=8)
plot(-log10(ci_rank[ci_rank < 0.01]), -log10(ci_pvals[ci_rank < 0.01]), xlab="-log10 Empirical Quantile", ylab="-log10 Theoretical Quantile",
     main=sprintf("QQ plot of CI p-values of permutations of %s, n=%dk, length=%d", dsname, iter/1000, n), 
     col="red", pch=16, xlim=c(0, maxy), ylim=c(0, maxy))
points(-log10(ci_rank[ci_rank < 0.01]), -log10(ci_vals$p.value[ci_rank < 0.01]), pch=16, col="blue")
lines(c(0,10), c(0,10), lty=2, col="grey")
legend(x="bottomright", legend=c("Exact Analytical Null", "Asymptotic Null"), col=c("red", "blue"), pch=c(16,16))
grid()
dev.off()


     
     
     