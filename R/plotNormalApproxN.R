plotNormalApproxN <-
function(N){
    
    plot.new()
    
    test <- nullCIDist(N,1,cumulative = 0, force_sym = 1)
    
    tt <- naive_approximate_normal_from_density(test)
    
    tt2 <- computeExpectedApproximation(N)
    
    plot(0:(length(test)-1), test, type='l', col='red')
    
    lines(0:(length(test)), dnorm(0:length(test), tt[1], tt[2]), col='blue')
    
    
    lines(0:(length(test)), dnorm(0:length(test), tt2[1], tt2[2]), col='green')
    
    
    legend("topleft", legend = c("True", "Fitted Approximation", "Expected Approximation"),fill = c("red", "blue", "green"))
}
