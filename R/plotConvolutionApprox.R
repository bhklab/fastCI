plotConvolutionApprox <-
function(n1, n2){
  
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
