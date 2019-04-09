computeExpectedApproximation <-
function(N){
  
  mean <- choose(N, 2)/2 
  
  var <- sqrt((2*N^3 + 3*N^2 - 5*N)/72)
  return(c("mean" = mean, "var"=var))
}
