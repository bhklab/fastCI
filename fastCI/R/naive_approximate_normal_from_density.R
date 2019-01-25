naive_approximate_normal_from_density <-
function(x){
  
  mean <- sum((seq_along(x) - 1)*x)
  
  sd = sqrt(sum((seq_along(x) - 1)^2*x) - sum((seq_along(x) - 1)*x)^2)
  
  return(c("mean" = mean, "sd"=sd))
  
}
