

profvis::profvis({
  for(i in 1:5000) test <- fastCI(rnorm(100), rnorm(100))
  }, interval = 0.005)
