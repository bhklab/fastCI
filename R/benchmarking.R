# Generate random matrices, compute naive CI 

# Compare fastCI against naiveCI for no ties, R code:
myiter <- 1000
naivevals <- numeric(myiter)
fastvals <- numeric(myiter)
for (i in seq(myiter)){
  if (i %% 100 == 0) {
    print(i)
  }
  mydim <- 100
  z <- rnorm(mydim)
  x <- sample(z, mydim/2)
  y <- setdiff(z, x)
  
  fastvals[i] <- fastCI(x,y)$cindex
  naivevals[i] <- naiveCI(x,y)$cindex
}


