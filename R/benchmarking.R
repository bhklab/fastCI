# Generate random matrices, compute naive CI 

# Compare run times for everything:


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
  
  fastvals[i] <- fastCI(x,y, C=FALSE, CPP=FALSE)$cindex
  naivevals[i] <- naiveCI(x,y)$cindex
}
plot(naivevals, fastvals, pch=16)
print(paste("Number of conflicts in", length(fastvals), "iterations, CI with no ties:", sum(fastvals != naivevals)))


# Compare fastCI against naiveCI with ties, various methods of counting:
myiter <- 1000
naive00 <- numeric(myiter)
fast00 <- numeric(myiter)
for (i in seq(myiter)){
  if (i %% 100 == 0) {
    print(i)
  }
  mydim <- 5
  x <- sample(2*mydim, mydim, replace=TRUE)
  y <- sample(2*mydim, mydim, replace=TRUE)
  
  fast00[i] <- fastCI(x,y, C=FALSE, CPP=FALSE, discardObsTies=FALSE, discardPredTies=FALSE)$cindex
  naive00[i] <- naiveCI(x,y, count_xties = TRUE, count_yties = TRUE)$cindex
  
  if (fast00[i] != naive00[i]) {
    break
  }
}
plot(naive00, fast00, pch=16)
print(paste("Number of conflicts in", length(fast00), "iterations, CI with ties, discard=00:", sum(fast00 != naive00)))