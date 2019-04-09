source("~/Code/fastCI/fastCI.R")
library(mCI)


nperm <- 100

ci.outx.false <- numeric(nperm)
ci.outx.true <- numeric(nperm)
ci.noised <- numeric(nperm)

for(i in 1:nperm){
  
  x <- sample.int(10, size=20, replace=TRUE)
  
  y <- x
  
  ## swap 5 random pairs
  
  for(j in 1:5){
    
    i1 <- sample(20,1)
    i2 <- sample(20,1)
    
    temp <- y[i1]
    y[i1] <- y[i2]
    y[i2] <- temp
      
  }
  
  ci.outx.false[i] <- paired.concordance.index(predictions = y, observations = x, outx = FALSE, delta.obs = 0, delta.pred = 0)[[1]]
  ci.outx.true[i] <- paired.concordance.index(predictions = y, observations = x, outx = TRUE, delta.obs = 0, delta.pred = 0)[[1]]
  
  dup.y <- duplicated(y)
  
  y[dup.y] <- y[dup.y] + rnorm(sum(dup.y),0,sqrt(.Machine$double.eps))
  
  
  dup.x <- duplicated(x)
  
  x[dup.x] <- x[dup.x] + rnorm(sum(dup.x),0,sqrt(.Machine$double.eps))
  
  ci.noised[i] <- paired.concordance.index(predictions = y, observations = x, outx = TRUE, delta.obs = 0, delta.pred = 0)[[1]]
  
}

plot(ci.outx.false-ci.noised)
plot(ci.outx.true-ci.noised)
cor(ci.outx.false, ci.noised)


# Testing effect on single vector

x <- sample.int(10, 20, replace=TRUE)

y <- x

## swap 5 random pairs

for(j in 1:5){
  
  i1 <- sample(20,1)
  i2 <- sample(20,1)
  
  temp <- y[i1]
  y[i1] <- y[i2]
  y[i2] <- temp
  
}


ci.outx.false <- numeric(nperm)
ci.outx.true <- numeric(nperm)
ci.noised <- numeric(nperm)

for(i in 1:100){
  
  x2 <- x
  y2 <- y
  
  ci.outx.false[i] <- paired.concordance.index(predictions = y2, observations = x2, outx = FALSE, delta.obs = 0, delta.pred = 0)[[1]]
  ci.outx.true[i] <- paired.concordance.index(predictions = y2, observations = x2, outx = TRUE, delta.obs = 0, delta.pred = 0)[[1]]
  
  dup.y2 <- duplicated(y2)
  
  y2[dup.y2] <- y2[dup.y2] + rnorm(sum(dup.y2),0,sqrt(.Machine$double.eps))
  
  
  dup.x2 <- duplicated(x2)
  
  x2[dup.x2] <- x2[dup.x2] + rnorm(sum(dup.x2),0,sqrt(.Machine$double.eps))
  
  ci.noised[i] <- paired.concordance.index(predictions = y2, observations = x2, outx = TRUE, delta.obs = 0, delta.pred = 0)[[1]]
  
}

hist(ci.noised, xlim = c(min(unlist(list(ci.noised, ci.outx.false, ci.outx.true))), max(unlist(list(ci.noised, ci.outx.false, ci.outx.true)))))
abline(v=ci.outx.false[1], col="blue", lwd=2)
abline(v=ci.outx.true[1], col="red", lwd=2)



