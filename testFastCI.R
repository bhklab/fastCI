source("~/Code/fastCI.R")
library(mCI)

resOuterList <- list()

for (j in 1:5){
  print(j)
  xList <- list(runif(100), runif(1000), runif(10000))#, runif(100000), runif(1000000))
  yList <- xList
  yList <- lapply(yList, function(x) return(rev(x + x*runif(length(x), -0.1, 0.1))))
  
  timeInnerListmCI <- numeric(5)
  
  for(i in 1:3) {
    timeInnerListmCI[i] <- system.time(paired.concordance.index(xList[[i]], yList[[i]], delta.pred = 0, delta.obs = 0))[[3]]
  }
  
  xList <- list(runif(100), runif(1000), runif(10000), runif(100000), runif(1000000))
  yList <- xList
  yList <- lapply(yList, function(x) return(rev(x + x*runif(length(x), -0.1, 0.1))))
  
  timeInnerListfastCI <- numeric(5)
  
  for(i in 1:5) {
    timeInnerListfastCI[i] <- system.time(fastCI(xList[[i]], yList[[i]]))[[3]]
  }
  resOuterList[[j]] <- list(mCI = timeInnerListmCI, fastCI = timeInnerListfastCI) 
} 

mCI.times <- lapply(resOuterList, function(x) return(data.frame("algo" = "mCI", "npoints" = c(100,1000,10000),"time" = x$mCI[1:3])))

mCI.times <- do.call(rbind, mCI.times)


fastCI.times <- lapply(resOuterList, function(x) return(data.frame("algo" = "fastCI", "npoints" = c(100,1000,10000, 100000, 1000000),"time" = x$fastCI[1:5])))

fastCI.times <- do.call(rbind, fastCI.times)


toPlot <- rbind(mCI.times, fastCI.times) 

require(ggplot2)

ggplot(toPlot, aes(npoints, time, colour=algo)) + geom_point() + scale_y_log10() + scale_x_log10() + geom_smooth(method="lm")


toPlot2 <- toPlot
toPlot2$npoints <- log10(toPlot2$npoints)
toPlot2$time <- log10(toPlot2$time)

summary(lm(time~npoints, subset(toPlot2, toPlot2$algo == "mCI")))

summary(lm(time~npoints, subset(toPlot2, toPlot2$algo == "fastCI")))



xList <- list(runif(1000), runif(1000), runif(1000), runif(1000), runif(1000))
yList <- xList
yList <- lapply(yList, function(x) return(rev(x + x*runif(length(x), -0.1, 0.1))))

for(i in 1:5) {
  fastCI.val <- fastCI(xList[[i]], yList[[i]])[[1]]
  mCI.val <- paired.concordance.index(xList[[i]], yList[[i]], delta.pred = 0, delta.obs = 0)[[1]]
  print(fastCI.val - mCI.val)
}


xList <- list(runif(1000), runif(1000), runif(1000), runif(1000), runif(1000))
yList <- list(runif(1000), runif(1000), runif(1000), runif(1000), runif(1000))
# yList <- lapply(yList, function(x) return(rev(x + x*runif(length(x), -0.1, 0.1))))

for(i in 1:5) {
  fastCI.val <- fastCI(xList[[i]], yList[[i]], outx=TRUE)[[1]]
  mCI.val <- paired.concordance.index(xList[[i]], yList[[i]], delta.pred = 0, delta.obs = 0, outx=TRUE)[[1]]
  print(fastCI.val - mCI.val)
}


for(i in 1:5) {
  fastCI.val <- fastCI(xList[[i]], yList[[i]])[3]
  mCI.val <- paired.concordance.index(xList[[i]], yList[[i]], delta.pred = 0, delta.obs = 0)[3]
  print(fastCI.val)
  print(mCI.val)
}

