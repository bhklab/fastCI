### A script to investigate the effect of repeated elements on the distribution of CIs


library(gtools)
source("~/Code/fastCI/fastCI.R")

my.nums <- c(1,1,1,1,5, 5, 5, 8, 9, 10, 11, 12)

for (i in 7:7){
  myx <- permutations(i, i)
  # print(myx)
  res <- apply(myx, 1, function(myx){
    input <- list(my.nums[myx], numeric(length(my.nums[myx])), rep(length(my.nums[myx])-1, length(my.nums[myx])))
    output <- merge_sort(input, TRUE)
    # browser()
    output_discordant <- output[[2]]
    return(sum(output_discordant)/2)
  })
  print(table(res))
  hist(res) 
  # print(res)
}


source("~/Code/fastCI/computeNullCI.R")


makeNullCIDist(7, c(1,3,3),cumulative = 0)*factorial(7)


## Need to figure out what the hell is the distibution for this: x = [1 2 2 3 3 4] and y = [5 5 6 6 7 8]

x <- c(1,2,2,3,3,4)
y <- c(5,5,6,6,7,8)
computeInversionsNumber <- function(x,y){
  # order.y <- order(y)
  # y <- y[order.y]
  # x <- x[order.y]
  numInv <- 0
  for(i in 1:length(x)){
    for(j in i:length(x)){
      if(x[i] > x[j] & y[i] != y[j]){
        numInv = numInv + 1
      }
    }
  }
  return(numInv)
}


for (i in 6:6){
  myx <- permutations(i, i)
  # print(myx)
  res <- apply(myx, 1, function(myx){
    cur.x <- x[myx]
    cur.y <- y[myx]
    return(computeInversionsNumber(cur.x, cur.y))
  })
  print(table(res))
  hist(res) 
  # print(res)
}




