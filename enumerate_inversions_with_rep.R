### A script to investigate the effect of repeated elements on the distribution of CIs


library(gtools)
source("~/Code/fastCI/fastCI.R")

my.nums <- c(1,1,1,4,5, 5, 5, 8, 9, 10, 11, 12)

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
