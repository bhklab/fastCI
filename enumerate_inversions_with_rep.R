### A script to investigate the effect of repeated elements on the distribution of CIs


library(gtools)
source("~/Code/fastCI/fastCI.R")

my.nums <- c(1,1,1,4,5, 6, 7, 8, 9, 10, 11, 12)

for (i in 1:8){
  myx <- permutations(i, i)
  # print(myx)
  res <- apply(myx, 1, function(myx){
    input <- list(my.nums[myx], numeric(length(my.nums[myx])), rep(length(my.nums[myx])-1, length(my.nums[myx])))
    output <- merge_sort(input, TRUE)
    output_discordant <- output[[2]]
    return(sum(output_discordant)/2)
  })
  print(table(res))
  hist(res) 
}




