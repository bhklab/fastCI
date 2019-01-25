source("~/Code/fastCI/computeNullCI.R")
library(memoise)
nullCIDist <- memoise(nullCIDist)

makeTableUpToN <- function(N){
  
  res <- list()
  
  for(i in seq_len(N)){
    
    res[[i]] <- nullCIDist(i, multvect = 1, force_sym = TRUE, cumulative = FALSE)
    
  }
  return(res)
}
