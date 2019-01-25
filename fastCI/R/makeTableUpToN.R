makeTableUpToN <-
function(N){
  nullCIDistMem <- memoise(nullCIDistMem)
  res <- list()
  
  for(i in seq_len(N)){
    
    res[[i]] <- nullCIDistMem(i, multvect = 1, force_sym = TRUE, cumulative = FALSE)
    
  }
  return(res)
}
