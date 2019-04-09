
require(matrixStats)

source("~/Code/fastCI/computeNullCI.R")
nullCIDist <- memoise::memoise(nullCIDist)

merge_two_sides <- function(left, right, outx){
  
  left_vals <- left[[1]]
  left_discordant <- left[[2]]
  left_pairs <- left[[3]]
  
  right_vals <- right[[1]]
  right_discordant <- right[[2]]
  right_pairs <- right[[3]]
  
  RLR <- 0
  LLL <- length(left_vals)
  
  LR <- length(right_vals)
  
  out_vals <- numeric(LLL + LR)
  out_discordant <- numeric(length(out_vals))
  out_pairs <- numeric(length(out_vals))
  
  Li <- 1
  Ri <- 1
  i <- 1
  while(i <= length(out_vals)){
    
    if(LLL == 0){
      #Break out of loop if left list is empty
      out_vals[i] <- right_vals[Ri]
      out_discordant[i] <- right_discordant[Ri] + LLL
      out_pairs[i] <- right_pairs[Ri]
      Ri <- Ri + 1
      i <- i + 1
      next
    }
    if(RLR == LR){
      #Break out of loop if right list is empty
      out_vals[i] <- left_vals[Li]
      out_discordant[i] <- left_discordant[Li] + RLR
      out_pairs[i] <- left_pairs[Li]
      Li <- Li + 1
      i <- i + 1
      next
    }
    
    if(left_vals[Li] < right_vals[Ri]) {
      out_vals[i] <- left_vals[Li]
      out_discordant[i] <- left_discordant[Li] + RLR
      out_pairs[i] <- left_pairs[Li]
      LLL <- LLL - 1
      Li <- Li + 1
      i <- i + 1
    } else if(left_vals[Li] > right_vals[Ri]) {
      out_vals[i] <- right_vals[Ri]
      out_discordant[i] <- right_discordant[Ri] + LLL
      out_pairs[i] <- right_pairs[Ri]
      RLR <- RLR + 1
      Ri <- Ri + 1
      i <- i + 1
    } else {
      # only case left is if the two values are equal.
      if(outx){
        out_vals[i] <- left_vals[Li]
        out_discordant[i] <- left_discordant[Li] + RLR 
        out_pairs[i] <- left_pairs[Li] - 1
        i <- i + 1
        out_vals[i] <- right_vals[Ri]
        out_discordant[i] <- right_discordant[Ri] + LLL - 1
        out_pairs[i] <- right_pairs[Ri] - 1
        LLL <- LLL - 1
        Li <- Li + 1
        RLR <- RLR + 1
        Ri <- Ri + 1
        i <- i + 1
      } else {
        # stop("Not implemented correctly?")
        out_vals[i] <- left_vals[Li]
        out_discordant[i] <- left_discordant[Li] + RLR + 0.5
        out_pairs[i] <- left_pairs[Li]
        i <- i + 1
        out_vals[i] <- right_vals[Ri]
        out_discordant[i] <- right_discordant[Ri] + LLL - 0.5
        out_pairs[i] <- right_pairs[Ri]
        LLL <- LLL - 1
        Li <- Li + 1
        RLR <- RLR + 1
        Ri <- Ri + 1
        i <- i + 1
      }
    }
  }
  
  return(list(out_vals, out_discordant, out_pairs))
  
}


merge_sort <- function(input, outx){
  if(length(input[[1]]) == 1){
    return(input)
  } else {
    input_vals <- input[[1]]
    input_discordant <- input[[2]]
    input_pairs <- input[[3]]
    split_idx <- floor(length(input_vals)/2)
    left <- list(input_vals[seq(1, split_idx)], 
                 input_discordant[seq(1, split_idx)],
                 input_pairs[seq(1, split_idx)])
    right <- list(input_vals[seq(split_idx+1, length(input_vals))], 
                  input_discordant[seq(split_idx+1, length(input_vals))],
                  input_pairs[seq(split_idx+1, length(input_vals))])
    left <- merge_sort(left, outx)
    right <- merge_sort(right, outx)
    output <- merge_two_sides(left, right, outx)
    return(output)
  }
}
## TODO: this code does not handle missing values well. 

fastCI <- function(observations, predictions, outx = TRUE, alpha = 0.05, alternative = c("two.sided", "greater", "less")){
  
  alternative = match.arg(alternative)
  if(!length(observations) == length(predictions)){
    stop("Size of vectors must be the same")
  }
  
  
  myCompleteCases <- complete.cases(observations, predictions)
  observations <- observations[myCompleteCases]
  predictions <- predictions[myCompleteCases]
  
  
  myorder <- order(observations)
  
  predictions <- predictions[myorder]
  
  input <- list(predictions, numeric(length(predictions)), rep(length(predictions)-1, length(predictions)))
  output <- merge_sort(input, outx)
  output_discordant <- output[[2]]
  output_pairs <- output[[3]]
  comppairs=10
  
  N <- length(predictions)
  D <- sum(output_discordant)
  Cvec <- output_pairs-output_discordant
  C <-  sum(Cvec)
  # browser()
  if (N < 3 || (C == 0 && D == 0)) {
    return(list("cindex"=NA, "p.value"=NA, "sterr"=NA, "lower"=NA, "upper"=NA, "relevant.pairs.no"=0))
  }
  cindex <- C/(C+D)
  
  nullDist <- nullCIDist(N)
  
 # browser() 
  p <- getCIPvals(nullDist, N, D/2, alternative = alternative)

  return(list("cindex"=cindex,
              "p.value"=p,
              "relevant.pairs.no"=(C + D) / 2))
}
