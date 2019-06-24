
# require(matrixStats)
# dyn.load("~/Code/fastCI/fastCI.so")
# library(Rcpp)
# sourceCpp("~/Code/fastCI/fastCI.cpp")


merge_two_sides <- function(left, right, outx, delta.pred, delta.obs, logic.operator){
  # browser()
  # Unpacking the elements of the lists for convenience. 
  left_observations <- left[[1]]
  left_predictions <- left[[2]]
  left_discordant <- left[[3]]
  left_pairs <- left[[4]]
  
  right_observations <- right[[1]]
  right_predictions <- right[[2]]
  right_discordant <- right[[3]]
  right_pairs <- right[[4]]

  #RLR = Right List Removed
  RLR <- 0
  #LLL = Left List Left
  LLL <- length(left_observations)
  
  # Length Right
  LR <- length(right_observations)
  
  # Create output vectors of right length to iterate through
  out_observations <- numeric(LLL + LR)
  out_predictions <- numeric(LLL + LR)
  out_discordant <- numeric(LLL + LR)
  out_pairs <- numeric(LLL + LR)
  
  #Left Index; Right Index, index (of output vector)
  Li <- 1
  Ri <- 1
  i <- 1

  while(i <= length(out_observations)){
    
    if(LLL == 0){
      ## If left list is empty the only things we can do is fill in the
      ## output with right list elements.
      # RCI: fixed.  base condition, don't need to augment in any way
      out_observations[i] <- right_observations[Ri]
      out_predictions[i] <- right_predictions[Ri]
      out_discordant[i] <- right_discordant[Ri] + LLL #LLL = 0, but for consistency leaving here
      out_pairs[i] <- right_pairs[Ri]
      Ri <- Ri + 1
      i <- i + 1
      next
    }
    if(RLR == LR){
      ## If all elements from the right list have been removed, we fill in from left list
      #  IS: couldn't we just copy all remaining left list into out_*?  I guess it's not C-like
      #  RCI: fixed
      out_observations[i] <- left_observations[Li]
      out_predictions[i] <- left_predictions[Li]
      out_discordant[i] <- left_discordant[Li]
      out_pairs[i] <- left_pairs[Li]
      Li <- Li + 1
      i <- i + 1
      next
    }
    if(left_predictions[Li] == right_predictions[Ri] || (left_observations[Li] == right_observations[Ri] && outx)){
      # Is this still wrong?
      ## This loop removes elements from the left list while they remain tied with the leftmost element of the right list 
      # No RCI cost in the inner loop, but outside.  Anything with a cost of right_discordant + LLL must be augmented
      while(LLL && (left_observations[Li] == right_observations[Ri] || left_predictions[Li] == right_predictions[Ri])){
        out_observations[i] <- left_observations[Li]
        out_predictions[i] <- left_predictions[Li]
        out_discordant[i] <- left_discordant[Li]
        out_pairs[i] <- left_pairs[Li] - 1
        i <- i + 1
        LLL <- LLL - 1
        Li <- Li + 1
      }
      out_observations[i] <- right_observations[Ri]
      out_predictions[i] <- right_predictions[Ri]
      ###  Add RCI cost here
      out_discordant[i] <- right_discordant[Ri] + LLL 
      out_pairs[i] <- right_pairs[Ri] - 1
      RLR <- RLR + 1
      Ri <- Ri + 1
      i <- i + 1
    } else if(left_observations[Li] < right_observations[Ri]) {
      #  RCI: fixed
      out_observations[i] <- left_observations[Li]
      out_predictions[i] <- left_predictions[Li]
      out_discordant[i] <- left_discordant[Li]
      out_pairs[i] <- left_pairs[Li]
      LLL <- LLL - 1
      Li <- Li + 1
      i <- i + 1
    } else if(left_observations[Li] > right_observations[Ri]) {
      out_observations[i] <- right_observations[Ri]
      out_predictions[i] <- right_predictions[Ri]
      ###  Add RCI cost here
      out_discordant[i] <- right_discordant[Ri] + LLL
      out_pairs[i] <- right_pairs[Ri]
      RLR <- RLR + 1
      Ri <- Ri + 1
      i <- i + 1
    } else {
      # only case left is if the two values are equal and outx is false
      # stop("Not implemented correctly?")
      out_observations[i] <- left_observations[Li]
      out_predictions[i] <- left_predictions[Li]
      out_discordant[i] <- left_discordant[Li] + RLR + 0.5
      out_pairs[i] <- left_pairs[Li]
      i <- i + 1
      out_observations[i] <- right_observations[Ri]
      out_predictions[i] <- right_predictions[Ri]
      out_discordant[i] <- right_discordant[Ri] + LLL - 0.5
      out_pairs[i] <- right_pairs[Ri]
      LLL <- LLL - 1
      Li <- Li + 1
      RLR <- RLR + 1
      Ri <- Ri + 1
      i <- i + 1
    }
  }
  
  return(list(out_observations, out_predictions, out_discordant, out_pairs))
}


rci_binary_update <- function(element, numvec, discordant, threshold){
  # inputs: element, a number; element is less than the smallest (first) element of numvec
  #         numvec, a sorted vector of numbers
  #         discordant, a vector of integers of the same length as numvec
  
  L <- length(numvec)
  idx <- get_count(element+threshold, numvec, L)
  discordant[seq_len(L - idx) + idx] <- discordant[seq_len(L - idx) + idx] + 1
  list(discordant=discordant, cost=L-idx)
}

get_count <- function(value, numvec, veclength){
  ## Find the largest index i such that numvec[i] < value, given that numvec is sorted
  ## Binary search method is slow in R
  # mid <- ceil(veclength/2)
  # if (veclength == 0){
  #   return(0)
  # } else if (numvec[mid] >= value) {
  #   return( get_count(value, numvec[1:(mid-1)], mid-1) )
  # } else {
  #   return( mid + get_count(value, numvec[(mid+1):veclength], veclength-mid) )
  # }
  return(sum(numvec < value))
}

rmerge_sort <- function(input, outx, delta.pred, delta.obs, logic.operator){
  if(length(input[[1]]) == 1){
    return(input)
  } else {
    input_observations <- input[[1]]
    input_predictions <- input[[2]]
    input_discordant <- input[[3]]
    input_pairs <- input[[4]]
    split_idx <- floor(length(input_observations)/2)
    left <- list(input_observations[seq(1, split_idx)],
                 input_predictions[seq(1, split_idx)], 
                 input_discordant[seq(1, split_idx)],
                 input_pairs[seq(1, split_idx)])
    right <- list(input_observations[seq(split_idx+1, length(input_observations))],
                  input_predictions[seq(split_idx+1, length(input_predictions))], 
                  input_discordant[seq(split_idx+1, length(input_observations))],
                  input_pairs[seq(split_idx+1, length(input_observations))])
    left <- rmerge_sort(left, outx, delta.pred, delta.obs, logic.operator)
    right <- rmerge_sort(right, outx, delta.pred, delta.obs, logic.operator)
    output <- rmerge_two_sides(left, right, outx, delta.pred, delta.obs, logic.operator)
    return(output)
  }
}
## Currently, the following code gives prediction intervals for new CIs of the same sample size. 

fastRCI <- function(observations, predictions, outx = TRUE, alpha = 0.05, 
                    delta.pred = 0.2, delta.obs = 0.2, logic.operator = c("and", "or"), 
                    alternative = c("two.sided", "greater", "less"), 
                    interval = c("confidence", "prediction"), noise.ties = FALSE, 
                    noise.eps = sqrt(.Machine$double.eps), C = FALSE, CPP = FALSE){

  alternative = match.arg(alternative)
  interval = match.arg(interval)
  logic.operator <- match.arg(logic.operator)
  
  if(!length(observations) == length(predictions)){
    stop("Size of vectors must be the same")
  }
  
  myCompleteCases <- complete.cases(observations, predictions)
  observations <- observations[myCompleteCases]
  predictions <- predictions[myCompleteCases]
  
  myorder <- order(predictions, method = "radix")
  
  predictions <- predictions[myorder]
  observations <- observations[myorder]
  
  if(noise.ties){
    dup.pred <- duplicated(predictions)
    dup.obs <- duplicated(observations)
    
    ## Being extra-precautious about possible duplicates from rnorm. (VERY UNLIKELY)
    
    while(any(dup.obs) || any(dup.pred)){
      predictions[dup.pred] <- predictions[dup.pred] + rnorm(sum(dup.pred), 0, noise.eps)
      observations[dup.obs] <- observations[dup.obs] + rnorm(sum(dup.obs), 0, noise.eps)
      
      dup.pred <- duplicated(predictions)
      dup.obs <- duplicated(observations)
    }
  }

  if(C){
    discordant <- numeric(length(predictions))
    pairs <- rep(length(predictions)-1, length(predictions))
    output <- .Call("rmerge_sort_c", observations,
         predictions,
         discordant,
         pairs,
         length(observations), outx)
  } else {
      input <- list(observations, predictions, numeric(length(predictions)), rep(length(predictions)-1, length(predictions)))
      output <- rmerge_sort(input, outx, delta.pred, delta.obs, logic.operator)
    
  }

  output_discordant <- output[[3]]
  output_pairs <- output[[4]]
  comppairs=10
  # N <- length(predictions)
  # D <- exp(logSumExp(log(output_discordant)))
  # C <- exp(logSumExp(log((N-1)-output_discordant)))
  # CC <- exp(logSumExp(log(C) + log(C-1)))
  # DD <- exp(logSumExp(log(D) + log(D-1)))
  # CD <- exp(logSumExp(log(C) + log(D)))

  N <- length(predictions)
  D <- sum(output_discordant)
  Cvec <- output_pairs-output_discordant
  C <-  sum(Cvec)
  CC <- sum(Cvec*(Cvec-1))
  DD <- sum(output_discordant*(output_discordant-1))
  CD <- sum(Cvec*output_discordant)
  # browser()
  # if (N < 3 || (C == 0 && D == 0)) {
  #   return(list("cindex"=NA, "p.value"=NA, "sterr"=NA, "lower"=NA, "upper"=NA, "relevant.pairs.no"=0))
  # }
  # if(C==0 || D==0 || C * (C - 1)==0 || D * (D - 1)==0 || C * D==0 || (C + D) < comppairs){
  #   return(list("cindex"=NA, "p.value"=NA, "sterr"=NA, "lower"=NA, "upper"=NA, "relevant.pairs.no"=(C + D) / 2))
  # }
  # cindex <- exp(C) / exp(logSumExp(c(C, D)))
  cindex <- C/(C+D)
  varp <- 4 * ((D ^ 2 * CC - 2 * C * D * CD + C ^ 2 * DD) / (C + D) ^ 4) * N * (N - 1) / (N - 2) 
  
  # varp <- 4 * ((exp(logSumExp(c(2*D + CC, 2*C + DD))) - 2 *exp(C + D + CD)) / exp(logSumExp(c(C, D)))^4) * N * (N - 1) / (N - 2)
  
  if (varp >= 0) {
    sterr <- sqrt(varp / (N-1))
    if(interval == "confidence"){
      ci <- qnorm(p = alpha / 2, lower.tail = FALSE) * sterr
      p <- pnorm((cindex - 0.5) / sterr)
    } else{
      ci <- qnorm(p = alpha / 2, lower.tail = FALSE) * sterr
      p <- pnorm((cindex - 0.5) / sterr)
      ci <- qnorm(p = alpha/2, lower.tail = FALSE) * sterr * sqrt(2)
    }

  } else {
    return(list("cindex"=cindex,
                "p.value"=1,
                "sterr"=NA,
                "lower"=0,
                "upper"=0,
                "relevant.pairs.no"=(C + D) / 2))
  }
  return(list("cindex"=cindex,
              "p.value"=switch(alternative, less=p, greater=1 - p, two.sided=2 * min(p, 1 - p)),
              "sterr"=sterr,
              "lower"=max(cindex - ci, 0),
              "upper"=min(cindex + ci, 1),
              "relevant.pairs.no"=(C + D) / 2))
}


justFastRCI <- function(observations, predictions, outx = TRUE, 
                    delta.pred = 0.2, delta.obs = 0.2, logic.operator = c("and", "or"), 
                    noise.ties = FALSE, noise.eps = sqrt(.Machine$double.eps), C = TRUE, CPP = TRUE){
  # Incomplete

  logic.operator <- match.arg(logic.operator)

  if(!length(observations) == length(predictions)){
    stop("Size of vectors must be the same")
  }
  
  myCompleteCases <- complete.cases(observations, predictions)
  if(!sum(myCompleteCases)){
    return(c("CI" = NA, "N" =0))
  }
  observations <- observations[myCompleteCases]
  predictions <- predictions[myCompleteCases]
  
  
  myorder <- order(predictions, method = "radix")
  
  predictions <- predictions[myorder]
  observations <- observations[myorder]
  
  if(noise.ties){
    
    dup.pred <- duplicated(predictions)
    dup.obs <- duplicated(observations)
    
    ## Being extra-precautious about possible duplicates from rnorm. (VERY UNLIKELY)
    
    while(any(dup.obs) || any(dup.pred)){
      predictions[dup.pred] <- predictions[dup.pred] + rnorm(sum(dup.pred), 0, noise.eps)
      observations[dup.obs] <- observations[dup.obs] + rnorm(sum(dup.obs), 0, noise.eps)
      
      dup.pred <- duplicated(predictions)
      dup.obs <- duplicated(observations)
      
    }
  }
  if(C){
    discordant <- numeric(length(predictions))
    pairs <- rep(length(predictions)-1, length(predictions))
    
    # output_observations <- observations
    # output_predictions <- predictions
    # output_discordant <- discordant
    # output_pairs <- pairs
    
    # cres <- .C("merge_sort_c", as.double(observations),
    #                    as.double(predictions),
    #                    as.double(discordant),
    #                    as.double(pairs),
    #                    as.double(output_observations),
    #                    as.double(output_predictions),
    #                    as.double(output_discordant),
    #                    as.double(output_pairs), as.integer(length(observations)), as.integer(outx))
    # output <- cres[5:8]
    
    output <- .Call("merge_sort_c", observations,
                    predictions,
                    discordant,
                    pairs,
                    length(observations), outx)
    
  } else {
    if(CPP){
      output <- merge_sort_c(observations, predictions, numeric(length(predictions)), rep(length(predictions)-1, length(predictions)), outx)
    } else{
      input <- list(observations, predictions, numeric(length(predictions)), rep(length(predictions)-1, length(predictions)))
      
      output <- merge_sort(input, outx)
      
    }
    
  }
  
  output_discordant <- output[[3]]
  output_pairs <- output[[4]]

  
  N <- length(predictions)
  D <- sum(output_discordant)
  Cvec <- output_pairs-output_discordant
  C <-  sum(Cvec)


  cindex <- C/(C+D)
  
  return(c("CI" = cindex, "N" = N))
}




