
# require(matrixStats)
# dyn.load("~/Code/fastCI/fastCI.so")
# library(Rcpp)
# sourceCpp("~/Code/fastCI/fastCI.cpp")

# (deprecated) outx {boolean} set to TRUE to not count pairs of observations tied on either observations or predictions as a relevant pair. 
# discardObsTies, discardPredTies {boolean} set to true to not count pairs with tied elements in observations and predictions
# respectively in either pairs or concordant/discordant fields.  If false, a tied pair results in 1/2 added to concordant and
# discordant pairs. 

merge_two_sides_noties <- function(left, right){
  left_observations <- left[[1]]
  left_predictions <- left[[2]]
  left_discordant <- left[[3]]
  left_pairs <- left[[4]]
  
  right_observations <- right[[1]]
  right_predictions <- right[[2]]
  right_discordant <- right[[3]]
  right_pairs <- right[[4]]

  LL <- length(left_observations)
  LR <- length(right_observations)

  out_observations <- numeric(LL + LR)
  out_predictions <- numeric(LL + LR)
  out_discordant <- numeric(LL + LR)
  out_pairs <- numeric(LL + LR)
  
  #Left Index; Right Index, index (of output vector)
  Li <- 1
  Ri <- 1
  i <- 1
  
  while(i <= (LL + LR)) {
    if (Li > LL){
      # Left list empty - fill in from right list
      out_observations[i] <- right_observations[Ri]
      out_predictions[i] <- right_predictions[Ri]
      out_discordant[i] <- right_discordant[Ri] + LL - Li + 1  #LL-Li+1 = 0, but for consistency
      out_pairs[i] <- right_pairs[Ri]
      Ri <- Ri + 1
      i <- i + 1
      next
    } else if (Ri > LR) {
      ## If all elements from the right list have been removed, we fill in from left list
      out_observations[i] <- left_observations[Li]
      out_predictions[i] <- left_predictions[Li]
      out_discordant[i] <- left_discordant[Li] + Ri - 1        #Ri - 1 = LR, but for consistency
      out_pairs[i] <- left_pairs[Li]
      Li <- Li + 1
      i <- i + 1
      next
    } else if (left_observations[Li] < right_observations[Ri]) {
      out_observations[i] <- left_observations[Li]
      out_predictions[i] <- left_predictions[Li]
      out_discordant[i] <- left_discordant[Li] + Ri - 1
      out_pairs[i] <- left_pairs[Li]
      Li <- Li + 1
      i <- i + 1
    } else if (left_observations[Li] > right_observations[Ri]) {
      out_observations[i] <- right_observations[Ri]
      out_predictions[i] <- right_predictions[Ri]
      out_discordant[i] <- right_discordant[Ri] + LL - Li + 1
      out_pairs[i] <- right_pairs[Ri]
      Ri <- Ri + 1
      i <- i + 1
    } else {
      error("There is a tie when none should exist")
    } 
  }
  return(list(out_observations, out_predictions, out_discordant, out_pairs))
}
  


merge_two_sides <- function(left, right, discardTies){
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

  discardObsTies <- discardTies[1]
  discardPredTies <- discardTies[2]
  
  # Handle ties in the predictions
  LPredMax <- max(left_predictions)
  R_ix <- which(right_predictions == LPredMax)
  if (length(R_ix) > 0) {
    L_ix <- which(left_predictions == LPredMax)
    if (discardPredTies){
      right_pairs[R_ix] <- right_pairs[R_ix] - length(L_ix)
      left_pairs[L_ix] <- left_pairs[L_ix] - length(R_ix)
    } else {
      right_discordant[R_ix] <- right_discordant[R_ix] + 1/2 * length(L_ix)
      left_discordant[L_ix] <- left_discordant[L_ix] + 1/2 * length(R_ix)
    }
  }

  LpredTieCount <- 0
  RpredTieCount <- 0
  
  LL <- length(left_observations)
  LR <- length(right_observations)

  ## Create output vectors of right length to iterate through
  out_observations <- numeric(LL + LR)
  out_predictions <- numeric(LL + LR)
  out_discordant <- numeric(LL + LR)
  out_pairs <- numeric(LL + LR)
  
  #Left Index; Right Index, index (of output vector)
  Li <- 1
  Ri <- 1
  i <- 1
  
  while(i <= (LL + LR)) {     # length(out_observations)){
    if(Li > LL){
      ## If left list is empty the only things we can do is fill in the output with right list elements.
      out_observations[i] <- right_observations[Ri]
      out_predictions[i] <- right_predictions[Ri]
      out_discordant[i] <- right_discordant[Ri] + LL - Li + 1   #LLL #LLL = 0, but for consistency leaving here
      out_pairs[i] <- right_pairs[Ri]
      Ri <- Ri + 1
      i <- i + 1
      next
    }
    if(Ri > LR){
      ## If all elements from the right list have been removed, we fill in from left list
      out_observations[i] <- left_observations[Li]
      out_predictions[i] <- left_predictions[Li]
      out_discordant[i] <- left_discordant[Li] + Ri - 1
      out_pairs[i] <- left_pairs[Li]
      Li <- Li + 1
      i <- i + 1
      next
    } else if (left_observations[Li] < right_observations[Ri]) {
      out_observations[i] <- left_observations[Li]
      out_predictions[i] <- left_predictions[Li]
      # Need to account for ties in predictions from right list
      if (left_predictions[Li] == LPredMax) {
        # Number of R insertions minus ties in Pred.  Note: discPredTies is irrelevant
        # because the discardPredTies conditional above the while loop ALREADY handled ALL pred ties
        # correctly.  We only need to count the number of discordant elements that AREN'T tied in pred. 
        out_discordant[i] <- left_discordant[Li] + (Ri - 1) - RpredTieCount
        LpredTieCount <- LpredTieCount + 1
      } else {
        out_discordant[i] <- left_discordant[Li] + (Ri - 1)
      }
      out_pairs[i] <- left_pairs[Li]
      Li <- Li + 1
      i <- i + 1
    } else if (left_observations[Li] > right_observations[Ri]) {
      out_observations[i] <- right_observations[Ri]
      out_predictions[i] <- right_predictions[Ri]
      if (right_predictions[Ri] == LPredMax) {
        # Account for ties in predictions from left list
        # Discordant elements are those remaining in left list minus any remaining left elements with predTies
        out_discordant[i] <- right_discordant[Ri] + (LL - Li + 1) - (length(L_ix) - LpredTieCount)
        RpredTieCount <- RpredTieCount + 1
      } else {
        out_discordant[i] <- right_discordant[Ri] + (LL - Li + 1)
      } 
      out_pairs[i] <- right_pairs[Ri]
      Ri <- Ri + 1
      i <- i + 1
    } elseif (left_observations[Li] == right_observations[Ri]) {
      # Tie in observations; need to count how elements are tied in obs
      L_start <- Li
      R_start <- Ri
      L_end <- get_elt_count(left_observations, Li, left_observations[Li])
      R_end <- get_elt_count(right_observations, Ri, right_observations[Ri])
      L_RangePredMax <- sum(left_predictions[L_start:L_end] == LPredMax)
      R_RangePredMax <- sum(right_predictions[R_start:R_end] == RPredMax)

      # There is a minor bug here: elements tied in BOTH pred and obs defer to how they are handled for pred
      # This probably shouldn't happen, but there aren't any guarantees at present
      if (discardObsTies) {
        while (Li <= L_end){
          out_observations[i] <- left_observations[Li]
          out_predictions[i] <- left_predictions[Li]
          if (left_predictions[Li] == LPredMax) {
            out_discordant[i] <- left_discordant[Li] + (Ri-1) - RpredTieCount
            LpredTieCount <- LpredTieCount + 1
            out_pairs[i] <- left_pairs[Li] - (R_end - R_start + 1 - R_RangePredMax)
          } else {
            out_discordant[i] <- left_discordant[Li] + (Ri - 1)
            out_pairs[i] <- left_pairs[Li] - (R_end - R_start + 1)
          } 
          Li <- Li + 1
          i <- i + 1
        } 
        while (Ri <= R_end) {
          out_observations[i] <- right_observations[Ri]
          out_predictions[i] <- right_predictions[Ri]
          if (right_predictions[Ri] == LPredMax) {
            out_discordant[i] <- right_discordant[Ri] + (LL - Li + 1) - (length(L_ix) - LpredTieCount)
            RpredTieCount <- RpredTieCount + 1
            out_pairs[i] <- right_pairs[Ri] - (L_end - L_start + 1 - L_RangePredMax)
          } else {
            out_discordant[i] <- right_discordant[Ri] + (LL - Li + 1) 
            out_pairs[i] <- right_pairs[Ri] - (L_end - L_start + 1)
          }
          Ri <- Ri + 1
          i <- i + 1
        }
      else { #discardObsTies = F
        while (Li <= L_end) {
          out_observations[i] <- left_observations[Li]
          out_predictions[i] <- left_predictions[Li]
          if (left_predictions[Li] == LPredMax) {
            out_discordant[i] <- left_discordant[Li] + (Ri - 1) - RpredTieCount + (1/2)*(R_end - R_start + 1) - R_RangePredMax
            LpredTieCount <- LpredTieCount + 1
            out_pairs[i] <- left_pairs[Li]
          } else {
            out_discordant[i] <- left_discordant[Li] + (Ri - 1) + (1/2) *(R_end - R_start + 1)
            out_pairs[i] <- left_pairs[Ri]
          }
          Li <- Li + 1
          i <- i + 1
        }
        while (Ri  <= R_end) {
          out_observations[i] <- right_observations[Ri]
          out_predictions[i] <- right_predictions[Ri]
          if (right_predictions[Ri] == LPredMax) {
            out_discordant[i] <- right_discordant[Ri] + (LL - Li + 1) - (length(L_ix) - LpredTieCount) + (1/2) * (L_end - L_start + 1) - L_RangePredMax
            RpredTieCount <- RpredTieCount + 1
            out_pairs[i] <- right_pairs[Ri]
          } else {
            out_discordant[i] <- right_discordant[Ri] + (LL - Li + 1) + (1/2) * (L_end - L_start + 1)
            out_pairs[i] <- right_pairs[Ri]
          }

          
      }
    }
  elseif ((left_observations[Li] == right_observations[Ri]) || right_predictions[Ri] <= Lmax_pred) {  
  #     # There is a tie somewhere; deal with it appropriately
  #     # Ri is smallest obs but has a tie in obs 
  #     if((left_predictions[Li] == right_predictions[Ri] && discardPredTies) || (left_observations[Li] == right_observations[Ri] && discardObsTies)){
  #     ## This loop removes elements from the left list while they remain tied with the leftmost element of the right list 
  #     while(LLL && (left_observations[Li] == right_observations[Ri] && discardPredTies) || (left_predictions[Li] == right_predictions[Ri] && discardObsTies))){
  #       out_observations[i] <- left_observations[Li]
  #       out_predictions[i] <- left_predictions[Li]
  #       out_discordant[i] <- left_discordant[Li] + RLR 
  #       out_pairs[i] <- left_pairs[Li] - 1
  #       i <- i + 1
  #       LLL <- LLL - 1
  #       Li <- Li + 1
  #     }
  #     out_observations[i] <- right_observations[Ri]
  #     out_predictions[i] <- right_predictions[Ri]
  #     out_discordant[i] <- right_discordant[Ri] + LLL 
  #     out_pairs[i] <- right_pairs[Ri] - 1
  #     RLR <- RLR + 1
  #     Ri <- Ri + 1
  #     i <- i + 1
  #     # the two values are equal and ties are not to be discarded
  #     
  #     # Count how many ties there are on the left side.  
  #     # This creates a problem because we don't know which of the observations are tied; all right_obs >= left_obs, but 
  #     # the lists themselves are unsorted.  One possible solution is to first count ties in the observed lists, but 
  #     # we don't know which elements are tied in BOTH observed and predictions.  
  #     
  #     if (right_observations[Ri] > max(left_observations[Li:Li+LLL]) {
  #        # Only ties in pred, which is a sorted list on each side.
  #        # Count how many elements in left_pred tie with right_pred[Ri]
  #        tiecount <- 0
  #        while (left_predictions[Li + tiecount] == right_predictions[Ri]){
  #          tiecount <- tiecount + 1
  #        } 
  #        
  #     } else {
  #        # Ties in either
  #     }
  #     ln <- count_front_duplicates(left_observations, left_predictions, Li, Ri, )
  #     print("In last case of loop - outx = FALSE")
  #     out_observations[i] <- left_observations[Li]
  #     out_predictions[i] <- left_predictions[Li]
  #     out_discordant[i] <- left_discordant[Li] + RLR + 0.5
  #     out_pairs[i] <- left_pairs[Li]
  #     i <- i + 1
  #     out_observations[i] <- right_observations[Ri]
  #     out_predictions[i] <- right_predictions[Ri]
  #     out_discordant[i] <- right_discordant[Ri] + LLL - 0.5
  #     out_pairs[i] <- right_pairs[Ri]
  #     LLL <- LLL - 1
  #     Li <- Li + 1
  #     RLR <- RLR + 1
  #     Ri <- Ri + 1
  #     i <- i + 1
  #   } else if(left_observations[Li] > right_observations[Ri]) {
  #     out_observations[i] <- right_observations[Ri]
  #     out_predictions[i] <- right_predictions[Ri]
  #     out_discordant[i] <- right_discordant[Ri] + LLL
  #     out_pairs[i] <- right_pairs[Ri]
  #     RLR <- RLR + 1
  #     Ri <- Ri + 1
  #     i <- i + 1
  #   }
  # }
  # 
  # return(list(out_observations, out_predictions, out_discordant, out_pairs))
}


get_elt_count <- function(myvector, init){
  final <- init
  while (myvector[final] == myvector[init] && final <= length(myvector)){
    final <- final + 1
  }
  final
}


merge_sort_noties <- function(input){
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
    left <- merge_sort_noties(left)
    right <- merge_sort_noties(right)
    output <- merge_two_sides_noties(left, right)
    return(output)
  }
}

merge_sort <- function(input, discardTies){
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
    left <- merge_sort(left, discardTies)
    right <- merge_sort(right, discardTies)
    output <- merge_two_sides(left, right, discardTies)
    return(output)
  }
}

## Currently, the following code gives prediction intervals for new CIs of the same sample size. 

fastCI <- function(observations, predictions, discardObsTies = TRUE, discardPredTies = TRUE, alpha = 0.05, 
                   alternative = c("two.sided", "greater", "less"), interval = c("confidence", "prediction"), 
                   noise.ties = FALSE, noise.eps = sqrt(.Machine$double.eps), C = TRUE, CPP = TRUE){

  discardTies = c(discardObsTies, discardPredTies)
  alternative = match.arg(alternative)
  interval = match.arg(interval)
  
  if(!length(observations) == length(predictions)){
    stop("Size of vectors must be the same")
  }
  
  myCompleteCases <- complete.cases(observations, predictions)
  observations <- observations[myCompleteCases]
  predictions <- predictions[myCompleteCases]
  if(!sum(myCompleteCases)){
    return(c("CI" = NA, "N" =0))
  }
  
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
    output <- .Call("merge_sort_c", observations,
          predictions,
          discordant,
          pairs,
          length(observations), discardTies)
  } else {
      input <- list(observations, predictions, numeric(length(predictions)), rep(length(predictions)-1, length(predictions)))
      if (length(unique(observations)) == length(unique(predictions)) && length(unique(observations)) == length(observations)){
        output <- merge_sort_noties(input)
      }
      else {
        output <- merge_sort(input, discardTies)
      }
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


justFastCI <- function(observations, predictions, discardObsTies = TRUE, discardPredTies=TRUE,
                   noise.ties = FALSE, noise.eps = sqrt(.Machine$double.eps), C = TRUE, CPP = TRUE){
  
  discardTies = c(discardObsTies, discardPredTies)
  if(!length(observations) == length(predictions)){
    stop("Size of vectors must be the same")
  }
  
  myCompleteCases <- complete.cases(observations, predictions)
  if(!sum(myCompleteCases)){
  if(!sum(myCompleteCases)){
    return(c("CI" = NA, "N" =0))
  }
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
      
      if (length(unique(observations)) == length(unique(predictions)) && length(unique(observations)) == length(observations)){
        output <- merge_sort_noties(input, discardTies)
      }
      else {
        output <- merge_sort(input, discardTies)
      }
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




