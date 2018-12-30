
require(matrixStats)


merge_two_sides <- function(left, right, outx){
  # browser()
  left_observations <- left[[1]]
  left_predictions <- left[[2]]
  left_discordant <- left[[3]]
  left_pairs <- left[[4]]
  
  right_observations <- right[[1]]
  right_predictions <- right[[2]]
  right_discordant <- right[[3]]
  right_pairs <- right[[4]]
  
  RLR <- 0
  LLL <- length(left_observations)
  
  LR <- length(right_observations)
  
  out_observations <- numeric(LLL + LR)
  out_predictions <- numeric(LLL + LR)
  out_discordant <- numeric(length(out_observations))
  out_pairs <- numeric(length(out_observations))
  
  Li <- 1
  Ri <- 1
  i <- 1
  while(i <= length(out_observations)){
    
    if(LLL == 0){
      #Break out of loop if left list is empty
      out_observations[i] <- right_observations[Ri]
      out_predictions[i] <- right_predictions[Ri]
      out_discordant[i] <- right_discordant[Ri] + LLL
      out_pairs[i] <- right_pairs[Ri]
      Ri <- Ri + 1
      i <- i + 1
      next
    }
    if(RLR == LR){
      #Break out of loop if right list is empty
      out_observations[i] <- left_observations[Li]
      out_predictions[i] <- left_predictions[Li]
      out_discordant[i] <- left_discordant[Li] + RLR
      out_pairs[i] <- left_pairs[Li]
      Li <- Li + 1
      i <- i + 1
      next
    }
    if(left_predictions[Li] == right_predictions[Ri] || (left_observations[Li] == right_observations[Ri] && outx)){
      ## This should be split into two cases, one that counts tied predictions, and one that counts tied observations. 
      current_observation <- left_observations[Li]
      current_prediction <- left_predictions[Li]
      ## TODO: This out_pair counting below is incorrect. Maybe we should just count up for each valid comparison made?
      while(LLL && (left_observations[Li] == current_observation || left_predictions[Li] == current_prediction)){
        out_observations[i] <- left_observations[Li]
        out_predictions[i] <- left_predictions[Li]
        out_discordant[i] <- left_discordant[Li] + RLR 
        out_pairs[i] <- left_pairs[Li] - 1
        i <- i + 1
        LLL <- LLL - 1
        Li <- Li + 1
      }
      out_observations[i] <- right_observations[Ri]
      out_predictions[i] <- right_predictions[Ri]
      out_discordant[i] <- right_discordant[Ri] + LLL 
      out_pairs[i] <- right_pairs[Ri] - 1
      RLR <- RLR + 1
      Ri <- Ri + 1
      i <- i + 1
    } else if(left_observations[Li] < right_observations[Ri]) {
      out_observations[i] <- left_observations[Li]
      out_predictions[i] <- left_predictions[Li]
      out_discordant[i] <- left_discordant[Li] + RLR
      out_pairs[i] <- left_pairs[Li]
      LLL <- LLL - 1
      Li <- Li + 1
      i <- i + 1
    } else if(left_observations[Li] > right_observations[Ri]) {
      out_observations[i] <- right_observations[Ri]
      out_predictions[i] <- right_predictions[Ri]
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


merge_sort <- function(input, outx){
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
    left <- merge_sort(left, outx)
    right <- merge_sort(right, outx)
    output <- merge_two_sides(left, right, outx)
    return(output)
  }
}
## Currently, the following code gives prediction intervals for new CIs of the same sample size. 

fastCI <- function(observations, predictions, outx = TRUE, alpha = 0.05, alternative = c("two.sided", "greater", "less"), interval = c("confidence", "prediction")){

  alternative = match.arg(alternative)
  interval = match.arg(interval)
  
  if(!length(observations) == length(predictions)){
    stop("Size of vectors must be the same")
  }
  
  
  myCompleteCases <- complete.cases(observations, predictions)
  observations <- observations[myCompleteCases]
  predictions <- predictions[myCompleteCases]

  
  myorder <- order(predictions, method = "radix")
  
  predictions <- predictions[myorder]
  observations <- observations[myorder]

  input <- list(observations, predictions, numeric(length(predictions)), rep(length(predictions)-1, length(predictions)))
  output <- merge_sort(input, outx)
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
