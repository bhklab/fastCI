calculateDistributionEstimate(n, costs){
  
  if(length(costs) != choose(n,2)){
    print("Please define a cost for each pair")
  }
  costDist <- table(costs)
   
}