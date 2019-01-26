#include <Rcpp.h>
using namespace Rcpp;

#include <algorithm>
#include <cmath>
#include <random>
// #include <Rcpp.h>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <vector>
#include <string>
#include <list>
#include <iostream>


/* function returns sign of argument */
int sgn(double x) {
  if (x > 0) {
    return(1);
  } else if (x == 0) {
    return(0);
  } else if (x < 0) {
    return(-1);
  } else {
    throw std::invalid_argument ("erfinv received an argument that was not a real number.");
  }
}

/* function computes inverse of error function for normal quantile calculation */
double erfinv(double x){
  int sgnx = sgn(x);
  x = 1 - x * x;
  double lnx = log(x);
  double a = 4.3307467508 + lnx / 2; // magic number for Winitzki approximation
  double b = 6.80272108844 * lnx; // other magic number for approximation
  return(sgnx * sqrt(sqrt(a * a - b) - a));
}


//function to do and, or
double logicOpF(bool x, bool y, std::string logicOp){
  if (logicOp.compare("and") == 0){
    return(x && y);
  }else if(logicOp.compare("or") == 0){
    return(x || y);
  }

  throw std::invalid_argument ("CI calculation failed.");

}




std::list<std::vector<double> > merge_two_sides(std::list<std::vector<double> > left, std::list<std::vector<double> > right, bool outx){


	std::vector<double> left_observations = left.front(); left.pop_front();
	std::vector<double> left_predictions = left.front(); left.pop_front();
	std::vector<double> left_discordant = left.front(); left.pop_front();
	std::vector<double> left_pairs = left.front(); left.pop_front();

	std::vector<double> right_observations = right.front(); right.pop_front();
	std::vector<double> right_predictions = right.front(); right.pop_front();
	std::vector<double> right_discordant = right.front(); right.pop_front();
	std::vector<double> right_pairs = right.front(); right.pop_front();


  //RLR = Right List Removed  
  int RLR = 0;
  //LLL = Left List Left
  int LLL = left_observations.size();
  
  // Length Right
  int LR = right_observations.size();
  
  std::list<std::vector<double> > out;

  // Create output vectors of right length to iterate through
  std::vector<double> out_observations(LLL + LR);
  std::vector<double> out_predictions(LLL + LR);
  std::vector<double> out_discordant(LLL + LR);
  std::vector<double> out_pairs(LLL + LR);
  

  //Left Index; Right Index, index (of output vector)
  int Li = 0;
  int Ri = 0;
  int i = 0;


  while(i < out_observations.size()){
    
    if(LLL == 0){
      //// If left list is empty the only things we can do is fill in the
      //// output with right list elements.
      out_observations[i] = right_observations[Ri];
      out_predictions[i] = right_predictions[Ri];
      out_discordant[i] = right_discordant[Ri] + LLL; //LLL = 0, but for consistency leaving here;
      out_pairs[i] = right_pairs[Ri];
      Ri = Ri + 1;
      i = i + 1;
      continue;
    }
    if(RLR == LR){
      //// If all elements from the right list have been removed, we fill in from left list
      out_observations[i] = left_observations[Li];
      out_predictions[i] = left_predictions[Li];
      out_discordant[i] = left_discordant[Li] + RLR;
      out_pairs[i] = left_pairs[Li];
      Li = Li + 1;
      i = i + 1;
      continue;
    }
    if(left_predictions[Li] == right_predictions[Ri] || (left_observations[Li] == right_observations[Ri] && outx)){
      // Is this still wrong?
      //// This loop removes elements from the left list while they remain tied with the leftmost element of the right list 
      while(LLL && (left_observations[Li] == right_observations[Ri] || left_predictions[Li] == right_predictions[Ri])){
        out_observations[i] = left_observations[Li];
        out_predictions[i] = left_predictions[Li];
        out_discordant[i] = left_discordant[Li] + RLR;
        out_pairs[i] = left_pairs[Li] - 1;
        i = i + 1;
        LLL = LLL - 1;
        Li = Li + 1;
      }
      out_observations[i] = right_observations[Ri];
      out_predictions[i] = right_predictions[Ri];
      out_discordant[i] = right_discordant[Ri] + LLL ;
      out_pairs[i] = right_pairs[Ri] - 1;
      RLR = RLR + 1;
      Ri = Ri + 1;
      i = i + 1;
    } else if(left_observations[Li] < right_observations[Ri]) {
      out_observations[i] = left_observations[Li];
      out_predictions[i] = left_predictions[Li];
      out_discordant[i] = left_discordant[Li] + RLR;
      out_pairs[i] = left_pairs[Li];
      LLL = LLL - 1;
      Li = Li + 1;
      i = i + 1;
    } else if(left_observations[Li] > right_observations[Ri]) {
      out_observations[i] = right_observations[Ri];
      out_predictions[i] = right_predictions[Ri];
      out_discordant[i] = right_discordant[Ri] + LLL;
      out_pairs[i] = right_pairs[Ri];
      RLR = RLR + 1;
      Ri = Ri + 1;
      i = i + 1;
    } else {
      // only case left is if the two values are equal and outx is false
      // stop("Not implemented correctly?")
      out_observations[i] = left_observations[Li];
      out_predictions[i] = left_predictions[Li];
      out_discordant[i] = left_discordant[Li] + RLR + 0.5;
      out_pairs[i] = left_pairs[Li];
      i = i + 1;
      out_observations[i] = right_observations[Ri];
      out_predictions[i] = right_predictions[Ri];
      out_discordant[i] = right_discordant[Ri] + LLL - 0.5;
      out_pairs[i] = right_pairs[Ri];
      LLL = LLL - 1;
      Li = Li + 1;
      RLR = RLR + 1;
      Ri = Ri + 1;
      i = i + 1;
    }
  }

  out.push_back(out_observations);
  out.push_back(out_predictions);
  out.push_back(out_discordant);
  out.push_back(out_pairs);

  return out;

}


std::list<std::vector<double> > merge_sort(std::list<std::vector<double> > input, bool outx){
  if(input.front().size() == 1){
    return(input);
  } else {


    std::vector<double> input_observations = input.front(); input.pop_front();
    std::vector<double> input_predictions = input.front(); input.pop_front();
    std::vector<double> input_discordant = input.front(); input.pop_front();
    std::vector<double> input_pairs = input.front(); input.pop_front();
    
    std::list<std::vector<double> > left;
    std::list<std::vector<double> > right;

    std::list<std::vector<double> > output;


    int split_idx = floor(input_observations.size()/2);
    
    int input_size = input_observations.size();

    std::vector<double> left_observations(&input_observations[0], &input_observations[split_idx]);
    std::vector<double> left_predictions(&input_predictions[0], &input_predictions[split_idx]);
    std::vector<double> left_discordant(&input_discordant[0], &input_discordant[split_idx]);
    std::vector<double> left_pairs(&input_pairs[0], &input_pairs[split_idx]);

    std::vector<double> right_observations(&input_observations[split_idx], &input_observations[input_size]);
    std::vector<double> right_predictions(&input_predictions[split_idx], &input_predictions[input_size]);
    std::vector<double> right_discordant(&input_discordant[split_idx], &input_discordant[input_size]);
    std::vector<double> right_pairs(&input_pairs[split_idx], &input_pairs[input_size]);
    
    left.push_back(left_observations);
    left.push_back(left_predictions);
    left.push_back(left_discordant);
    left.push_back(left_pairs);

    right.push_back(right_observations);
    right.push_back(right_predictions);
    right.push_back(right_discordant);
    right.push_back(right_pairs);
    

    left = merge_sort(left, outx);
    right = merge_sort(right, outx);
    output = merge_two_sides(left, right, outx);
    return output;
  }
}

// [[Rcpp::export]]
List merge_sort_c(std::vector<double> observations,
                  std::vector<double> predictions,
                  std::vector<double> discordant,
                  std::vector<double> pairs,
                  bool outx){

  std::list<std::vector<double> > input;
  std::list<std::vector<double> > output;

  List ret;

  input.push_back(observations);
  input.push_back(predictions);
  input.push_back(discordant);
  input.push_back(pairs);


  output = merge_sort(input, 1);

  ret["observations"] = output.front(); output.pop_front();
  ret["predictions"] = output.front(); output.pop_front();
  ret["discordant"] = output.front(); output.pop_front();
  ret["pairs"] = output.front(); output.pop_front();

  return ret;

}
// 
// int main(){
// 
//   std::list<std::vector<double> > input;
// 
//   std::list<std::vector<double> > output;
// 
//   double ci;
// 
//   std::vector<double> observations{1,2,3,4,6,5};
//   std::vector<double> predictions{1,2,3,4,5,6};
//   std::vector<double> discordant{0,0,0,0,0,0};
//   std::vector<double> pairs{5,5,5,5,5,5};
// 
//   std::vector<double> out_observations(observations.size());
//   std::vector<double> out_predictions(observations.size());
//   std::vector<double> out_discordant(observations.size());
//   std::vector<double> out_pairs(observations.size());
// 
//   input.push_back(observations);
//   input.push_back(predictions);
//   input.push_back(discordant);
//   input.push_back(pairs);
// 
//   output = merge_sort(input, 1);
// 
//   out_observations = output.front(); output.pop_front();
//   out_predictions = output.front(); output.pop_front();
//   out_discordant = output.front(); output.pop_front();
//   out_pairs = output.front(); output.pop_front();
// 
//   ci = std::accumulate(out_discordant.begin(), out_discordant.end(), 0.0)/2.0;
// 
//   std::cout << ci << " ";
// 
// }







