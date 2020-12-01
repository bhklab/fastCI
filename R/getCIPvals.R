getCIPvals <- function(nullCIDist, discpairs, alternative=c("two.sided", "greater", "less")){
  if (nullCIDist[length(nullCIDist)] < 0.99) {
    return("Error: nullCIDist is not a CDF.  Run nullCIDist with argument cumulative=1")
  } else {
    alternative <- match.arg(alternative)
    ix_twosided <- min(discpairs, length(nullCIDist) - 1 - discpairs)
    
    switch (alternative, 
            "two.sided" = ifelse(ix_twosided == (length(nullCIDist)-1)/2, 
                                 1, 
                                 2*nullCIDist[ix_twosided+1]),
            "greater" = 1 - nullCIDist[discpairs+1],
            "less" = nullCIDist[discpairs+1]
    )
  }
}
