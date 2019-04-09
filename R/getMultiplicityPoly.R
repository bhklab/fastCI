getMultiplicityPoly <-
function(elements, multiplicity){
  num <- getSimplePolyProduct(elements, multiplicity)
  denom <- getSimplePolyProduct(multiplicity, multiplicity)
  
  return(as.numeric(as.polynomial(num)/denom))
}
