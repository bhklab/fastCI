pad <-
function(x, len, location = c("end", "start", "sym")){
  location = match.arg(location)
  if(length(x) > len){
    stop("Cannot do negative padding.")
  }
  if (location == "end"){
    return(c(x, numeric(len - length(x))))
  }
  if (location == "start"){
    return(c(numeric(len - length(x)), x))
  }
}
