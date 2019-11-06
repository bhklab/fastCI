# library(Rcpp)

# sourceCpp("main.cpp", verbose=TRUE)


# dyn.load("main.so")

rCI <- function(x, y, dx, dy){
	dtx <- 0
	dty <- 0
	.Call("rCI", x, y, dx, dy, length(y), dtx, dty)
}