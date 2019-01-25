\name{getSimplePolyProduct}
\alias{getSimplePolyProduct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getSimplePolyProduct(elts, range, norm = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{elts}{
%%     ~~Describe \code{elts} here~~
}
  \item{range}{
%%     ~~Describe \code{range} here~~
}
  \item{norm}{
%%     ~~Describe \code{norm} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (elts, range, norm = 0) 
{
    if (range > elts) {
        return("Error: Range (or number of polynomials) is greater than the input number of elements")
    }
    if (range < 1 | range != round(range)) {
        return("Error: Range must be a positive integer")
    }
    retpoly <- c(rep(1, (elts - range + 1)), rep(0, range/2 * 
        (2 * elts - range - 1) + 1 - (elts - range + 1)))
    count <- (elts - range + 1)
    if (range > 1) {
        for (k in (elts - range + 2):elts) {
            retpoly[1:(count + k)] = (cumsum(c(retpoly[1:count], 
                rep(0, k))) - cumsum(c(rep(0, k), retpoly[1:count])))
            count <- count + k - 1
        }
    }
    if (norm != 0) {
        retpoly <- retpoly/sum(retpoly)
    }
    return(retpoly[1:count])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
