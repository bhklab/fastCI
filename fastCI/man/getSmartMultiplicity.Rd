\name{getSmartMultiplicity}
\alias{getSmartMultiplicity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getSmartMultiplicity(elements, multiplicity, norm = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{elements}{
%%     ~~Describe \code{elements} here~~
}
  \item{multiplicity}{
%%     ~~Describe \code{multiplicity} here~~
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
function (elements, multiplicity, norm = 1) 
{
    if (multiplicity == 1) {
        return(getSimplePolyProduct(elements, elements, norm = norm))
    }
    numlist <- getSimplePolyList(elements, multiplicity, norm = norm)
    for (ii in multiplicity:2) {
        x <- rep(1/ii, ii)
        for (jj in 1:length(numlist)) {
            if (sum(abs(as.numeric(as.polynomial(numlist[[jj]])\%\%x))) < 
                10^-14) {
                numlist <- c(numlist, list(as.numeric(divide.p(numlist[[jj]], 
                  x))))
                numlist[jj] <- c()
                break
            }
        }
    }
    if (sum(unlist(lapply(numlist, function(x) length(x) - 1))) != 
        (elements - multiplicity) * multiplicity) {
    }
    return(numlist)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
