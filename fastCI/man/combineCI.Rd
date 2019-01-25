\name{combineCI}
\alias{combineCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
combineCI(CI1, CI2, N1, N2, multvect1, multvect2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CI1}{
%%     ~~Describe \code{CI1} here~~
}
  \item{CI2}{
%%     ~~Describe \code{CI2} here~~
}
  \item{N1}{
%%     ~~Describe \code{N1} here~~
}
  \item{N2}{
%%     ~~Describe \code{N2} here~~
}
  \item{multvect1}{
%%     ~~Describe \code{multvect1} here~~
}
  \item{multvect2}{
%%     ~~Describe \code{multvect2} here~~
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
function (CI1, CI2, N1, N2, multvect1, multvect2) 
{
    null1 <- nullCIDist(N1, multvect1, cumulative = FALSE)
    null2 <- nullCIDist(N2, multvect2, cumulative = FALSE)
    if (N1 < N2) {
        smaller <- null1
        bigger <- null2
    }
    else {
        smaller <- null2
        bigger <- null1
    }
    res.null <- convolve(bigger, smaller, type = "o")
    res.CI <- (choose(N1, 2) * CI1 + choose(N2, 2) * CI2)/(choose(N1, 
        2) + choose(N2, 2))
    res.pairs <- (choose(N1, 2) * CI1 + choose(N2, 2) * CI2)
    res.p <- getCIPvals(cumsum(res.null), res.pairs)
    return(c(res.CI, res.p))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
