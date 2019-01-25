\name{plotNormalApproxN}
\alias{plotNormalApproxN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotNormalApproxN(N)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
%%     ~~Describe \code{N} here~~
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
function (N) 
{
    plot.new()
    test <- nullCIDist(N, 1, cumulative = 0, force_sym = 1)
    tt <- naive_approximate_normal_from_density(test)
    tt2 <- computeExpectedApproximation(N)
    plot(0:(length(test) - 1), test, type = "l", col = "red")
    lines(0:(length(test)), dnorm(0:length(test), tt[1], tt[2]), 
        col = "blue")
    lines(0:(length(test)), dnorm(0:length(test), tt2[1], tt2[2]), 
        col = "green")
    legend("topleft", legend = c("True", "Fitted Approximation", 
        "Expected Approximation"), fill = c("red", "blue", "green"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
