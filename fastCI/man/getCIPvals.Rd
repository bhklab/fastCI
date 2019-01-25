\name{getCIPvals}
\alias{getCIPvals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getCIPvals(nullCIDist, discpairs, alternative = c("two.sided", "greater", "less"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nullCIDist}{
%%     ~~Describe \code{nullCIDist} here~~
}
  \item{discpairs}{
%%     ~~Describe \code{discpairs} here~~
}
  \item{alternative}{
%%     ~~Describe \code{alternative} here~~
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
function (nullCIDist, discpairs, alternative = c("two.sided", 
    "greater", "less")) 
{
    if (nullCIDist[length(nullCIDist)] < 0.99) {
        return("Error: nullCIDist is not a CDF.  Run nullCIDist with argument cumulative=1")
    }
    else {
        alternative <- match.arg(alternative)
        ix_twosided <- min(discpairs, length(nullCIDist) - discpairs)
        switch(alternative, two.sided = ifelse(discpairs == length(nullCIDist)/2, 
            1, 2 * nullCIDist[ix_twosided + 1]), greater = 1 - 
            nullCIDist[discpairs + 1], less = nullCIDist[discpairs + 
            1])
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
