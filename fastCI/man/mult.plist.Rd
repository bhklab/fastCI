\name{mult.plist}
\alias{mult.plist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mult.plist(plist, outOrder)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plist}{
%%     ~~Describe \code{plist} here~~
}
  \item{outOrder}{
%%     ~~Describe \code{outOrder} here~~
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
function (plist, outOrder) 
{
    if (missing(outOrder)) {
        outOrder <- sum(unlist(lapply(plist, length))) - length(plist) + 
            1
    }
    if (length(plist) == 1) {
        return(plist[[1]])
    }
    product <- fft(pad(plist[[1]], outOrder))
    for (ii in 2:length(plist)) {
        product <- product * fft(pad(plist[[ii]], outOrder))
    }
    return(fft(1/length(product) * product, inverse = TRUE))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
