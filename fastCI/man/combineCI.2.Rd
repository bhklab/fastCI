\name{combineCI.2}
\alias{combineCI.2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
combineCI.2(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function (x) 
{
    x <- x[, complete.cases(t(x))]
    total.N <- sum(x[2, ])
    if (total.N >= 100) {
        vec1 <- x[, 1]
        total.out <- choose(vec1[2], 2)
        conc.out <- vec1[1] * total.out
        pars.out <- computeExpectedApproximation(total.out)
        x <- x[, -1, drop = FALSE]
        while (NCOL(x)) {
            vec1 <- x[, 1]
            pars.new <- computeExpectedApproximation(vec1[2])
            pars.out[1] <- pars.out[1] + pars.new[1]
            pars.out[2] <- sqrt(pars.out[2]^2 + pars.new[2]^2)
            conc.out <- conc.out + vec1[1] * choose(vec1[2], 
                2)
            total.out <- total.out + choose(vec1[2], 2)
            x <- x[, -1, drop = FALSE]
        }
        if (conc.out < pars.out[1]) {
            prob <- pnorm(conc.out + 1, mean = pars.out[1], sd = pars.out[2])
            p.out <- 2 * prob
        }
        else {
            prob <- pnorm(conc.out - 1, mean = pars.out[1], sd = pars.out[2], 
                lower.tail = FALSE)
            p.out <- 2 * prob
        }
        CI.out <- conc.out/total.out
    }
    else {
        nullTable <- makeTableUpToN(max(x[2, ]))
        vec1 <- x[, 1]
        x <- x[, -1, drop = FALSE]
        nullOut <- nullTable[[vec1[2]]]
        conc.out <- vec1[1] * choose(vec1[2], 2)
        total.out <- choose(vec1[2], 2)
        while (NCOL(x)) {
            vec1 <- x[, 1]
            nullOut <- convolve(nullOut, nullTable[[vec1[2]]], 
                type = "o")
            total.out <- total.out + choose(vec1[2], 2)
            conc.out <- conc.out + vec1[1] * choose(vec1[2], 
                2)
            x <- x[, -1, drop = FALSE]
        }
        CI.out <- conc.out/total.out
        p.out <- getCIPvals(cumsum(nullOut), total.out - conc.out)
    }
    return(c(CI = CI.out[[1]], p = p.out[[1]]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
