\name{nullCIDist}
\alias{nullCIDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
nullCIDist(n, multvect = c(1), makeplot = 0, outdir = "", cumulative = 1, return_polylist = 0, force_sym = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{multvect}{
%%     ~~Describe \code{multvect} here~~
}
  \item{makeplot}{
%%     ~~Describe \code{makeplot} here~~
}
  \item{outdir}{
%%     ~~Describe \code{outdir} here~~
}
  \item{cumulative}{
%%     ~~Describe \code{cumulative} here~~
}
  \item{return_polylist}{
%%     ~~Describe \code{return_polylist} here~~
}
  \item{force_sym}{
%%     ~~Describe \code{force_sym} here~~
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
structure(function (n, multvect = c(1), makeplot = 0, outdir = "", 
    cumulative = 1, return_polylist = 0, force_sym = 0) 
{
    mc <- match.call()
    encl <- parent.env(environment())
    called_args <- as.list(mc)[-1]
    default_args <- Filter(function(x) !identical(x, quote(expr = )), 
        as.list(formals()))
    default_args <- default_args[setdiff(names(default_args), 
        names(called_args))]
    args <- c(lapply(called_args, eval, parent.frame()), lapply(default_args, 
        eval, envir = environment()))
    hash <- encl$`_cache`$digest(c(as.character(body(encl$`_f`)), 
        args, lapply(encl$`_additional`, function(x) eval(x[[2L]], 
            environment(x)))))
    if (encl$`_cache`$has_key(hash)) {
        res <- encl$`_cache`$get(hash)
    }
    else {
        mc[[1L]] <- encl$`_f`
        res <- withVisible(eval(mc, parent.frame()))
        encl$`_cache`$set(hash, res)
    }
    if (res$visible) {
        res$value
    }
    else {
        invisible(res$value)
    }
  }, memoised = TRUE, class = c("memoised", "function"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
