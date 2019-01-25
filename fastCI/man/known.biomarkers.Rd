\name{known.biomarkers}
\alias{known.biomarkers}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("known.biomarkers")}
\format{
  A data frame with 37 observations on the following 8 variables.
  \describe{
    \item{\code{Drug}}{a factor with levels \code{17-AAG} \code{AZD0530} \code{AZD6244} \code{Crizotinib} \code{Erlotinib} \code{lapatinib} \code{Nilotinib} \code{Nutlin-3} \code{PD-0325901} \code{PD-0332991} \code{PHA-665752} \code{PLX4720} \code{Sorafenib} \code{TAE684}}
    \item{\code{Gene}}{a factor with levels \code{ALK} \code{BCR_ABL} \code{BRAF} \code{CDK4} \code{CDK6} \code{EGFR} \code{ERBB2} \code{FLT1} \code{FLT3} \code{FLT4} \code{HGF} \code{HSP90AA1} \code{HSP90AB1} \code{KDR} \code{KIT} \code{MAP2K1} \code{MAP2K2} \code{MDM2} \code{MET} \code{NQO1} \code{PDGFRA} \code{RAF1} \code{SRC}}
    \item{\code{Type}}{a factor with levels \code{amplification} \code{expression} \code{fusion} \code{mutation}}
    \item{\code{GDSC.effect.size}}{a numeric vector}
    \item{\code{GDSC.pvalue}}{a numeric vector}
    \item{\code{CCLE.effect.size}}{a numeric vector}
    \item{\code{CCLE.pvalue}}{a numeric vector}
    \item{\code{Reproducible}}{a factor with levels \code{} \code{NO} \code{NS} \code{YES}}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(known.biomarkers)
## maybe str(known.biomarkers) ; plot(known.biomarkers) ...
}
\keyword{datasets}
