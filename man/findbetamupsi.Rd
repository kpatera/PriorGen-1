\name{findbetamupsi}
\alias{findbetamupsi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The findbetamupsi function
%%  ~~function to do ... ~~
}
\description{A function to estimate (a) the  parameters of a Beta distribution for the expected mean of a proportion - usually the prevalence of disease/infection for the units in an area/region and (b) the parameters of a Gamma distribution expressing our prior beleif about the variability of the prevalence estimates across the units of the area/region under consideration.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
findbetamupsi(themean, percentile=0.95, lower.v=T,
percentile.value,  psi.percentile=0.90, percentile.median,
percentile95value)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{themean}{specify your prior beleif about the mean. It takes a value between 0 and 1.
%%     ~~Describe \code{x} here~~
}

  \item{percentile}{specify the level of confidence that the true value of the mean is greater or lower than the percentile.value. It takes a value between 0 and 1 and the default is 0.95.
%%     ~~Describe \code{x} here~~
}

  \item{lower.v}{logical, if TRUE the specified percentile.value is the upper limit for the mean at the specified confidence level (percentile). If {FALSE} the specified percentile.value is the lower limit for the mean at the specified confedence level (percentile).The default is TRUE.
%%     ~~Describe \code{x} here~~
}

\item{percentile.value}{specify the upper or lower limit for the mean at the specified level of confidence (percentile). It takes a value between 0 and 1.
%%     ~~Describe \code{x} here~~
}

  \item{psi.percentile}{specify the level of confidence that a certain fraction of the units under study has a prevalence less than the percentile.median. It takes a value between 0 and 1 and the default is 0.90.
%%     ~~Describe \code{x} here~~
}

  \item{percentile.median}{specify the median value that corresponds to the defined psi.percentile. It takes a value between 0 and 1 and has to be higher than both themean and the percentile.
%%     ~~Describe \code{x} here~~
}

  \item{percentile95value}{specify the value that the percentile.median does not exceed with 95\% confidence. It takes a value between 0 and 1 and has to be higher than the percentile.median.
%%     ~~Describe \code{x} here~~
}

}

\references{
Branscum, A. J., Gardner, I. A., & Johnson, W. O. (2004): Bayesian modeling of animal-- and herd--level prevalences. Preventive Veterinary Medicine, \bold{66}, 101--112.
%% ~put references to the literature/web site here ~
}


\examples{
##Example 1
##The mean prevalence of a disease/infection for the units within an area/region
##is thought to be 0.20 and we are 99\% confident that it is not more than 0.40.
##Within this area, we are also confident that 90\% of all units have a prevalence
##less or equal to 0.50 and we are 95\% certain that it does not exceed 0.60

findbetamupsi(themean=0.20, percentile=0.99, lower.v=TRUE,
percentile.value=0.30, psi.percentile=0.90,
percentile.median=0.50, percentile95value=0.60)

##The output is:
##[1] "The desired Beta distribution that satisfies the specified conditions
##about the mean of the prevalence 'mu' is: Beta( 20.27 81.07 )"
##[1] "The desired Gamma distribution that satisfies the specified conditions
##about the variability 'psi' of the prevalence is: Gamma( 8.97 2.79 )"
##[1] "The plot gives the specified prior beleif on the prevalence distribution."
##[1] "Descriptive statistics for this distrubiton are:"
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##0.00000 0.03589 0.13164 0.20156 0.30799 0.99971
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
