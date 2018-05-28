\name{findbeta}
\alias{findbeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The findbeta function
%%  ~~function to do ... ~~
}
\description{A function to estimate the parameters alpha and beta of a Beta distribution based on the existing prior beliefs (data and/or expert opinion). Information is provided about the mean (or the median or the mode) and whether it is lower or greater that a certain value with a pre-specified certainty (usually set at 95\%).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
findbeta(themean=NULL, themedian=NULL, themode=NULL,
percentile=0.95,lower.v=F, percentile.value)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{themean}{specify your prior beleif about the mean. It takes a value between 0 and 1. Not to be specified if a value has been given for the median or the mode.
%%     ~~Describe \code{x} here~~
}
  \item{themedian}{specify your prior beleif about the median. It takes a value between 0 and 1. Not to be specified if a value has been given for the mean or the mode.
%%     ~~Describe \code{x} here~~
}

  \item{themode}{specify your prior beleif about the mode. It takes a value between 0 and 1. Not to be specified if a value has been given for the mean or the median.
%%     ~~Describe \code{x} here~~
}

  \item{percentile}{specify the level of confidence that the true value of the mean (or the median or the mode) is greater or lower than the percentile.value. It takes a value between 0 and 1 and the default =0.95.
%%     ~~Describe \code{x} here~~
}

\item{lower.v}{logical, if TRUE the specified percentile.value is the upper limit for the mean (or the median or the mode) at the specified confidence level (percentile). If {FALSE} the specified percentile.value is the lower limit for the mean (or the median or the mode) at the specified confidence level (percentile). The default is FALSE.
%%     ~~Describe \code{x} here~~
}

\item{percentile.value}{specify the upper or lower limit for the mean (or the median or the mode) at the specified level of confidence (percentile). It takes a value between 0 and 1.
%%     ~~Describe \code{x} here~~
}
}

\references{
Branscum, A. J., Gardner, I. A., & Johnson, W. O. (2005): Estimation of diagnostic test sensitivity and specificity through Bayesian modeling. Preventive veterinary medicine, \bold{68}, 145--163.
%% ~put references to the literature/web site here ~
}

\examples{
##Example 1
##Based on the available literature the mean value for the sensitivity of a test
##is expected to be 0.90 and we can be 95\% sure that it is higher than 0.80.

findbeta(themean=0.90, percentile=0.95,lower.v=FALSE, percentile.value=0.80)

## The output is:
##[1] "The desired Beta distribution that satisfies the specified conditions is:
##Beta( 27.79 3.09 )"
##[1] "Here is a plot of the specified distribution."
##[1] "Descriptive statistics for this distribution are:"
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
## 0.6192  0.8688  0.9077  0.8996  0.9391  0.9974
##[1] "Verification: The percentile value 0.8 corresponds to the 0.05th
##percentile


##Example 2
##Based on the available literature the median value for the specificity of a
##test is expected to be 0.99 and we can be 95\% sure that it is higher than
##0.90.

findbeta(themedian=0.99, percentile=0.95,lower.v=FALSE, percentile.value=0.90)

## The output is:
##[1] "The desired Beta distribution that satisfies the specified conditions is:
##Beta( 18.97 0.52 )"
##[1] "Here is a plot of the specified distribution."
##[1] "Descriptive statistics for this distribution are:"
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
## 0.5634  0.9637  0.9871  0.9735  0.9970  1.0000
##[1] "Verification: The percentile value 0.9 corresponds to the 0.05th
##percentile"

##Example 3
##The most probable value (mode) for the prevalence of a disease/infection in a
##population is expected to be 0.15 and we are 90\% sure that it is less that
##0.40.

findbeta(themode=0.15, percentile=0.90,lower.v=TRUE, percentile.value=0.40)

## The output is:
##[1] "The desired Beta distribution that satisfies the specified conditions is:
##Beta( 2.15 7.52 )"
##[1] "Here is a plot of the specified distribution."
##[1] "Descriptive statistics for this distribution are:"
##    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
##0.001939 0.125475 0.204776 0.223724 0.300515 0.785156
##[1] "Verification: The percentile value 0.4 corresponds to the 0.9th
##percentile"




}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
