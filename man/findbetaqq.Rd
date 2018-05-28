\name{findbetaqq}
\alias{findbetaqq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The findbetaqq function
%%  ~~function to do ... ~~
}
\description{A function to estimate the parameters alpha and beta of a Beta distribution based on the existing prior belief (data and/or expert opinion) about the values of two distinct percentiles.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
findbetaqq(percentile.value1,percentile1,percentile.value2,percentile2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{percentile.value1}{specify the value for the first percentile. It takes a value between 0 and 1.
%%     ~~Describe \code{x} here~~
}
  \item{percentile1}{specify which is the percentile that corrersponds to percentile.value1. It takes a value between 0 and 1.
%%     ~~Describe \code{x} here~~
}

\item{percentile.value2}{specify the value for the second percentile. It takes a value between 0 and 1.
%%     ~~Describe \code{x} here~~
}
  \item{percentile2}{specify which is the percentile that corrersponds to percentile.value2. It takes a value between 0 and 1.
%%     ~~Describe \code{x} here~~
}
}

\references{
Kostoulas, P., Nielsen, S. S., Branscum, A. J., Johnson, W. O., Dendukuri, N., Dhand, N. K., Toft, N., Gardner, I. A. (2017): Reporting guidelines for diagnostic accuracy studies that use Bayesian latent class models (STARD--BLCM). Statistics in medicine, \bold{23}, 3603--3604.
%% ~put references to the literature/web site here ~
}

\examples{
##Example 1
##We believe that 20\% of the units in an area/region have a prevalence of
##disease/infection less than or equal to 0.30 while at the same time we are 90\%
##certain that the prevalence is less than 0.60

findbetaqq(percentile.value1=0.30,percentile1=0.20,percentile.value2=0.60,percentile2=0.90)

##The output is:
##[1] "The desired Beta distribution that satisfies the specified conditions is:
##Beta( 5.19 7.17 )"
##[1] "Here is a plot of the specified distribution."
##[1] "Descriptive statistics for this distribution are:"
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##0.02999 0.31999 0.41419 0.41974 0.51413 0.85487
##[1] "Verification: The first percentile value 0.3 corresponds to the 0.2 th
##percentile"
##[1] "Verification: The second percentile value 0.6 corresponds to the 0.9 th
##percentile"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
