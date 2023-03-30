
# PriorGen - Ver 2.0 beta

PriorGen — Generates Prior Distributions for Proportions

## Installation

To install the current source from GitHub use:

    install.packages(c("devtools", "remotes"))
    require(devtools)
    require(remotes)
    remotes::install_github("kpatera/PriorGen-1")

To install a stable version from the drat repository use:

    ## Will be added once a stable version is available

## Basic functions of PriorGen package

    require(PriorGen)

To run elicitation analysis on an example dataset:

    out<-findbeta(themode=0.15, percentile=0.90,lower.v=TRUE, percentile.value=0.40, silent=TRUE)
    out

To create a plot of the analysed data:

    findbeta_plot(out)

The basic two functions of PriorGen are findbeta() and findbetamupsi().
For help on these functions type:

    ?findbeta
    ?findbetamupsi
    ?findbeta_plot

## Troubleshooting

In case an error during download occur try the following

    remotes::install_github("kpateras/PriorGen-1", force = TRUE, dependecies = TRUE)

# Updates from Version 1 -\> Version 2.0 Beta

## Summary of updated and candidated functions

-   Updated basic findbeta (original function with percentiles)
    -   Mean, Median, Mode
-   Raw findbeta (standard location and scale measures
    (mean,median,mode,variance,range))
    -   Mean, Median, Mode
-   Abstract findbeta (general statements of how low/high are the mean
    and variance)
    -   Mean
-   Panel findbeta (Multiple experts or sources of information
    contribute to define this prior, based on simple averaging)
    -   Mean, Median, Mode
-   Updated basic findbetaqq (original functions with 2 percentiles)
    -   Percentiles
-   Updated basic findbetamupsi (original function with percentiles -
    two non-linear root functions)
    -   Mean
-   Raw findbetamupsi (standard location and scale measures (mean,
    variance))
    -   Mean
-   Abstract findbetamupsi (general statements of how low/high are the
    mean and variance)
    -   Mean
-   Panel findbeta (Multiple experts or sources of information
    contribute to define this prior, based on simple averaging)
    -   Mean
-   Plot all the above samples function
    -   Generic and applicable to all functions above.
-   New root and optim functions
    -   RootSolve alterantive
    -   Optim alternative @KP
