
# PriorGen - Ver 1.2 beta

PriorGen — Generates Prior Distributions for Proportions

## Installation

To install the current source from GitHub use:

    install.packages(c("devtools", "remotes"))
    require(devtools)
    require(remotes)
    remotes::install_github("kpatera/PriorGen-1")

To install a stable version from the drat repository use:

    install.packages("PriorGen")

## Basic functions of PriorGen package

    require(PriorGen)

To run elicitation analysis on an example dataset:

    out<-findbeta(themode=0.15, percentile=0.90,lower.v=TRUE, percentile.value=0.40, silent=TRUE)
    out

To create a plot of the analyzed data as a PriorGen class object:

    plot(out)

The basic two functions of PriorGen are findbeta() and findbetamupsi().
For help on these functions type:

    ?findbeta
    ?findbetamupsi
    ?findbeta_panel

## Troubleshooting

In case an error during download occur try the following

    remotes::install_github("kpateras/PriorGen-1", force = TRUE, dependecies = TRUE)
@KP
