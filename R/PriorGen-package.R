#' Generates Prior Distributions for Proportions
#'
#' Translates beliefs into prior information in the form
#' of Beta and Gamma distributions. It can be mainly used
#' for the generation of priors on the prevalence of disease
#' and the sensitivity/specificity of diagnostic tests.
#'
#' @docType package
#' @name PriorGen
#' @import graphics
#' @import rootSolve
#' @import nleqslv
#' @importFrom stats dbeta
#' @importFrom stats optim
#' @importFrom stats qbeta
#' @importFrom stats qgamma
#' @importFrom stats qnorm
#' @importFrom stats rbeta
#' @importFrom stats rgamma
#' @importFrom stats runif
#' @importFrom stats var
NULL
