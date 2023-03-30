require(devtools)
require(remotes)
remotes::install_github("kpatera/PriorGen-1",force = T)
library(PriorGen)

# Check example coverage + documentation ####
install.packages("covr")
require(covr)
devtools::document()
package_coverage(path = "C:/Users/kosta/Documents/GitHub/PriorGen-1", commentDonttest = FALSE, commentDontrun = FALSE)
package_coverage(path = "C:/Users/kosta/Documents/GitHub/PriorGen-1",type = c("examples", "vignettes"), commentDonttest = FALSE, commentDontrun = FALSE)


# Check readme file rmarkdown #####
usethis::use_readme_rmd() # New readme file
rmarkdown::render("C:/Users/kosta/Documents/GitHub/PriorGen-1/README.md", output_format = rmarkdown::github_document())

# Examples running ####
devtools::run_examples(run_dontrun = TRUE, run_donttest = TRUE)

# Test ####
library(testthat)
library(PriorGen)

test_check("PriorGen")


out<-findbeta(themode=0.15, percentile=0.90,lower.v=TRUE, percentile.value=0.40)
out2<-findbeta(themode=0.20, percentile=0.90,lower.v=TRUE, percentile.value=0.40)
plot(out,type="l")
print(out)

plot(out, type="l",lwd=3,ylim=c(0,3.7))
lines(out2, type="l",lwd=3,col="red")

temp=findbetamupsi(themean=0.10, percentile=0.79, 
                   lower.v=TRUE, percentile.value=0.26, 
                   psi.percentile=0.95, percentile.median=0.28, 
                   percentile95value=0.3)
temp$param_gamma
temp$param_beta

# N(logit(0.2), 0.3 2)
# LO c + 𝜃 + 𝜏Z

plot(density(rnorm(n = 10000,mean = log(0.2/0.8),sd = 0.3^2)),xlim=c(-2,2))
lines(density(rnorm(n = 10000,mean = log(0.2/0.8)+1+0+0.05*rnorm(10000),sd = 0.3^2)))
