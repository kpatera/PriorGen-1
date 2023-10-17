# load funs#### 
require(devtools)
require(remotes)
remotes::install_github("kpatera/PriorGen-1",force = T)
library(PriorGen)
#install.packages("spelling")
#install.packages("styler")
#install.packages("lintr")
#install.packages("covr")
require(covr)
require(spelling)
require(styler)
require(lintr)

sessioninfo::session_info(include_base = TRUE)

path="C:/Users/kosta/Documents/GitHub/PriorGen-1"
# Check example coverage + documentation ####
devtools::document()
#package_coverage(path = path, commentDonttest = FALSE, commentDontrun = FALSE)
package_coverage(path = path,type = c("examples", "vignettes"), commentDonttest = FALSE, commentDontrun = FALSE)

pkgdown::build_site(pkg = path)
# Check readme file rmarkdown #####
usethis::use_readme_rmd() # New readme file
#rmarkdown::render("C:/Users/kosta/Documents/GitHub/PriorGen-1/README.md", output_format = rmarkdown::github_document())

# Examples running ####
devtools::run_examples(run_dontrun = TRUE, run_donttest = TRUE)


# URL checker #####
urlchecker::url_check(path = path)

# Spelling mistakes #####
spelling::spell_check_package(pkg = path)

# Valid HTML5 for pages ####
roxygen2::roxygenise(package.dir = path)

# Check across multiple settings #####
rcmdcheck::rcmdcheck(path = path)

# Style of code  ####
styler::style_pkg(pkg = path)

# Detecting lints ####
lintr::lint_package(path)

# Exclude directories ####
usethis::use_build_ignore(c("local"))
usethis::use_pkgdown() 

# Check dependencies #####
dependsOnPkgs("rJava")
# Test ####
library(testthat)
library(PriorGen)

test_check("PriorGen")
# Build vignette #####
devtools::build_vignettes()
# Check commands ####
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

data.frame()
x_i1=c(5,	7,	8,	8,	8,	10,	10,	12,	13,	14,	14,	14,15)
d_i1=c(1,	1,	1,	0,	0,	0,	0,	1,	1,	0,	0,	0,	0)


x_i2=c(	3,	4,	5,	6,	7,	8,	10,	11, 12,	14,	15)
d_i2=c(1,	1,	1,	1,	1,	1,	0,	1,	1,	0,	0)
AG<-data.frame(x=c(x_i1,x_i2),d=c(d_i1,d_i2),group=c(rep(1,length(x_i1)),rep(2,length(x_i2))))

library(survival)
plot(survfit(data = AG,formula = Surv(time = AG$x,event = AG$d)~group))
# N(logit(0.2), 0.3 2)
# LO c + 𝜃 + 𝜏Z

plot(density(rnorm(n = 10000,mean = log(0.2/0.8),sd = 0.3^2)),xlim=c(-2,2))
lines(density(rnorm(n = 10000,mean = log(0.2/0.8)+1+0+0.05*rnorm(10000),sd = 0.3^2)))
