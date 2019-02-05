library(devtools)
library(roxygen2)

#/Users/matt/OneDrive/regsim
setwd("/Users/matt/OneDrive/")
# initialze and create package
#usethis::create_package("regsim")

setwd("/Users/matt/OneDrive/regsim")
#document()
setwd("..")
install("regsim")
library(regsim)

setwd("/Users/matt/OneDrive/regsim")
check()

# library(devtools)
# library(roxygen2)
# library(testthat)
#
# setwd("/Users/matt/OneDrive/giftedCalcs")
# document()
#
# setwd("..")
# install("giftedCalcs")
#
# setwd("/Users/matt/OneDrive/giftedCalcs")
# check()
# devtools::test()
#
# test(filter="var_mean")
# test(filter="r_identified")
# test(filter="d_identified")
# test(filter="d_identified_v")
# test(filter="q_identified")
# test(filter="p_identified")
# test(filter="mean_identified")
# test(filter="sd_identified")
# test(filter="marginal_psychometrics")
# test(filter="conditional_moments")
# test(filter="conditional_p_id")
# test(filter="estimate_valid")
# test(filter="boot_estimate_valid")
# test(filter="estimate_performance")
#
#
#
#
#
#
#
#
#
#
# setwd("..")
# install("giftedCalcs")
#
#
# detach("package:giftedCalcs", unload=T)
# remove.packages("giftedCalcs")
#
# devtools::install_github("mcbeem/giftedCalcs", force=T)
# library(giftedCalcs)
# help(package="giftedCalcs")

