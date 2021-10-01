#!/usr/bin/env Rscript

library(devtools)
devtools::load_all()

library(testthat)
library(zaminfluence)

#test_check("zaminfluence", reporter="summary")
test_file("testthat/test_derivs.R")
test_file("testthat/test_base_values.R")
test_file("testthat/test_influence.R")
test_file("testthat/test_catchall.R")
