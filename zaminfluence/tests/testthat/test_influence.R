#!/usr/bin/env Rscript
#
# Test the manual derivatives using numerical differentiation.
# Effectively, this tests GetIVSEDerivs and GetRegressionSEDerivs with
# both grouped and ungrouped standard errors.

library(AER)
library(zaminfluence)
library(numDeriv)
library(sandwich)
library(testthat)
library(tidyverse)

context("zaminfluence")


TestInfluence <- function(test_instance) {

}
