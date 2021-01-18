#!/usr/bin/env Rscript

library(AER)
library(zaminfluence)
library(sandwich)
library(testthat)
library(tidyverse)

context("zaminfluence")


test_that("se groups can be non-ordered", {
  df <- GenerateIVRegressionData(100, 0.5, num_groups=10)
  iv_res <- ivreg(y ~ x1 + 1 | z1 + 1, data=df, x=TRUE, y=TRUE)
  reg_res <- lm(y ~ x1 + 1, data=df, x=TRUE, y=TRUE)

  ComputeRegressionResults(reg_res)
  ComputeIVRegressionResults(iv_res)
})
