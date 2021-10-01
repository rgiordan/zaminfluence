
library(AER)
library(zaminfluence)
library(testthat)
library(tidyverse)
library(purrr)

context("zaminfluence")



test_that("GetWeightVector_correct", {
  num_obs <- 23
  w <- runif(num_obs) + 0.5
  drop_inds <- c(1, 2, 4)

  wdrop <- w
  wdrop[drop_inds] <- 0

  # Test when the base weight is specified
  wtest <- GetWeightVector(drop_inds, orig_weights=w)
  AssertNearlyEqual(wtest, wdrop)

  wtest <- GetWeightVector(drop_inds, orig_weights=w, num_obs=num_obs)
  AssertNearlyEqual(wtest, wdrop)

  wtest <- GetWeightVector(drop_inds, orig_weights=w, bool=TRUE)
  AssertNearlyEqual(wtest, wdrop != 0)

  wtest <- GetWeightVector(drop_inds, orig_weights=w, bool=TRUE, invert=TRUE)
  AssertNearlyEqual(wtest, wdrop == 0)

  wkeep <- rep(0, num_obs)
  wkeep[drop_inds] <- w[drop_inds]
  wtest <- GetWeightVector(drop_inds, orig_weights=w, invert=TRUE)
  AssertNearlyEqual(wtest, wkeep)

  # Test when the base weight is not specified
  w <- rep(1, num_obs)
  wdrop <- w
  wdrop[drop_inds] <- 0

  wtest <- GetWeightVector(drop_inds, num_obs=num_obs)
  AssertNearlyEqual(wtest, wdrop)

  wtest <- GetWeightVector(drop_inds, num_obs=num_obs, bool=TRUE)
  AssertNearlyEqual(wtest, wdrop != 0)

  wtest <- GetWeightVector(drop_inds, num_obs=num_obs, bool=TRUE, invert=TRUE)
  AssertNearlyEqual(wtest, wdrop == 0)

  wkeep <- rep(0, num_obs)
  wkeep[drop_inds] <- w[drop_inds]
  wtest <- GetWeightVector(drop_inds, num_obs=num_obs, invert=TRUE)
  AssertNearlyEqual(wtest, wkeep)
})
