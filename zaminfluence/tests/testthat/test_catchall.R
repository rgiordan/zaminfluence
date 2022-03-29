
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

  # Test failure
  expect_error(GetWeightVector(drop_inds, num_obs=num_obs - 1, orig_weights=w))
  expect_error(GetWeightVector(c(1, -4), num_obs=num_obs))
  expect_error(GetWeightVector(c(1, num_obs + 1), num_obs=num_obs))
  expect_error(GetWeightVector(c(1, 0), num_obs=num_obs))
  expect_error(GetWeightVector(drop_inds))
})


test_that("NAN_in_APIP_issue25", {
  # This is for Issue 25 on github
  # https://github.com/rgiordan/zaminfluence/issues/25
  #
  # In short, check that zaminfluence can gracefully handle cases where
  # the desired change cannot be produced.

  library(zaminfluence)
  library(testthat)

  set.seed(42)
  x_dim <- 1

  # We should not be able to reverse the sign of such a large parameter.
  param_true <- 1000
  df <- GenerateRegressionData(10, param_true, num_groups=NULL)

  # Fit a regression model.
  x_names <- sprintf("x%d", 1:x_dim)
  reg_form <- formula(sprintf("y ~ %s - 1", paste(x_names, collapse=" + ")))
  fit_object <- lm(data=df, formula=reg_form, x=TRUE, y=TRUE)

  # Get influence and reruns.
  # Derivatives are only computed for keep_pars, which can be in any order.
  model_grads <-
      ComputeModelInfluence(fit_object) %>%
      AppendTargetRegressorInfluence("x1")
  signals <- GetInferenceSignals(model_grads)

  expect_error()
  signals$x1$sign$apip
  testthat::expect_true(is.na(signals$x1$sign$apip$inds))
  testthat::expect_false(signals$x1$sign$apip$success)

  preds <- PredictForSignals(signals, model_grads)
  reruns <- RerunForSignals(signals, model_grads)
  testthat::expect_true(is.null(preds$x1$sign))
  testthat::expect_true(is.null(reruns$x1$sign))
})



test_that("keep_pars_works", {

set.seed(42)
x_dim <- 3
num_obs <- 10

param_true <- 1:x_dim / x_dim
df <- GenerateRegressionData(num_obs, param_true, num_groups=NULL)

# Fit a regression model.
x_names <- sprintf("x%d", 1:x_dim)
reg_form <- formula(sprintf("y ~ %s - 1", paste(x_names, collapse=" + ")))
fit_object <- lm(data=df, formula=reg_form, x=TRUE, y=TRUE)

# Check the behavior of keep_pars.
model_grads_full <- ComputeModelInfluence(
  fit_object, keep_pars=c("x1", "x2", "x3"))
model_grads_null <- ComputeModelInfluence(fit_object)
model_grads_2 <- ComputeModelInfluence(fit_object, keep_pars="x2")
model_grads_23 <- ComputeModelInfluence(fit_object, keep_pars=c("x2", "x3"))
model_grads_32 <- ComputeModelInfluence(fit_object, keep_pars=c("x3", "x2"))

TestKeepPars <- function(m, grad_pars, full_m=model_grads_full) {
    k <- length(grad_pars)

    # Make sure the kept parameter names are correct
    expect_equivalent(m$parameter_names, grad_pars)

    # The ModelFit should have the whole parameter vector
    expect_equivalent(m$model_fit$parameter_names, c("x1", "x2", "x3"))

    # The gradients should only be computed for the kept parameters
    # and should match the corresponding gradients in the full model.
    expect_equivalent(dim(m$param_grad), c(k, num_obs))
    expect_equivalent(dim(m$se_grad), c(k, num_obs))
    for (par in grad_pars) {
        i <- GetParameterIndex(m, par)
        full_i <- GetParameterIndex(full_m, par)
        expect_equivalent(m$param_grad[i, ], full_m$param_grad[full_i, ])
        expect_equivalent(m$se_grad[i, ], full_m$se_grad[full_i, ])
    }
}

TestKeepPars(model_grads_full, c("x1", "x2", "x3"))
TestKeepPars(model_grads_null, c("x1", "x2", "x3"))
TestKeepPars(model_grads_2, c("x2"))
TestKeepPars(model_grads_23, c("x2", "x3"))
TestKeepPars(model_grads_32, c("x3", "x2"))

# Check that it fails if you request a parameter that's not present
expect_error(
  AppendTargetRegressorInfluence(model_grads_23, "x4"), "x4 not found")
})
