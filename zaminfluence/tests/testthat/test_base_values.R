#!/usr/bin/env Rscript

# Test the manually computed OLS and IV solutions in ols_iv_grads_lib.R

library(AER)
library(zaminfluence)
library(sandwich)
library(testthat)
library(tidyverse)

context("zaminfluence")


# Test that ComputeModelInfluence, AppendTargetRegressorInfluence, and
# RerunFun give the same answers as R on the original data.
TestConfiguration <- function(model_fit, se_group) {
  model_grads <-
    ComputeModelInfluence(model_fit, se_group)%>%
    AppendTargetRegressorInfluence("x1")

  # Test that the coefficient estimates and standard errors in model_grads
  # match what we expect from R.
  AssertNearlyEqual(
    model_grads$model_fit$paramhat, coefficients(model_fit), desc="paramhat equal")
  se_r <- GetFitCovariance(model_fit, se_group) %>% diag() %>% sqrt()
  AssertNearlyEqual(
    model_grads$model_fit$se, se_r, desc="std error equal")
  testthat::expect_equivalent(
    model_grads$model_fit$n_obs, length(model_fit$y), info="num obs")
  # testthat::expect_equivalent(
  #   model_grads$weights, model_fit$weights, info="weights")
  testthat::expect_equivalent(
    model_grads$model_fit$parameter_names, names(coefficients(model_fit)),
    info="column names")

  # Test that the base values in param_infl are correct.
  param_infl <- model_grads$param_infls[["x1"]]
  target_index <- param_infl$target_index
  testthat::expect_equivalent(
    "x1", names(coefficients(model_fit))[target_index], info="target index")
  testthat::expect_equivalent(
    param_infl$param$base_value,
    coefficients(model_fit)[target_index],
    info="param base value")
  testthat::expect_equivalent(
    param_infl$param_pzse$base_value,
    coefficients(model_fit)[target_index] +
      param_infl$sig_num_ses * se_r[target_index],
    info="param_pzse base value")
  testthat::expect_equivalent(
    param_infl$param_mzse$base_value,
    coefficients(model_fit)[target_index] -
      param_infl$sig_num_ses * se_r[target_index],
    info="param_mzse base value")

  # Test that if we re-run we get the same answer.
  rerun <- model_grads$RerunFun(model_fit$weights)
  AssertNearlyEqual(
    rerun$paramhat, coefficients(model_fit), desc="rerun paramhat equal")
  AssertNearlyEqual(
    rerun$se, se_r, desc="rerun std error equal")
}



test_that("regression works", {
  TestRegressionConfiguration <- function(num_groups, weights) {
    df <- GenerateRegressionData(100, 0.5, num_groups=num_groups)
    lm_result <- lm(y ~ x1 + 1, df, x=TRUE, y=TRUE, weights=weights)
    TestConfiguration(lm_result, se_group=df[["se_group"]])
  }

  TestRegressionConfiguration(num_groups=NULL, weights=NULL)
  TestRegressionConfiguration(num_groups=10, weights=NULL)
  TestRegressionConfiguration(num_groups=NULL, weights=runif(100))
  TestRegressionConfiguration(num_groups=10, weights=runif(100))

  TestIVRegressionConfiguration <- function(num_groups, weights) {
    df <- GenerateIVRegressionData(100, 0.5, num_groups=num_groups)
    iv_res <- ivreg(y ~ x1 + 1 | z1 + 1,
                    data=df, x=TRUE, y=TRUE, weights=weights)
    TestConfiguration(iv_res, se_group=df[["se_group"]])
  }

  TestIVRegressionConfiguration(num_groups=NULL, weights=NULL)
  TestIVRegressionConfiguration(num_groups=10, weights=NULL)
  TestIVRegressionConfiguration(num_groups=NULL, weights=runif(100))
  TestIVRegressionConfiguration(num_groups=10, weights=runif(100))

})


test_that("se groups can be non-ordered", {
  num_obs <- 100
  df <- GenerateIVRegressionData(num_obs, 0.5, num_groups=10)
  iv_res <- ivreg(y ~ x1 + 1 | z1 + 1, data=df, x=TRUE, y=TRUE)
  reg_res <- lm(y ~ x1 + 1, data=df, x=TRUE, y=TRUE)

  TestSEGroup <- function(se_group) {
    reg_zam <- ComputeRegressionResults(reg_res, se_group=se_group)
    iv_zam <- ComputeIVRegressionResults(iv_res, se_group=se_group)

    # The coefficients shouldn't depend on se_group, but test for good measure
    AssertNearlyEqual(reg_zam$paramhat, reg_res$coefficients)
    AssertNearlyEqual(iv_zam$paramhat, iv_res$coefficients)

    AssertNearlyEqual(reg_zam$se_mat, GetFitCovariance(reg_res, se_group=se_group))
    AssertNearlyEqual(iv_zam$se_mat, GetFitCovariance(iv_res, se_group=se_group))
  }

  TestSEGroup(NULL)
  ordered_groups <- rep(1:20, each=5)
  TestSEGroup(ordered_groups)
  TestSEGroup(ordered_groups + 50)
  TestSEGroup((ordered_groups + 2) * 2)
  TestSEGroup(ordered_groups[sample(num_obs, replace=TRUE)])
  TestSEGroup(ordered_groups[sample(num_obs)])
})


# Check that rerun matches R with left-out observations.
test_that("rerun works", {
  # Generate base data.
  num_obs <- 100
  df <- GenerateIVRegressionData(num_obs, 0.5, num_groups=10)
  df$w <- runif(num_obs) + 0.5

  iv_fit <- ivreg(y ~ x1 + 1 | z1 + 1, data=df, x=TRUE, y=TRUE, weights=df$w)
  reg_fit <- lm(y ~ x1 + 1, data=df, x=TRUE, y=TRUE, weights=df$w)

  # Use Rerun to get fits using our code.  Check that our results
  # match R's results.  (Note that this is only an extra sanity check
  # here --- this is principally tested above in TestConfiguration)
  zam_iv_fit <- ComputeIVRegressionResults(
    iv_fit, weights=df$w, se_group=df$se_group)
  iv_vcov <- GetFitCovariance(iv_fit, se_group=df$se_group)
  AssertNearlyEqual(iv_fit$coefficients, zam_iv_fit$paramhat)
  AssertNearlyEqual(iv_vcov, zam_iv_fit$se_mat)

  zam_reg_fit <- ComputeRegressionResults(
    reg_fit, weights=df$w, se_group=df$se_group)
  reg_vcov <- GetFitCovariance(reg_fit, se_group=df$se_group)
  AssertNearlyEqual(reg_fit$coefficients, zam_reg_fit$paramhat)
  AssertNearlyEqual(reg_vcov, zam_reg_fit$se_mat)

  # Test that rerun works with left-out observations.  Generate a weight
  # vector with randomly left-out observations.
  w_bool <- rep(TRUE, num_obs)
  w_bool[sample(100, 10)] <- FALSE
  new_w <- df$w

  # Note that this test will fail if the weights are exactly zero,
  # since vcovCL is actually discontinuous when weights are set
  # to exactly zero.  To avoid this, make the weights very small instead.
  new_w[!w_bool] <- 1e-6

  # Make sure all the groups are still present
  stopifnot(length(unique(df$se_group[w_bool])) == length(unique(df$se_group)))

  # Re-run using OLS or IV, and grouped standard errors or not, and check
  # that our results match R.
  for (use_iv in c(TRUE, FALSE)) {
    for (use_se_group in c(TRUE, FALSE)) {
      if (use_se_group) {
        se_group <- df$se_group
      } else {
        se_group <- NULL
      }
      if (use_iv) {
        new_fit <-
          ivreg(y ~ x1 + 1 | z1 + 1, data=df %>%
            mutate(w=!!new_w), x=TRUE, y=TRUE, weights=w)
        zam_fit <- ComputeIVRegressionResults(iv_fit, new_w, se_group=se_group)
      } else {
        new_fit <-
          lm(y ~ x1 + 1, data=df %>%
            mutate(w=!!new_w), x=TRUE, y=TRUE, weights=w)
        zam_fit <- ComputeRegressionResults(reg_fit, new_w, se_group=se_group)
      }
      new_vcov <- GetFitCovariance(new_fit, se_group=se_group)

      AssertNearlyEqual(new_fit$coefficients, zam_fit$paramhat)
      AssertNearlyEqual(new_vcov, zam_fit$se_mat)
    }
  }
})
