#!/usr/bin/env Rscript

# Test the manually computed OLS and IV solutions in ols_iv_grads_lib.R

library(AER)
library(zaminfluence)
library(sandwich)
library(testthat)
library(tidyverse)

context("zaminfluence")


vcovWrap <- function(obj, cluster=NULL) {
  vcovCL(obj, cluster=cluster, type="HC0", cadjust=FALSE)
}


GetFitCovariance <- function(fit, se_group=NULL) {
  # Get a version of the sandwich covariance that should match
  # our computations.
  if (is.null(se_group)) {
    return(vcov(fit))
  } else {
    return(vcovCL(fit, cluster=se_group, type="HC0", cadjust=FALSE))
  }
}



# reg_infl should be the output of ComputeRegressionInfluence
TestConfiguration <- function(model_fit, se_group) {
  reg_infl <- ComputeModelInfluence(model_fit, se_group)

  # grad_df <- GetTargetRegressorGrads(reg_infl, "x1")
  # influence_dfs <- SortAndAccumulate(grad_df)
  #
  # Check that the regressions match.
  reg <- tidy(reg_infl$model_fit)
  testthat::expect_equivalent(reg$estimate, reg_infl$betahat)

  # Check that the target index is correct..
  # target_index <- attr(grad_df, "target_index")
  # base_vals <- attr(grad_df, "base_vals")
  # testthat::expect_equivalent(reg$estimate[target_index], base_vals["beta"])

  # Check that the standard errors match...
  if (is.null(se_group)) {
      se <- reg$std.error[target_index]
  } else {
      vcov_se_cov <- vcovCL(model_fit, cluster=se_group,
                            type="HC0", cadjust=FALSE)
      se <- sqrt(diag(vcov_se_cov))[target_index]
  }
  testthat::expect_equivalent(se, base_vals["se"])

  # Check the "scale" of the influence function is correct (note that this
  # may not match the actual standard errors for grouping).
  n_obs <- length(model_fit$y)
  vcov_se_cov <- vcovCL(model_fit, cluster=1:n_obs, type="HC0", cadjust=FALSE)
  robust_se <- sqrt(diag(vcov_se_cov))[target_index]

  # infl_scale <- GetInfluenceScale(grad_df$beta_grad)
  # testthat::expect_equivalent(sqrt(n_obs) * robust_se, infl_scale)

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
    AssertNearlyEqual(reg_zam$betahat, reg_res$coefficients)
    AssertNearlyEqual(iv_zam$betahat, iv_res$coefficients)

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


test_that("rerun works", {
  num_obs <- 100
  df <- GenerateIVRegressionData(num_obs, 0.5, num_groups=10)
  df$w <- runif(num_obs) + 0.5

  iv_fit <- ivreg(y ~ x1 + 1 | z1 + 1, data=df, x=TRUE, y=TRUE, weights=df$w)
  reg_fit <- lm(y ~ x1 + 1, data=df, x=TRUE, y=TRUE, weights=df$w)

  zam_iv_fit <- RerunIVRegression(rep(TRUE, num_obs), iv_fit, se_group=df$se_group)
  iv_vcov <- GetFitCovariance(iv_fit, se_group=df$se_group)
  AssertNearlyEqual(iv_fit$coefficients, zam_iv_fit$betahat)
  AssertNearlyEqual(iv_vcov, zam_iv_fit$se_mat)

  zam_reg_fit <- RerunRegression(rep(TRUE, num_obs), reg_fit, se_group=df$se_group)
  reg_vcov <- GetFitCovariance(reg_fit, se_group=df$se_group)
  AssertNearlyEqual(reg_fit$coefficients, zam_reg_fit$betahat)
  AssertNearlyEqual(reg_vcov, zam_reg_fit$se_mat)


  w_bool <- rep(TRUE, num_obs)
  w_bool[sample(100, 10)] <- FALSE
  new_w <- df$w
  # Note that this will fail if the weights are exactly zero,
  # since vcovCL is actually discontinuous when weights are set
  # to exactly zero.
  new_w[!w_bool] <- 1e-6
  # Make sure all the groups are still present
  stopifnot(length(unique(df$se_group[w_bool])) == length(unique(df$se_group)))
  se_group <- NULL

  zam_reg_fit <- ComputeRegressionResults(reg_fit, new_w, se_group=se_group)
  new_reg_fit <- lm(y ~ x1 + 1, data=df %>% mutate(w=!!new_w), weights=w)
  new_reg_vcov <- GetFitCovariance(new_reg_fit, se_group=se_group)
  AssertNearlyEqual(new_reg_fit$coefficients, zam_reg_fit$betahat, tol=1e-6)
  AssertNearlyEqual(new_reg_vcov, zam_reg_fit$se_mat)


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

      AssertNearlyEqual(new_fit$coefficients, zam_fit$betahat)
      AssertNearlyEqual(new_vcov, zam_fit$se_mat)
    }
  }
})
