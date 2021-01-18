#!/usr/bin/env Rscript

library(AER)
library(zaminfluence)
library(sandwich)
library(testthat)
library(tidyverse)

context("zaminfluence")

AssertNearlyEqual <- function(x, y, tol=1e-9) {
  diff_norm <- max(abs(x - y))
  info_str <- sprintf("%e > %e", diff_norm, tol)
  expect_true(diff_norm < tol, info=info_str)
}

GetSandwichCov <- function(fit, se_group=NULL) {
  if (is.null(se_group)) {
    return(vcov(fit))
  } else {
    return(vcovCL(fit, cluster=se_group, type="HC0", cadjust=FALSE))
  }
}


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

    AssertNearlyEqual(reg_zam$se_mat, GetSandwichCov(reg_res, se_group=se_group))
    AssertNearlyEqual(iv_zam$se_mat, GetSandwichCov(iv_res, se_group=se_group))
  }

  TestSEGroup(NULL)
  ordered_groups <- rep(1:20, each=5)
  TestSEGroup(ordered_groups)
  TestSEGroup(ordered_groups + 50)
  TestSEGroup((ordered_groups + 2) * 2)
  TestSEGroup(ordered_groups[sample(num_obs, replace=TRUE)])
  TestSEGroup(ordered_groups[sample(num_obs)])
})



# test_that("rerun works", {
num_obs <- 100
df <- GenerateIVRegressionData(num_obs, 0.5, num_groups=10)
df$w <- runif(num_obs) + 0.5
w_bool <- rep(TRUE, num_obs)
w_bool[sample(100, 10)] <- FALSE
new_w <- df$w
new_w[!w_bool] <- 0

iv_fit <- ivreg(y ~ x1 + 1 | z1 + 1, data=df, x=TRUE, y=TRUE, weights=df$w)
reg_fit <- lm(y ~ x1 + 1, data=df, x=TRUE, y=TRUE, weights=df$w)

for (use_iv in c(TRUE, FALSE)) {
  for (use_se_group in c(TRUE, FALSE)) {
    if (use_se_group) {
      se_group <- df$se_group
    } else {
      se_group <- NULL
    }
    if (use_iv) {
      new_fit <- ivreg(y ~ x1 + 1 | z1 + 1, data=df, x=TRUE, y=TRUE, weights=new_w)
      zam_fit <- RerunIVRegression(w_bool, iv_fit, se_group=se_group)
    } else {
      new_fit <- lm(y ~ x1 + 1, data=df, x=TRUE, y=TRUE, weights=new_w)
      zam_fit <- RerunRegression(w_bool, reg_fit, se_group=se_group)
    }
    new_vcov <- GetSandwichCov(new_fit, se_group=se_group)

    new_vcov
    zam_fit$se_cov

    cat(use_iv, use_se_group, new_vcov - zam_fit$se_mat, "\n\n")      
    print(new_vcov)
    print(zam_fit$se_mat)
    cat('------------------\n\n\n')
    # AssertNearlyEqual(new_fit$coefficients, zam_fit$betahat)
    # AssertNearlyEqual(new_vcov, zam_fit$se_cov)
  }
}
#})
