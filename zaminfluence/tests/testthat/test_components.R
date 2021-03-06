#!/usr/bin/env Rscript

library(AER)
library(zaminfluence)
library(sandwich)
library(testthat)
library(tidyverse)



context("zaminfluence")

vcovWrap <- function(obj, cluster=NULL) {
  vcovCL(obj, cluster=cluster, type="HC0", cadjust=FALSE)
}

AssertNearlyEqual <- function(x, y, tol=1e-9) {
  diff_norm <- max(abs(x - y))
  info_str <- sprintf("%e > %e", diff_norm, tol)
  expect_true(diff_norm < tol, info=info_str)
}

GetSandwichCov <- function(fit, se_group=NULL) {
  # Get a version of the sandwich covariance that should match
  # our computations.
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


test_that("rerun works", {
  num_obs <- 100
  df <- GenerateIVRegressionData(num_obs, 0.5, num_groups=10)
  df$w <- runif(num_obs) + 0.5
  
  iv_fit <- ivreg(y ~ x1 + 1 | z1 + 1, data=df, x=TRUE, y=TRUE, weights=df$w)
  reg_fit <- lm(y ~ x1 + 1, data=df, x=TRUE, y=TRUE, weights=df$w)
  
  zam_iv_fit <- RerunIVRegression(rep(TRUE, num_obs), iv_fit, se_group=df$se_group)
  iv_vcov <- GetSandwichCov(iv_fit, se_group=df$se_group)
  AssertNearlyEqual(iv_fit$coefficients, zam_iv_fit$betahat)
  AssertNearlyEqual(iv_vcov, zam_iv_fit$se_mat)
  
  zam_reg_fit <- RerunRegression(rep(TRUE, num_obs), reg_fit, se_group=df$se_group)
  reg_vcov <- GetSandwichCov(reg_fit, se_group=df$se_group)
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
  new_reg_vcov <- GetSandwichCov(new_reg_fit, se_group=se_group)
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
        new_fit <- ivreg(y ~ x1 + 1 | z1 + 1, data=df %>% mutate(w=!!new_w), x=TRUE, y=TRUE, weights=w)
        zam_fit <- ComputeIVRegressionResults(iv_fit, new_w, se_group=se_group)
      } else {
        new_fit <- lm(y ~ x1 + 1, data=df %>% mutate(w=!!new_w), x=TRUE, y=TRUE, weights=w)
        zam_fit <- ComputeRegressionResults(reg_fit, new_w, se_group=se_group)
      }
      new_vcov <- GetSandwichCov(new_fit, se_group=se_group)
  
      AssertNearlyEqual(new_fit$coefficients, zam_fit$betahat)
      AssertNearlyEqual(new_vcov, zam_fit$se_mat)
    }
  }
})



if (FALSE) {
  ###########################################################
  # Some extra debugging, involving the fact that 
  # there is discontinuous behavior with the weights in vcovCL
  
  RegressWithWeights <- function(w, se_group=NULL) {
    df$w <- w
    reg_fit <- lm(y ~ x1 + 1, data=df, x=TRUE, y=TRUE, weights=df$w)
    reg_vcov <- GetSandwichCov(reg_fit, se_group=se_group)
    return(list(betahat=reg_fit$coefficients, se_mat=reg_vcov))
  }
  
  
  eps_vec <- c(0, 1e-6, 2e-6)
  foo <- rep(NA, length(eps_vec))
  for (i in 1:length(eps_vec)) {
    w_new <- df$w
    w_new[1] <- eps_vec[i]
    foo[i] <- RegressWithWeights(w_new, se_group=df$se_group)$se_mat[1, 1]
  }
  foo[2] - foo[1]
  foo[3] - foo[2]
  
  foo[2] / foo[1]
  (num_obs - 2) / (num_obs - 3)
  
  
  ########################
  
  library(devtools)
  devtools::load_all("/home/rgiordan/Documents/git_repos/zaminfluence/zaminfluence")
  
  num_obs <- 100
  df <- GenerateIVRegressionData(num_obs, 0.5, num_groups=10)
  df$w <- runif(num_obs) + 0.5
  #df$w[1:10] <- 1e-6
  w0 <- df$w
  w0[1:10] <- 1e-6  # This will fail with weights that are exactly zero due to nonlinearity in vcovCL.
  reg_fit <- lm(y ~ x1 + 1, data=df, x=TRUE, y=TRUE, weights=w)
  reg_fit0 <- lm(y ~ x1 + 1, data=df %>% mutate(w=!!w0), x=TRUE, y=TRUE, weights=w)
  se_group <- NULL
  
  GetRegressionVariables(reg_fit)
  
  zam_reg_fit <- ComputeRegressionResults(reg_fit, weights=df$w)
  zam_reg_fit0 <- ComputeRegressionResults(reg_fit, weights=w0)
  zam_reg_fit00 <- ComputeRegressionResults(reg_fit0, weights=w0)
  
  # zaminfluence matches itself
  AssertNearlyEqual(zam_reg_fit0$betahat, zam_reg_fit00$betahat)
  AssertNearlyEqual(zam_reg_fit0$se_mat, zam_reg_fit00$se_mat)
  
  reg_vcov <- GetSandwichCov(reg_fit, se_group=se_group)
  AssertNearlyEqual(reg_fit$coefficients, zam_reg_fit$betahat)
  AssertNearlyEqual(reg_vcov, zam_reg_fit$se_mat)
  
  max(abs(reg_vcov0 - reg_vcov)) # Check something actually changed
  reg_vcov0 <- GetSandwichCov(reg_fit0, se_group=se_group)
  AssertNearlyEqual(reg_vcov0, zam_reg_fit0$se_mat)
  AssertNearlyEqual(reg_vcov0, zam_reg_fit00$se_mat)
  
  
  ########################
  
  num_obs <- 100
  df <- GenerateIVRegressionData(num_obs, 0.5, num_groups=10)
  df$w <- runif(num_obs) + 0.5
  reg_fit <- lm(y ~ x1 + 1, data=df, x=TRUE, y=TRUE, weights=w)
  
  w0 <- w1 <- w2 <- df$w
  w0[1:10] <- 0
  w1[1:10] <- 1e-2
  w2[1:10] <- 2e-2
  
  reg0 <- lm(y ~ x1 + 1, data=df %>% mutate(w=!!w0), x=TRUE, y=TRUE, weights=w)
  reg1 <- lm(y ~ x1 + 1, data=df %>% mutate(w=!!w1), x=TRUE, y=TRUE, weights=w)
  reg2 <- lm(y ~ x1 + 1, data=df %>% mutate(w=!!w2), x=TRUE, y=TRUE, weights=w)
  weights(reg0)[1]
  weights(reg1)[1]
  weights(reg2)[1]
  
  zam_reg_fit0 <- ComputeRegressionResults(reg_fit, weights=w0)
  zam_reg_fit1 <- ComputeRegressionResults(reg_fit, weights=w1)
  zam_reg_fit2 <- ComputeRegressionResults(reg_fit, weights=w2)
  
  (coefficients(reg0) - zam_reg_fit0$betahat) %>% max() %>% abs()
  (coefficients(reg1) - zam_reg_fit1$betahat) %>% max() %>% abs()
  (coefficients(reg2) - zam_reg_fit2$betahat) %>% max() %>% abs()
  
  # All match
  GetSandwichCov(reg0)[1, 1] - zam_reg_fit0$se_mat[1, 1]
  GetSandwichCov(reg1)[1, 1] - zam_reg_fit1$se_mat[1, 1]
  GetSandwichCov(reg2)[1, 1] - zam_reg_fit2$se_mat[1, 1]
  
  # Why do the vcovCL results differ?
  
  # Here is part of the problem.
  vcovCL(reg0)[1, 1] - vcovCL(reg1)[1, 1]
  vcovCL(reg1)[1, 1] - vcovCL(reg2)[1, 1]
  
  meat(reg0)[1, 1] - meat(reg1)[1, 1]
  meat(reg1)[1, 1] - meat(reg2)[1, 1]
  
  max(abs(estfun(reg0) - estfun(reg1)))
  max(abs(estfun(reg1) - estfun(reg2)))
  
  # It's the bread
  
  bread(reg0)[1, 1] - bread(reg1)[1, 1]
  bread(reg1)[1, 1] - bread(reg2)[1, 1]
  
  nobs(reg0)
  nobs(reg1)
  nobs(reg2)
  
  vcov(reg0) / vcov(reg1)
  vcov(reg1) / vcov(reg2)
  
  
  ###################
}