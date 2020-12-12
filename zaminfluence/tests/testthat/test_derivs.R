#!/usr/bin/env Rscript

library(AER)
library(zaminfluence)
library(numDeriv)
library(sandwich)
library(testthat)
library(tidyverse)


context("zaminfluence")
source("utils.R")

AssertNearlyZero <- function(x, tol=1e-15) {
  if (max(abs(x)) > tol) {
      stop(sprintf("Diff %0.12e > tol = %0.12e", max(abs(x)), tol))
  }
}

test_that("derivs_correct", {
  # Generate data.
  x_dim <- 3
  beta_true <- runif(x_dim)
  df <- GenerateIVRegressionData(10, beta_true, num_groups=NULL)
  df$weights <- runif(nrow(df)) + 1

  # Fit an IV model.
  x_names <- sprintf("x%d", 1:x_dim)
  z_names <- sprintf("z%d", 1:x_dim)
  iv_form <- formula(sprintf("y ~ %s - 1 | %s - 1",
                             paste(x_names, collapse=" + "),
                             paste(z_names, collapse=" + ")))
  iv_fit <- ivreg(data=df, formula=iv_form, x=TRUE, y=TRUE, weights=weights)

  # Get influence.
  iv_infl <- ComputeModelInfluence(iv_fit)
  grad_df <- GetTargetRegressorGrads(iv_infl, "x1")
  influence_dfs <- SortAndAccumulate(grad_df)

  LocalGetIVSEDerivs <- function(w=df$weights, beta=iv_fit$coefficients) {
    GetIVSEDerivs(
      x=df[, x_names] %>% as.matrix(),
      z=df[, z_names] %>% as.matrix(),
      y=df$y,
      beta=beta,
      w0=w, testing=TRUE)
  }

  iv_se_list <- LocalGetIVSEDerivs()
  num_obs <- nrow(df)
  AssertNearlyZero(iv_se_list$betahat - iv_fit$coefficients, tol=1e-11)
  AssertNearlyZero(iv_se_list$sig2_hat - iv_fit$sigma^2)
  AssertNearlyZero(iv_se_list$se_mat - vcov(iv_fit), tol=1e-11)
  AssertNearlyZero(iv_se_list$se - vcov(iv_fit) %>% diag() %>% sqrt(), tol=1e-11)

  ########################
  # Test the derivatives

  # Test the partial derivatives
  w0 <- df$weights
  beta <- iv_fit$coefficients
  dsig2_hat_dw_num <-
      numDeriv::jacobian(function(w) {
        LocalGetIVSEDerivs(w=w)$sig2_hat
      }, w0) %>%
      as.numeric()
  AssertNearlyZero(dsig2_hat_dw_num - iv_se_list$dsig2_hat_dw_partial, tol=1e-8)

  dsand_mat_diag_dw_num <-
      numDeriv::jacobian(function(w) {
          LocalGetIVSEDerivs(w=w)$sand_mat %>% diag()
        }, w0)
  AssertNearlyZero(dsand_mat_diag_dw_num - iv_se_list$dsand_mat_diag_dw_partial, tol=1e-6)

  dse_mat_diag_dw_num <-
      numDeriv::jacobian(function(w) {
          LocalGetIVSEDerivs(w=w)$se_mat %>% diag()
        }, w0)
  AssertNearlyZero(dse_mat_diag_dw_num - iv_se_list$dse_mat_diag_dw_partial, tol=5e-6)

  dsig2_hat_dbeta_num <-
      numDeriv::jacobian(function(beta) {
          LocalGetIVSEDerivs(beta=beta)$sig2_hat
        }, beta)
  AssertNearlyZero(dsig2_hat_dbeta_num - iv_se_list$dsig2_hat_dbeta, tol=1e-8)

  #############################
  # Test the full derivatives
  GetIVTestResults <- function(w) {
      df_test <- df
      df_test$weights <- w
      beta_test <- ivreg(data=df_test, formula=iv_form,
                         x=TRUE, y=TRUE, weights=weights)$coefficients
      iv_se_test_list <- LocalGetIVSEDerivs(beta=beta_test, w=w)
      return(
          list(beta=iv_se_test_list$betahat,
               se=iv_se_test_list$se,
               sig2_hat=iv_se_test_list$sig2_hat,
               se_mat_diag=diag(iv_se_test_list$se_mat))
      )
  }

  dbetahat_dw_num <-
      numDeriv::jacobian(function(w) { GetIVTestResults(w)$beta }, w0)
  AssertNearlyZero(dbetahat_dw_num - iv_se_list$dbetahat_dw, tol=1e-8)
  #plot(dbetahat_dw_num, iv_se_list$dbetahat_dw); abline(0, 1)

  dsig2_hat_num <-
      numDeriv::jacobian(function(w) { GetIVTestResults(w)$sig2_hat }, w0)
  #plot(dsig2_hat_num, iv_se_list$dsig2_hat_dw); abline(0, 1)
  AssertNearlyZero(dsig2_hat_num - iv_se_list$dsig2_hat_dw, tol=1e-6)


  dse_mat_diag_dw_num <-
      numDeriv::jacobian(function(w) { GetIVTestResults(w)$se_mat_diag }, w0)
  #plot(dse_mat_diag_dw_num, iv_se_list$dse_mat_diag_dw); abline(0, 1)

  # Check the relative error for this one
  AssertNearlyZero((dse_mat_diag_dw_num - iv_se_list$dse_mat_diag_dw) /
                       (iv_se_list$dse_mat_diag_dw + 1e-3), tol=1e-6)

  dse_dw_num <- numDeriv::jacobian(function(w) { GetIVTestResults(w)$se }, w0)
  #plot(dse_dw_num, iv_se_list$dse_dw); abline(0, 1)
  # Check the relative error for this one
  AssertNearlyZero((dse_dw_num - iv_se_list$dse_dw) / (iv_se_list$dse_dw + 1e-3), tol=1e-6)


})
