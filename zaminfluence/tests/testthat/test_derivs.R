#!/usr/bin/env Rscript
#
# Test the manual derivatives using numerical differentiation.

library(AER)
library(zaminfluence)
library(numDeriv)
library(sandwich)
library(testthat)
library(tidyverse)

context("zaminfluence")
source("utils.R")


AssertNearlyZero <- function(x, tol=1e-15) {
  x_norm <- max(abs(x))
  info_str <- sprintf("%e > %e", x_norm, tol)
  expect_true(x_norm < tol, info=info_str)
}


TestRegressionDerivatives <- function(do_iv) {
  # Generate data.
  x_dim <- 3
  beta_true <- runif(x_dim)

  if (do_iv) {
    df <- GenerateIVRegressionData(10, beta_true, num_groups=NULL)
  } else {
    df <- GenerateRegressionData(10, beta_true, num_groups=NULL)
  }
  df$weights <- runif(nrow(df)) + 1

  if (do_iv) {
    x_names <- sprintf("x%d", 1:x_dim)
    z_names <- sprintf("z%d", 1:x_dim)
    iv_form <- formula(sprintf("y ~ %s - 1 | %s - 1",
                               paste(x_names, collapse=" + "),
                               paste(z_names, collapse=" + ")))
    reg_fit <- ivreg(data=df, formula=iv_form, x=TRUE, y=TRUE, weights=weights)

    LocalGetRegressionSEDerivs <- function(w=df$weights,
                                           beta=reg_fit$coefficients) {
      GetIVSEDerivs(
        x=df[, x_names] %>% as.matrix(),
        z=df[, z_names] %>% as.matrix(),
        y=df$y,
        beta=beta,
        w0=w, testing=TRUE)
      }
    reg_se_list <- LocalGetRegressionSEDerivs()

  } else {
    x_names <- sprintf("x%d", 1:x_dim)
    reg_form <- formula(sprintf("y ~ %s - 1",
                               paste(x_names, collapse=" + ")))
    reg_fit <- lm(data=df, formula=reg_form, x=TRUE, y=TRUE, weights=weights)

    LocalGetRegressionSEDerivs <- function(w=df$weights,
                                           beta=reg_fit$coefficients) {
      GetRegressionSEDerivs(
          x=df[, x_names] %>% as.matrix(),
          y=df$y,
          beta=beta,
          w0=w,
          testing=TRUE)
      }
    reg_se_list <- LocalGetRegressionSEDerivs()
  }

  if (do_iv) {
    reg_sigma <- reg_fit$sigma
  } else {
    reg_sigma <- sigma(reg_fit)
  }
  AssertNearlyZero(reg_se_list$betahat - reg_fit$coefficients, tol=1e-11)
  AssertNearlyZero(reg_se_list$sig2_hat - reg_sigma^2, tol=1e-11)
  AssertNearlyZero(reg_se_list$se_mat - vcov(reg_fit), tol=1e-11)
  AssertNearlyZero(reg_se_list$se -
                   vcov(reg_fit) %>% diag() %>% sqrt(), tol=1e-11)

  ########################
  # Test the derivatives

  # Test the partial derivatives
  w0 <- df$weights
  beta <- reg_fit$coefficients
  dsig2_hat_dw_num <-
      numDeriv::jacobian(function(w) {
        LocalGetRegressionSEDerivs(w=w)$sig2_hat
      }, w0) %>%
      as.numeric()
  AssertNearlyZero(dsig2_hat_dw_num -
                   reg_se_list$dsig2_hat_dw_partial, tol=1e-8)

  dsand_mat_diag_dw_num <-
      numDeriv::jacobian(function(w) {
          LocalGetRegressionSEDerivs(w=w)$sand_mat %>% diag()
        }, w0)
  AssertNearlyZero(dsand_mat_diag_dw_num -
                   reg_se_list$dsand_mat_diag_dw_partial, tol=1e-6)

  dse_mat_diag_dw_num <-
      numDeriv::jacobian(function(w) {
          LocalGetRegressionSEDerivs(w=w)$se_mat %>% diag()
        }, w0)
  AssertNearlyZero(dse_mat_diag_dw_num -
                   reg_se_list$dse_mat_diag_dw_partial, tol=5e-6)

  dsig2_hat_dbeta_num <-
      numDeriv::jacobian(function(beta) {
          LocalGetRegressionSEDerivs(beta=beta)$sig2_hat
        }, beta)
  AssertNearlyZero(
    dsig2_hat_dbeta_num - reg_se_list$dsig2_hat_dbeta, tol=1e-7)

  #############################
  # Test the full derivatives
  if (do_iv) {
    GetRegTestResults <- function(w) {
        df_test <- df
        df_test$weights <- w
        beta_test <- ivreg(data=df_test, formula=iv_form,
                           x=TRUE, y=TRUE, weights=weights)$coefficients
        iv_se_test_list <- LocalGetRegressionSEDerivs(beta=beta_test, w=w)
        return(
            list(beta=iv_se_test_list$betahat,
                 se=iv_se_test_list$se,
                 sig2_hat=iv_se_test_list$sig2_hat,
                 se_mat_diag=diag(iv_se_test_list$se_mat))
        )
    }
  } else {
    GetRegTestResults <- function(w) {
        df_test <- df
        df_test$weights <- w
        beta_test <- lm(data=df_test, formula=reg_form,
                           x=TRUE, y=TRUE, weights=weights)$coefficients
        reg_se_test_list <- LocalGetRegressionSEDerivs(beta=beta_test, w=w)
        return(
            list(beta=reg_se_test_list$betahat,
                 se=reg_se_test_list$se,
                 sig2_hat=reg_se_test_list$sig2_hat,
                 se_mat_diag=diag(reg_se_test_list$se_mat))
        )
    }
  }

  dbetahat_dw_num <-
      numDeriv::jacobian(function(w) { GetRegTestResults(w)$beta }, w0)
  AssertNearlyZero(dbetahat_dw_num - reg_se_list$dbetahat_dw, tol=1e-8)

  dsig2_hat_num <-
      numDeriv::jacobian(function(w) { GetRegTestResults(w)$sig2_hat }, w0)
  AssertNearlyZero(dsig2_hat_num - reg_se_list$dsig2_hat_dw, tol=1e-6)


  dse_mat_diag_dw_num <-
      numDeriv::jacobian(function(w) { GetRegTestResults(w)$se_mat_diag }, w0)

  # Check the relative error for this one
  AssertNearlyZero((dse_mat_diag_dw_num - reg_se_list$dse_mat_diag_dw) /
                       (reg_se_list$dse_mat_diag_dw + 1e-3), tol=1e-6)

  dse_dw_num <- numDeriv::jacobian(function(w) {
    GetRegTestResults(w)$se }, w0)
  # Check the relative error for this one
  AssertNearlyZero((dse_dw_num - reg_se_list$dse_dw) /
                   (reg_se_list$dse_dw + 1e-3), tol=1e-6)

}

TestGroupedRegressionDerivatives <- function(do_iv) {
  x_dim <- 3
  beta_true <- runif(x_dim)

  if (do_iv) {
      df <- GenerateIVRegressionData(20, beta_true, num_groups=5)
  } else {
      df <- GenerateRegressionData(20, beta_true, num_groups=5)
  }
  df$weights <- runif(nrow(df)) + 1

  # Fit a model.
  if (do_iv) {
      # IV:
      x_names <- sprintf("x%d", 1:x_dim)
      z_names <- sprintf("z%d", 1:x_dim)
      iv_form <- formula(sprintf("y ~ %s - 1 | %s - 1",
                                 paste(x_names, collapse=" + "),
                                 paste(z_names, collapse=" + ")))
      reg_fit <- ivreg(data=df, formula = iv_form,
                       x=TRUE, y=TRUE, weights=weights)
      # reg_cov <- sandwich::vcovCL(
      #   reg_fit, cluster=df$se_group, type="HC0", cadjust=FALSE)
  } else {
      # Regression:
      x_names <- sprintf("x%d", 1:x_dim)
      reg_form <- formula(sprintf("y ~ %s - 1",
                                  paste(x_names, collapse=" + ")))
      reg_fit <- lm(data=df, formula=reg_form,
                    x=TRUE, y=TRUE, weights=weights)
      # reg_cov <- sandwich::vcovCL(
      #   reg_fit, cluster=df$se_group, type="HC0", cadjust=FALSE)
  }

  # Get influence.

  # if (do_iv) {
  #     reg_se_list <- GetIVSEDerivs(
  #         x=df[, x_names] %>% as.matrix(),
  #         z=df[, z_names] %>% as.matrix(),
  #         y=df$y,
  #         beta=reg_fit$coefficients,
  #         w0=df$weights,
  #         se_group=df$se_group,
  #         testing=TRUE)
  # } else {
  #     reg_se_list <- GetRegressionSEDerivs(
  #         x=df[, x_names] %>% as.matrix(),
  #         y=df$y,
  #         beta=reg_fit$coefficients,
  #         w0=df$weights,
  #         se_group=df$se_group,
  #         testing=TRUE)
  # }

  LocalGetRegressionSEDerivs <- function(w=df$weights,
                                         beta=reg_fit$coefficients) {
      if (do_iv) {
          return(GetIVSEDerivs(
              x=df[, x_names] %>% as.matrix(),
              z=df[, z_names] %>% as.matrix(),
              y=df$y,
              beta=beta,
              w0=w,
              se_group=df$se_group,
              testing=TRUE))
      } else {
          return(GetRegressionSEDerivs(
              x=df[, x_names] %>% as.matrix(),
              y=df$y,
              beta=beta,
              w0=w,
              se_group=df$se_group,
              testing=TRUE))
      }
  }

  reg_se_list <- LocalGetRegressionSEDerivs()

  w0 <- df$weights
  beta <- reg_fit$coefficients

  AssertNearlyZero(reg_se_list$betahat - reg_fit$coefficients, tol=1e-9)
  AssertNearlyZero(colMeans(reg_se_list$s_mat), tol=1e-9)
  num_groups <- max(df$se_group) + 1
  AssertNearlyZero(cov(reg_se_list$s_mat) * (num_groups - 1) / num_groups -
                   reg_se_list$v_mat, tol=1e-9)

  vcov_se_cov <- vcovCL(reg_fit, cluster=df$se_group, type="HC0", cadjust=FALSE)
  AssertNearlyZero(vcov_se_cov - reg_se_list$se_mat, tol=1e-9)
  AssertNearlyZero(sqrt(diag(vcov_se_cov)) - reg_se_list$se, tol=1e-9)

  # Test that the s_mat_expanded worked correctly.
  for (n in 1:length(df$se_group)) {
      gn <- df$se_group[n] + 1
      AssertNearlyZero(reg_se_list$s_mat_expanded[n, ] -
                       reg_se_list$s_mat[gn, ], tol=1e-9)
  }

  #######

  ddiag_semat_dw_partial_num <-
      numDeriv::jacobian(function(w) {
          LocalGetRegressionSEDerivs(w=w)$se_mat %>% diag()
      }, w0)
  AssertNearlyZero(ddiag_semat_dw_partial_num -
                   reg_se_list$ddiag_semat_dw_partial, tol=1e-9)

  ddiag_semat_dbeta_partial_num <-
      numDeriv::jacobian(function(beta) {
          LocalGetRegressionSEDerivs(beta=beta)$se_mat %>% diag()
      }, beta)
  AssertNearlyZero(ddiag_semat_dbeta_partial_num -
                   reg_se_list$ddiag_semat_dbeta_partial, tol=1e-9)

  #######################
  # Full derivatives
  GetRegTestResults <- function(w) {
      df_test <- df
      df_test$weights <- w
      beta_test <- LocalGetRegressionSEDerivs(w=w)$betahat
      reg_se_test_list <- LocalGetRegressionSEDerivs(beta=beta_test, w=w)
      return(
          list(beta=reg_se_test_list$betahat,
               se=reg_se_test_list$se,
               se_mat_diag=diag(reg_se_test_list$se_mat))
      )
  }

  dbetahat_dw_num <-
      numDeriv::jacobian(function(w) { GetRegTestResults(w)$beta }, w0)
  AssertNearlyZero(dbetahat_dw_num - reg_se_list$dbetahat_dw, tol=1e-9)

  dse_mat_diag_dw_num <-
      numDeriv::jacobian(function(w) { GetRegTestResults(w)$se_mat_diag }, w0)
  AssertNearlyZero(dse_mat_diag_dw_num - reg_se_list$dse_mat_diag_dw, tol=1e-9)

  dse_dw_num <-
      numDeriv::jacobian(function(w) { GetRegTestResults(w)$se }, w0)
  AssertNearlyZero(dse_dw_num - reg_se_list$dse_dw, tol=1e-9)
}


test_that("ungrouped_derivs_correct", {
  set.seed(12611)
  TestRegressionDerivatives(do_iv=TRUE)
  set.seed(12611)
  TestRegressionDerivatives(do_iv=FALSE)
})


test_that("grouped_derivs_correct", {
  set.seed(12611)
  TestGroupedRegressionDerivatives(do_iv=TRUE)
  set.seed(12611)
  TestGroupedRegressionDerivatives(do_iv=FALSE)
})
