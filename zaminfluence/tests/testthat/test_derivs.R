#!/usr/bin/env Rscript

# Numerically test the derivatives.

library(AER)
library(zaminfluence)
library(testthat)
library(numDeriv)

context("zaminfluence")

test_that("derivatives work", {
  TestDerivs <- function(model_grads) {
    RerunCoeff <- function(w) {
      rerun_fit <- model_grads$RerunFun(w)
      return(rerun_fit$param)
    }

    RerunSe <- function(w) {
      rerun_fit <- model_grads$RerunFun(w)
      return(rerun_fit$se)
    }

    betahat_grad <- numDeriv::jacobian(
      RerunCoeff, model_grads$model_fit$weights)
    se_grad <- numDeriv::jacobian(RerunSe, model_grads$model_fit$weights)

    for (par in model_grads$parameter_names) {
      fit_ind <- GetParameterIndex(model_grads$model_fit, par)
      grad_ind <- GetParameterIndex(model_grads, par)
      AssertNearlyEqual(
        model_grads$param_grad[grad_ind, ], betahat_grad[fit_ind, ])
      AssertNearlyEqual(model_grads$se_grad[grad_ind, ], se_grad[fit_ind, ])
    }
  }

  TestRegressionConfigurationDerivs <- function(
        num_groups, weights, keep_pars, do_iv) {
    if (do_iv) {
      df <- GenerateIVRegressionData(
        num_obs, c(0.5, -0.5, 0.0), num_groups=num_groups)
      fit_object <- ivreg(y ~ x1 + x2 + x3 + 1 | z1 + z2 + z3 + 1,
                      data=df, x=TRUE, y=TRUE, weights=weights)
    } else {
      df <- GenerateRegressionData(
        num_obs, c(0.5, -0.5, 0.0), num_groups=num_groups)
      fit_object <-
        lm(y ~ x1 + x2 + x3 + 1, df, x=TRUE, y=TRUE, weights=weights)
    }
    model_grads <-
      ComputeModelInfluence(
        fit_object, se_group=df[["se_group"]], keep_pars=keep_pars)
    TestDerivs(model_grads)
  }

  test_configs <- expand.grid(
    num_groups=c(10, -1),
    random_weights=c(TRUE, FALSE),
    do_iv=c(TRUE, FALSE),
    keep_pars=c(1, 2, 3)
  )

  num_obs <- 20

  w_rand <- runif(num_obs) + 1
  w_ones <- rep(1, num_obs)

  for (n in 1:nrow(test_configs)) {
    config <- test_configs[n, ]
    cat(" Numerically testing derivatives for config ",
        paste(as.character(config), collapse=", "), "\n")
    weights <- if (config$random_weights) w_rand else w_ones
    keep_pars <-
      if (config$keep_pars == 1) {
        c("x1")
      } else if (config$keep_pars == 2) {
        c("x2", "x1")
      } else {
        NULL
      }
    num_groups <- if (config$num_groups == -1) NULL else config$num_groups
    TestRegressionConfigurationDerivs(
      num_groups, weights, keep_pars, config$do_iv)
  }
})
