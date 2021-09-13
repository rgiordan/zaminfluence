#!/usr/bin/env Rscript
#
# Test the manual derivatives using numerical differentiation.
# Effectively, this tests GetIVSEDerivs and GetRegressionSEDerivs with
# both grouped and ungrouped standard errors.

library(AER)
library(zaminfluence)
#library(numDeriv)
#library(sandwich)
library(testthat)
library(tidyverse)
library(purrr)

context("zaminfluence")


GenerateTestInstance <- function(do_iv, do_grouping) {
    x_dim <- 1
    beta_true <- 0.1
    num_obs <- 500

    GenerateFun <- if (do_iv)
      GenerateIVRegressionData else GenerateRegressionData
    if (do_grouping) {
        df <- GenerateFun(num_obs, beta_true, num_groups=5)
    } else {
        df <- GenerateFun(num_obs, beta_true)
    }

    df$weights <- runif(nrow(df)) + 1

    # Fit a model.
    if (do_iv) {
        # IV:
        x_names <- sprintf("x%d", 1:x_dim)
        z_names <- sprintf("z%d", 1:x_dim)
        reg_form <- formula(sprintf("y ~ %s - 1 | %s - 1",
                                    paste(x_names, collapse=" + "),
                                    paste(z_names, collapse=" + ")))
        model_fit <- ivreg(data=df, formula = reg_form,
                           x=TRUE, y=TRUE, weights=weights)
    } else {
        # Regression:
        x_names <- sprintf("x%d", 1:x_dim)
        reg_form <- formula(sprintf("y ~ %s - 1",
                                    paste(x_names, collapse=" + ")))
        model_fit <- lm(data=df, formula=reg_form,
                        x=TRUE, y=TRUE, weights=weights)
    }

    se_group <- if (do_grouping) df$se_group else NULL

    model_grads <-
        ComputeModelInfluence(model_fit) %>%
        AppendTargetRegressorInfluence("x1")
    signals <-
        GetInferenceSignals(model_grads$param_infl_list[["x1"]]) %>%
        RerunForTargetChanges(model_grads)

    return(list(
        model_grads=model_grads,
        signals=signals,
        model_fit=model_fit,
        se_group=se_group,
        df=df
    ))
}



# Sanity check an APIP.
# - Influence scores should match the sign
# - Cumulative influence scores should decreas or increase according to the sign
# - Lengths should match
TestAPIP <- function(qoi,  sign) {
    apip <- qoi[[sign]]
    apip_infl <- qoi$infl[apip$infl_inds]
    if (sign == "pos") {
        expect_true(all(apip_infl > 0),
          info="positive influence is positive")
        expect_true(all(diff(apip$infl_cumsum) > 0),
          info="positive influence is increasing")
    } else if (sign == "neg") {
        expect_true(all(apip_infl < 0), info="negative influence is negative")
        expect_true(all(diff(apip$infl_cumsum) < 0),
          info="negative influence is decreasing")
    } else {
        stop("Bad sign passed to test")
    }
    expect_true(length(apip$infl_inds) == length(apip$infl_cumsum),
      info="Inds len == cumsum len")
}


# Check the validity of a prediction when leaving out a small number of
# influential points.
# We check that the relative error in the difference is less than 100 * tol %.
TestPredictions <- function(
      model_grads, param_infl, qoi_name, sign, num_leave_out=2, tol=0.07) {
    # Get the indices to drop.
    qoi <- param_infl[[qoi_name]]
    apip <- qoi[[sign]]
    drop_inds <- apip$infl_inds[1:num_leave_out]
    w_bool <- GetWeightVector(drop_inds, num_obs=model_grads$n_obs, bool=TRUE)

    # The original values
    base_values <- GetBaseValues(param_infl)

    # Rerun
    rerun <- model_grads$RerunFun(model_grads$model_fit, w_bool)
    rerun_base_values <- GetRerunBaseValues(rerun, param_infl)
    diff_rerun <- rerun_base_values - base_values

    # Prediction
    diff_pred <-
        map_dbl(names(base_values),
                ~ PredictChange(param_infl[[.]], drop_inds))

    # Check the maximum relative error amongst beta, beta_mzse, and beta_pzse
    max_rel_err <- max(abs((diff_pred - diff_rerun) / diff_rerun))
    expect_true(max_rel_err < tol,
                info=sprintf("%s %s prediction error: %f",
                             qoi_name, sign, max_rel_err))
}


# Test that the AMIP at a signal's AMIS produces the expected change in
# the target metric.
TestSignalPrediction <- function(param_infl, signals, signal_name) {
    signal <- signals[[signal_name]]
    qoi_name <- signal$qoi_name

    # Form the prediction
    drop_inds <- signal$apip$inds
    if (is.na(drop_inds[1])) {
      # There is nothing to test.
      return()
    }
    base_value <- GetBaseValues(param_infl)[qoi_name]
    pred_diff <- PredictChange(param_infl[[qoi_name]], drop_inds)
    pred_value <- base_value + pred_diff

    # Assert that a sign change took place.
    # Everything we look for is a sign change.
    # cat(signal_name, ": ", base_value, " ", pred_value, "\n")
    expect_true(
        sign(base_value) != sign(pred_value),
        sprintf("%s predicted to achieve sign change",
                signal_name))

    if (length(drop_inds) > 1) {
        pred_diff <- PredictChange(param_infl[[qoi_name]],
                                   drop_inds[1:(length(drop_inds) - 1)])
        pred_value <- base_value + pred_diff
        expect_true(
            sign(base_value) == sign(base_value + pred_diff),
            sprintf("%s predicted to fail to sign change with one fewer point",
                    signal_name))

    }
}


# Basic sanity checks on the number of points returned by GetAMIS
TestGetAMIS <- function(qoi) {
  n_drops <- c(0, 0.5, 1, 10, 10000)
  for (sign in c("pos", "neg")) {
      for (n_drop in n_drops) {
          suppressWarnings(amis <- GetAMIS(qoi, sign=sign, n_drop=n_drop))
          if (n_drop == 0) {
              expect_equivalent(amis, NULL)
          } else if (n_drop == 0.5) {
              expect_equivalent(amis, qoi[[sign]]$infl_inds[1])
          } else if (n_drop == 1) {
              expect_equivalent(amis, qoi[[sign]]$infl_inds[1])
          } else if (n_drop == 10) {
              expect_equivalent(amis, qoi[[sign]]$infl_inds[1:10])
          } else if (n_drop == 10000) {
              expect_equivalent(amis, qoi[[sign]]$infl_inds)
          } else {
              stop(sprintf("Bad number of drops: %d", n_drop))
          }
      }
  }
}


# Run all tests on a particular test instance.
TestInfluence <- function(test_instance) {
  model_fit <- test_instance$model_fit
  model_grads <- test_instance$model_grads
  param_infl <- model_grads$param_infl_list[["x1"]]

  # Check the validity of the influence scores.
  qoi_names <- c("beta", "beta_mzse", "beta_pzse")
  for (qoi_name in qoi_names) {
      for (sign in c("pos", "neg")) {
          TestAPIP(param_infl[[qoi_name]], sign)
          TestPredictions(model_grads, param_infl, qoi_name, sign)
      }
  }

  # Check that the APIP predicts the appropriate change, and that
  # one fewer point does not.
  signals <- GetInferenceSignals(param_infl)
  for (signal_name in c("sign", "sig", "both")) {
      TestSignalPrediction(param_infl, signals, signal_name)

      # Just test that this runs.
      GetSignalDataFrame(signals[[signal_name]])
  }

  # Test GetAMIS for a single QOI
  TestGetAMIS(param_infl$beta_mzse)
}


test_that("influence_computations_correct", {
  set.seed(42)
  for (do_iv in c(TRUE, FALSE)) {
    for (do_grouping in c(TRUE, FALSE)) {
      sprintf("Running %s %s",
          if (do_iv) "IV" else "OLS",
          if (do_grouping) "grouped" else "ungrouped", "\n") %>% cat()
      GenerateTestInstance(do_iv, do_grouping) %>%
        TestInfluence()
    }
  }
})
