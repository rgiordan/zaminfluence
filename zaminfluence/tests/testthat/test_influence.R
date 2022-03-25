#!/usr/bin/env Rscript
#
# Test the manual derivatives using numerical differentiation.
# Effectively, this tests GetIVSEDerivs and GetRegressionSEDerivs with
# both grouped and ungrouped standard errors.

library(AER)
library(zaminfluence)
library(testthat)
library(tidyverse)
library(purrr)

#context("zaminfluence")


# DELTEME
library(devtools)
devtools::load_all("/home/rgiordan/Documents/git_repos/zaminfluence/zaminfluence")


GenerateTestInstance <- function(do_iv, do_grouping) {
    x_dim <- 1
    param_true <- 0.1
    num_obs <- 500

    GenerateFun <- if (do_iv)
      GenerateIVRegressionData else GenerateRegressionData
    if (do_grouping) {
        df <- GenerateFun(num_obs, param_true, num_groups=5)
    } else {
        df <- GenerateFun(num_obs, param_true)
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
        fit_object <- ivreg(data=df, formula = reg_form,
                           x=TRUE, y=TRUE, weights=weights)
    } else {
        # Regression:
        x_names <- sprintf("x%d", 1:x_dim)
        reg_form <- formula(sprintf("y ~ %s - 1",
                                    paste(x_names, collapse=" + ")))
        fit_object <- lm(data=df, formula=reg_form,
                        x=TRUE, y=TRUE, weights=weights)
    }

    se_group <- if (do_grouping) df$se_group else NULL

    model_grads <-
        ComputeModelInfluence(fit_object) %>%
        AppendTargetRegressorInfluence("x1")
    signals <- GetInferenceSignals(model_grads)

    return(list(
        model_grads=model_grads,
        signals=signals,
        #reruns=reruns,
        fit_object=fit_object,
        se_group=se_group,
        df=df
    ))
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
    w_new <- GetWeightVector(drop_inds, orig_weights=model_grads$model_fit$weights)

    # The original values
    base_values <- GetBaseValues(param_infl)

    # Rerun
    rerun <- model_grads$RerunFun(w_new)
    rerun_base_values <- GetParameterInferenceQOIs(
      model_fit=rerun,
      target_parameter=param_infl$target_parameter,
      sig_num_ses=param_infl$sig_num_ses)
    diff_rerun <-
      unlist(rerun_base_values)[names(base_values)] %>% as.numeric() -
      base_values[names(base_values)]
    names(diff_rerun) <- names(base_values)

    # Prediction
    diff_pred <-
        map_dbl(names(base_values),
                ~ PredictChange(param_infl[[.]], drop_inds))
    names(diff_pred) <- names(base_values)
    rel_error <-
      (diff_pred - diff_rerun) /
      ifelse(abs(diff_rerun) > 1e-2, abs(diff_rerun), 1)

    # cat("\n\n",
    #   names(base_values), "\n",
    #   diff_pred, "\n",
    #   diff_rerun, "\n",
    #   rel_error, "\n"
    #   )
    # if (max(abs(rel_error)) > tol) {
    #   cat("\n***************************\n")
    # }

    # Check the maximum relative error amongst param, param_mzse, and param_pzse
    max_rel_err <- max(rel_error)
    expect_true(max_rel_err < tol,
                info=sprintf("%s %s prediction error: %f",
                             qoi_name, sign, max_rel_err))
}


# Test that the AMIP at a signal's AMIS produces the expected change in
# the target metric.
TestSignalPrediction <- function(param_infl, signals, signal_name) {
    signal <- signals[[signal_name]]
    qoi_name <- signal$qoi$name

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
  param_infl <- model_grads$param_infls[["x1"]]

  # Check the validity of the influence scores.
  qoi_names <- c("param", "param_mzse", "param_pzse")
  for (qoi_name in qoi_names) {
      validate_QOIInfluence(param_infl[[qoi_name]])
      for (sign in c("pos", "neg")) {
          # TestAPIP(param_infl[[qoi_name]], sign)
          TestPredictions(model_grads, param_infl, qoi_name, sign)
      }
  }

  # Check that the APIP predicts the appropriate change, and that
  # one fewer point does not.
  signals <- GetInferenceSignalsForParameter(param_infl)
  for (signal_name in c("sign", "sig", "both")) {
      TestSignalPrediction(param_infl, signals, signal_name)

      # Just test that this runs.
      #as.data.frame(signals[[signal_name]])
      as.data.frame(signals[[signal_name]])
  }

  # Test GetAMIS for a single QOI
  TestGetAMIS(param_infl$param_mzse)

  # Just that summary functions run
  signals <- GetInferenceSignals(model_grads)
  reruns <- RerunForSignals(signals, model_grads)
  preds <- PredictForSignals(signals, model_grads)

  reruns_df <- GetSignalsAndRerunsDataframe(signals, reruns, model_grads)
  preds_df <- GetSignalsAndRerunsDataframe(signals, preds, model_grads)

  PlotSignal(model_grads, signals, "x1", "sign", reruns=reruns, apip_max=0.03)

  return(invisible(test_instance))
}


#debug(TestPredictions)

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
