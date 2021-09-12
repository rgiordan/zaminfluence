#!/usr/bin/env Rscript


# Test the manually computed derivatives in ols_iv_grads_lib.R

library(AER)
library(zaminfluence)
library(sandwich)
library(testthat)
library(tidyverse)

context("zaminfluence")


test_that("plotting works", {
  df <- GenerateRegressionData(100, 0.5)
  lm_result <- lm(y ~ x1 + 1, df, x=TRUE, y=TRUE)
  influence_dfs <-
    ComputeRegressionInfluence(lm_result) %>%
    GetTargetRegressorGrads("x1") %>%
    SortAndAccumulate()

  PlotInfluence(influence_dfs$sign, "num_removed", 10)

  target_change <- GetRegressionTargetChange(influence_dfs, "num_removed")
  PlotInfluence(influence_dfs$sign, "num_removed", 10, target_change)

  # Just so that there is a test statement and it's not described as skipped.
  testthat::expect_true(TRUE)
})


test_that("rerun works", {
  set.seed(42)

  test_rerun <- function(model_fit, paired) {
    if (paired) {
      grad_df <-
        ComputeModelInfluence(model_fit) %>%
        GetTargetRegressorGrads("x1") %>%
        bind_cols(df[c("arm", "group")])
      influence_dfs <-
          PairInfluenceScores(grad_df,
                              assignment_col="arm",
                              group_cols="group",
                              level0=0,
                              level1=1)
    } else {
      influence_dfs <-
          ComputeModelInfluence(model_fit) %>%
          GetTargetRegressorGrads("x1") %>%
          SortAndAccumulate()
    }

    target_change <- GetRegressionTargetChange(influence_dfs, "num_removed")
    base_vals <- SafeGetBaseVals(influence_dfs$sign)

    # Check the sign change
    sign_alpha_target <- filter(
      target_change, change == "sign")$num_removed
    direction <- filter(target_change, change == "sign")$direction
    infl_df <- influence_dfs$sign[[direction]]
    rerun_result <- RerunTargetModelForAlpha(
        infl_df, model_fit, "num_removed", sign_alpha_target)
    testthat::expect_true(sign(rerun_result$beta) != sign(base_vals["beta"]))

    # Check the significance change
    sig_alpha_target <- filter(
      target_change, change == "significance")$num_removed
    direction <- filter(target_change, change == "significance")$direction
    infl_df <- influence_dfs$sig[[direction]]
    rerun_result <- RerunTargetModelForAlpha(
        infl_df, model_fit, "num_removed", sig_alpha_target)
    rerun_sig <- sign(rerun_result$beta_pzse) == sign(rerun_result$beta_mzse)
    base_sig <- sign(base_vals["beta_pzse"]) == sign(base_vals["beta_mzse"])
    testthat::expect_true(rerun_sig != base_sig)

    # Check the sign and significance change
    sign_sig_alpha_target <- filter(
      target_change, change == "sign and significance")$num_removed
    direction <- filter(
      target_change, change == "sign and significance")$direction
    infl_df <- influence_dfs$sig[[direction]]
    rerun_result <- RerunTargetModelForAlpha(
        infl_df, model_fit, "num_removed", sign_sig_alpha_target)
    rerun_sig <- sign(rerun_result$beta_pzse) == sign(rerun_result$beta_mzse)
    testthat::expect_true(sign(rerun_result$beta) != sign(base_vals["beta"]))
    testthat::expect_true(rerun_sig != base_sig)

    # For now just check that this runs
    RerunForTargetChanges(influence_dfs, target_change, model_fit)
  }

  num_obs <- 1000

  # Simulate with a small x scale to make for high influence.
  x <- runif(num_obs) * 0.01
  x <- x - mean(x)

  df <- GenerateRegressionData(num_obs, 0.5, x=matrix(x, ncol=1))
  df$arm <- as.integer(runif(num_obs) < 0.5)  # Control or treatment
  df$group <- sample.int(3, num_obs, replace=TRUE)  # Some subgrouping
  lm_result <- lm(y ~ x1 + 1, df, x=TRUE, y=TRUE)
  test_rerun(lm_result, paired=FALSE)
  test_rerun(lm_result, paired=TRUE)

  df$z1 <- df$x1 + 0.01 * rnorm(num_obs)
  df$arm <- as.integer(runif(num_obs) < 0.5)  # Control or treatment
  df$group <- sample.int(3, num_obs, replace=TRUE)  # Some subgrouping
  iv_result <- ivreg(y ~ x1 + 1 | z1 + 1, data=df, x=TRUE, y=TRUE)
  test_rerun(iv_result, paired=FALSE)
  test_rerun(lm_result, paired=TRUE)
})
