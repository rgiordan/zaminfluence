#!/usr/bin/env Rscript

library(AER)
library(zaminfluence)
library(sandwich)
library(testthat)
library(tidyverse)

context("zaminfluence")
source("utils.R")

# This requires the venv to have been set up precisely as described in the
# README.md.
venv_bin <- "../../../venv/bin/python3"

test_that("python runs", {
    if (!file.exists(venv_bin)) {
        print(paste0("Virtual environment for testing is missing (",
                     venv_bin,
                     "). ",
                     "Please follow the installation directions in the ",
                     "README.md file."))
    }
    zaminfluence::InitializePython(venv_bin)
})



test_that("regression works", {
  # reg_infl should be the output of ComputeRegressionInfluence
  TestConfiguration <- function(model_fit, se_group) {
    reg_infl <- ComputeModelInfluence(model_fit, se_group)

    grad_df <- GetTargetRegressorGrads(reg_infl, "x1")
    influence_dfs <- SortAndAccumulate(grad_df)

    # Check that the regressions match.
    reg <- tidy(reg_infl$model_fit)
    testthat::expect_equivalent(reg$estimate, reg_infl$betahat)

    # Check that the target index is correct..
    target_index <- attr(grad_df, "target_index")
    base_vals <- attr(grad_df, "base_vals")
    testthat::expect_equivalent(reg$estimate[target_index], base_vals["beta"])

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

    infl_scale <- GetInfluenceScale(grad_df$beta_grad)
    testthat::expect_equivalent(sqrt(n_obs) * robust_se, infl_scale)
  }

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
    iv_res <- ivreg(y ~ x1 + 1 | z1 + 1, data=df, x=TRUE, y=TRUE, weights=weights)
    TestConfiguration(iv_res, se_group=df[["se_group"]])
  }

  TestIVRegressionConfiguration(num_groups=NULL, weights=NULL)
  #TestIVRegressionConfiguration(num_groups=10, weights=NULL)
  TestIVRegressionConfiguration(num_groups=NULL, weights=runif(100))
  #TestIVRegressionConfiguration(num_groups=10, weights=runif(100))

})


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



test_that("pairing works", {
  set.seed(42)

  n_obs <- 1000
  df <- GenerateRegressionData(n_obs, 0.5)
  df$assignment <- as.integer(runif(n_obs) < 0.3)
  df$group1 <- c("a", "b")[as.integer(runif(n_obs) < 0.5) + 1]
  df$group2 <- c("c", "d")[as.integer(runif(n_obs) < 0.5) + 1]
  df$row <- 1:nrow(df)
  lm_result <- lm(y ~ x1 + 1, df, x=TRUE, y=TRUE)

  grad_df <-
      ComputeRegressionInfluence(lm_result) %>%
      GetTargetRegressorGrads("x1")
  paired_grad_df <-
    inner_join(grad_df, df, by="row") %>%
    CopyGradAttributes(grad_df)

  influence_dfs <-
      PairInfluenceScores(paired_grad_df,
                          assignment_col="assignment",
                          group_cols=c("group1", "group2"),
                          level0=0,
                          level1=1)
  for (change in c("sign", "sig")) {
      for (direction in c("pos", "neg")) {
          infl_df <- influence_dfs[[change]][[direction]]
          infl_df <- filter(infl_df, num_removed > 0)
          df_0 <- df[infl_df$row_0, ]
          df_1 <- df[infl_df$row_1, ]
          testthat::expect_true(all(df_0$assignment == 0))
          testthat::expect_true(all(df_1$assignment  == 1))
          testthat::expect_true(all(df_0$group1 == df_1$group1))
          testthat::expect_true(all(df_0$group2 == df_1$group2))
          testthat::expect_true(all(infl_df$obs_per_row == 2))
      }
  }
})


test_that("aggregated pairing works", {
  set.seed(42)

  n_obs <- 1000
  df <- GenerateRegressionData(n_obs, 0.5)
  df$arm <- as.integer(runif(n_obs) < 0.3)
  df$group <- sample.int(floor(n_obs / 10), n_obs, replace=TRUE)
  df$row <- 1:nrow(df)
  lm_result <- lm(y ~ x1 + 1, df, x=TRUE, y=TRUE)

  base_grad_df <-
      ComputeRegressionInfluence(lm_result) %>%
      GetTargetRegressorGrads("x1")

  grad_df <-
      base_grad_df %>%
      bind_cols(df[c("arm", "group")]) %>%
      group_by(arm, group) %>%
      summarize(se_grad=sum(se_grad),
                beta_grad=sum(beta_grad),
                beta_pzse_grad=sum(beta_pzse_grad),
                beta_mzse_grad=sum(beta_mzse_grad),
                obs_per_row=sum(obs_per_row),
                .groups="drop") %>%
      mutate(grouped_row=1:n()) %>%
      CopyGradAttributes(base_grad_df)

  # There is no longer a column that matches rows to the original data.
  attr(grad_df, "data_row_cols") <- "grouped_row"

  # However, we still need to record how many distinct gradients are being
  # sorted so we can index into grad_df using our results.
  attr(grad_df, "n_grad_rows") <- nrow(grad_df)

  influence_dfs <-
      PairInfluenceScores(grad_df,
                          assignment_col="arm",
                          group_cols="group",
                          level0=0,
                          level1=1)

  for (change in c("sign", "sig")) {
      for (direction in c("pos", "neg")) {
          infl_df <- influence_dfs[[change]][[direction]]
          infl_df <- filter(infl_df, num_removed > 0)
          df_0 <- df[infl_df$row_0, ]
          df_1 <- df[infl_df$row_1, ]
          testthat::expect_true(all(df_0$assignment == 0))
          testthat::expect_true(all(df_1$assignment  == 1))
          testthat::expect_true(all(df_0$group == df_1$group))
          testthat::expect_error(
            RerunForTargetChanges(influence_dfs, target_change, reg_fit),
            ".*")
      }
  }
})


test_that("regression moment conditions works", {
  set.seed(42)

  df <- GenerateRegressionData(100, 0.5)
  lm_result <- lm(y ~ x1 + 1, df, x=TRUE, y=TRUE)
  reg_moment_sens <- ComputeRegressionMomentSensitivity(lm_result)
  new_offset <- runif(ncol(lm_result$x))

  new_reg <- RegressWithOffset(lm_result, new_offset)
  actual_change <- new_reg$betahat - lm_result$coefficients
  pred_change <- as.numeric(reg_moment_sens$beta_grad %*% new_offset)

  # The dependence is actually linear, so the prediction is exact.
  testthat::expect_equivalent(actual_change, pred_change)
  stopifnot(max(abs(actual_change - pred_change)) < 1e-8)
})
