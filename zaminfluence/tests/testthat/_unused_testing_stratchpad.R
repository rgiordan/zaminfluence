#!/usr/bin/env Rscript

library(AER)
library(zaminfluence)
library(sandwich)
library(testthat)
library(tidyverse)



context("zaminfluence")


GetFitCovariance <- function(fit, se_group=NULL) {
  # Get a version of the sandwich covariance that should match
  # our computations.
  if (is.null(se_group)) {
    return(vcov(fit))
  } else {
    return(vcovCL(fit, cluster=se_group, type="HC0", cadjust=FALSE))
  }
}


###########################################################
# Some extra debugging, involving the fact that
# there is discontinuous behavior with the weights in vcovCL

RegressWithWeights <- function(w, se_group=NULL) {
  df$w <- w
  reg_fit <- lm(y ~ x1 + 1, data=df, x=TRUE, y=TRUE, weights=df$w)
  reg_vcov <- GetFitCovariance(reg_fit, se_group=se_group)
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

reg_vcov <- GetFitCovariance(reg_fit, se_group=se_group)
AssertNearlyEqual(reg_fit$coefficients, zam_reg_fit$betahat)
AssertNearlyEqual(reg_vcov, zam_reg_fit$se_mat)

max(abs(reg_vcov0 - reg_vcov)) # Check something actually changed
reg_vcov0 <- GetFitCovariance(reg_fit0, se_group=se_group)
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
GetFitCovariance(reg0)[1, 1] - zam_reg_fit0$se_mat[1, 1]
GetFitCovariance(reg1)[1, 1] - zam_reg_fit1$se_mat[1, 1]
GetFitCovariance(reg2)[1, 1] - zam_reg_fit2$se_mat[1, 1]

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
          testthat::expect_true(all(infl_df$arm_0 == 0))
          testthat::expect_true(all(infl_df$arm_1 == 1))

          # Check that the rows in the paired dataframe match the original.
          influence_cols <- c("se_grad", "beta_grad", "beta_pzse_grad", "beta_mzse_grad")
          AssertNearlyZero(
            grad_df[infl_df$grouped_row_0, influence_cols] -
              infl_df[, paste(influence_cols, "0", sep="_")]
          )
          AssertNearlyZero(
            grad_df[infl_df$grouped_row_1, influence_cols] -
              infl_df[, paste(influence_cols, "1", sep="_")]
          )

          # Check that the pairs worked as expected.
          testthat::expect_true(all(grad_df[infl_df$grouped_row_0, "arm"] == 0))
          testthat::expect_true(all(grad_df[infl_df$grouped_row_1, "arm"] == 1))
          testthat::expect_true(all(grad_df[infl_df$grouped_row_0, "group"] ==
                                    grad_df[infl_df$grouped_row_1, "group"]))
      }
  }
})




# Now that we're not supporting the Python stuff by default this
# we would need a separate test that only runs when Python is set up.
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
