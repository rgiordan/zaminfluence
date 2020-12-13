#base_dir <- Sys.getenv("REPO")
base_dir  <- "/home/rgiordan/Documents/git_repos/zaminfluence"
setwd(base_dir)
source(file.path(base_dir, "zaminfluence/tests/testthat/utils.R"))

library(ggplot2)
library(dplyr)
library(reshape2)
library(gridExtra)
library(sandwich)
library(zaminfluence)
library(AER)
library(numDeriv)

library(devtools)
devtools::load_all("zaminfluence")

py_main <- InitializePython(file.path(base_dir, "venv/bin/python3"))
compare <- function(x, y) { return(max(abs(x - y))) }
check_equivalent  <- function(x, y) { stopifnot(compare(x, y) < 1e-8) }

n_obs <- 10000

# The test utilities can simulate data.
source(file.path(base_dir, "zaminfluence/tests/testthat/utils.R"))

set.seed(42)


#######################
# Let's do some checks


x_dim <- 3
beta_true <- runif(x_dim)
#df <- GenerateIVRegressionData(10, beta_true, num_groups=NULL)
df <- GenerateRegressionData(10, beta_true, num_groups=NULL)
df$weights <- runif(nrow(df)) + 1

# Fit a model.
# x_names <- sprintf("x%d", 1:x_dim)
# z_names <- sprintf("z%d", 1:x_dim)
# iv_form <- formula(sprintf("y ~ %s - 1 | %s - 1",
#                            paste(x_names, collapse=" + "),
#                            paste(z_names, collapse=" + ")))
# reg_fit <- ivreg(data = df, formula = iv_form, x=TRUE, y=TRUE, weights=weights)
x_names <- sprintf("x%d", 1:x_dim)
reg_form <- formula(sprintf("y ~ %s - 1",
                           paste(x_names, collapse=" + ")))
reg_fit <- lm(data=df, formula=reg_form, x=TRUE, y=TRUE, weights=weights)

# Get influence.
# iv_infl <- ComputeModelInfluence(reg_fit)
# grad_df <- GetTargetRegressorGrads(iv_infl, "x1")
# influence_dfs <- SortAndAccumulate(grad_df)

reg_se_list <- GetRegressionSEDerivs(
    x=df[, x_names] %>% as.matrix(),
    y=df$y,
    beta=reg_fit$coefficients,
    w0=df$weights, testing=TRUE)



##################
# New code

AssertNearlyZero <- function(x, tol=1e-15) {
    stopifnot(max(abs(x)) < tol)
}


###########################
# Load and test

AssertNearlyZero(reg_se_list$betahat - reg_fit$coefficients, tol=1e-11)
AssertNearlyZero(reg_se_list$sig2_hat - sigma(reg_fit)^2)
AssertNearlyZero(reg_se_list$se_mat - vcov(reg_fit), tol=1e-11)
AssertNearlyZero(reg_se_list$se - vcov(reg_fit) %>% diag() %>% sqrt(), tol=1e-11)



LocalGetRegressionSEDerivs <- function(w=df$weights, beta=reg_fit$coefficients) {
    GetRegressionSEDerivs(
        x=df[, x_names] %>% as.matrix(),
        y=df$y,
        beta=beta,
        w0=w,
        testing=TRUE)
}


# Test the partial derivatives
w0 <- df$weights
beta <- reg_fit$coefficients
dsig2_hat_dw_num <-
    numDeriv::jacobian(function(w) {
        LocalGetRegressionSEDerivs(w=w)$sig2_hat
    }, w0) %>%
    as.numeric()
AssertNearlyZero(dsig2_hat_dw_num - iv_se_list$dsig2_hat_dw_partial, tol=1e-8)

dsand_mat_diag_dw_num <-
    numDeriv::jacobian(function(w) {
        LocalGetRegressionSEDerivs(w=w)$sand_mat %>% diag()
    }, w0)
AssertNearlyZero(dsand_mat_diag_dw_num - iv_se_list$dsand_mat_diag_dw_partial, tol=1e-6)

dse_mat_diag_dw_num <-
    numDeriv::jacobian(function(w) {
        LocalGetRegressionSEDerivs(w=w)$se_mat %>% diag()
    }, w0)
AssertNearlyZero(dse_mat_diag_dw_num - iv_se_list$dse_mat_diag_dw_partial, tol=5e-6)

dsig2_hat_dbeta_num <-
    numDeriv::jacobian(function(beta) {
        LocalGetRegressionSEDerivs(beta=beta)$sig2_hat
    }, beta)
AssertNearlyZero(dsig2_hat_dbeta_num - iv_se_list$dsig2_hat_dbeta, tol=1e-8)

#############################
# Test the full derivatives
GetIVTestResults <- function(w) {
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

