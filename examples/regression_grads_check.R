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

py_main <- InitializePython(file.path(base_dir, "venv/bin/python3"))
compare <- function(x, y) { return(max(abs(x - y))) }
check_equivalent  <- function(x, y) { stopifnot(compare(x, y) < 1e-8) }

n_obs <- 10000

# The test utilities can simulate data.
source(file.path(base_dir, "zaminfluence/tests/testthat/utils.R"))

set.seed(42)

x_dim <- 3
beta_true <- runif(x_dim)
df <- GenerateRegressionData(n_obs, beta_true, num_groups=NULL)

# Fit a regression model.
x_names <- sprintf("x%d", 1:x_dim)
reg_form <- formula(sprintf("y ~ %s - 1", paste(x_names, collapse=" + ")))
reg_fit <- lm(data = df, formula = reg_form, x=TRUE, y=TRUE)


#######################
# Let's do some checks
# 
# # Generate data.
# x_dim <- 3
# beta_true <- runif(x_dim)
# df <- GenerateIVRegressionData(10, beta_true, num_groups=NULL)
# 
# # Fit an IV model.
# x_names <- sprintf("x%d", 1:x_dim)
# z_names <- sprintf("z%d", 1:x_dim)
# iv_form <- formula(sprintf("y ~ %s - 1 | %s - 1",
#                            paste(x_names, collapse=" + "),
#                            paste(z_names, collapse=" + ")))
# df$weights <- runif(nrow(df)) + 1
# iv_fit <- ivreg(data = df, formula = iv_form, x=TRUE, y=TRUE, weights=weights)
# 
# # Get influence.
# iv_infl <- ComputeModelInfluence(iv_fit)
# grad_df <- GetTargetRegressorGrads(iv_infl, "x1")
# influence_dfs <- SortAndAccumulate(grad_df)
# 


beta_true <- runif(x_dim)
df <- GenerateIVRegressionData(10, beta_true, num_groups=NULL)

# Fit an IV model.
x_names <- sprintf("x%d", 1:x_dim)
z_names <- sprintf("z%d", 1:x_dim)
iv_form <- formula(sprintf("y ~ %s - 1 | %s - 1",
                           paste(x_names, collapse=" + "),
                           paste(z_names, collapse=" + ")))
df$weights <- runif(nrow(df)) + 1
iv_fit <- ivreg(data = df, formula = iv_form, x=TRUE, y=TRUE, weights=weights)

# Get influence.
iv_infl <- ComputeModelInfluence(iv_fit)
grad_df <- GetTargetRegressorGrads(iv_infl, "x1")
influence_dfs <- SortAndAccumulate(grad_df)

iv_se_list <- GetIVSEDerivs(
    x=df[, x_names] %>% as.matrix(),
    z=df[, z_names] %>% as.matrix(),
    y=df$y,
    beta=iv_fit$coefficients,
    w0=df$weights, testing=TRUE)



##################
# New code

AssertNearlyZero <- function(x, tol=1e-15) {
    stopifnot(max(abs(x)) < tol)
}

# 
# x <- iv_fit$x$regressors
# num_obs <- nrow(x)
# z <- iv_fit$x$instruments
# y <- as.numeric(iv_fit$y)
# beta <- as.numeric(iv_fit$coefficients)
# if (is.null(iv_fit$weights)) {
#     w0 <- rep(1.0, num_obs)
# } else {
#     w0 <- iv_fit$weights
# }
# AssertNearlyZero(w0 - df$weights)


###########################
# Load and test

#source(file.path(base_dir, "zaminfluence/R/ols_iv_grads_lib.R"))

iv_se_list <- GetIVSEDerivs(
    x=df[, x_names] %>% as.matrix(),
    z=df[, z_names] %>% as.matrix(),
    y=df$y,
    beta=iv_fit$coefficients,
    w0=df$weights, testing=TRUE)


AssertNearlyZero(iv_se_list$betahat - iv_fit$coefficients, tol=1e-11)
AssertNearlyZero(iv_se_list$sig2_hat - iv_fit$sigma^2)
AssertNearlyZero(iv_se_list$se_mat - vcov(iv_fit), tol=1e-11)
AssertNearlyZero(iv_se_list$se - vcov(iv_fit) %>% diag() %>% sqrt(), tol=1e-11)

########################
# Test the derivatives

# Test the partial derivatives

library(numDeriv)

dsig2_hat_dw_num <-
    numDeriv::jacobian(function(w) {
        GetIVSEDerivs(x=x, z=z, y=y, beta=beta, w0=w, testing=TRUE)$sig2_hat }, w0) %>%
    as.numeric()
AssertNearlyZero(dsig2_hat_dw_num - iv_se_list$dsig2_hat_dw_partial, tol=1e-8)

dsand_mat_diag_dw_num <-
    numDeriv::jacobian(function(w) {
        GetIVSEDerivs(x=x, z=z, y=y, beta=beta, w0=w, testing=TRUE)$sand_mat %>% diag() }, w0)
AssertNearlyZero(dsand_mat_diag_dw_num - iv_se_list$dsand_mat_diag_dw_partial, tol=1e-6)

dse_mat_diag_dw_num <-
    numDeriv::jacobian(function(w) {
        GetIVSEDerivs(x=x, z=z, y=y, beta=beta, w0=w, testing=TRUE)$se_mat %>% diag() }, w0)
AssertNearlyZero(dse_mat_diag_dw_num - iv_se_list$dse_mat_diag_dw_partial, tol=5e-6)

dsig2_hat_dbeta_num <-
    numDeriv::jacobian(function(beta) {
        GetIVSEDerivs(x=x, z=z, y=y, beta=beta, w0=w0, testing=TRUE)$sig2_hat }, beta)
AssertNearlyZero(dsig2_hat_dbeta_num - iv_se_list$dsig2_hat_dbeta, tol=1e-8)


# Test the full derivatives

GetIVTestResults <- function(w) {
    df_test <- df
    df_test$weights <- w
    beta_test <- ivreg(data=df_test, formula=iv_form,
                       x=TRUE, y=TRUE, weights=weights)$coefficients
    iv_se_test_list <- GetIVSEDerivs(x=x, z=z, y=y, beta=beta_test, w0=w, testing=TRUE)
    return(
        list(beta=iv_se_test_list$betahat,
             se=iv_se_test_list$se,
             sig2_hat=iv_se_test_list$sig2_hat,
             se_mat_diag=diag(iv_se_test_list$se_mat))
    )
}


GetIVTestResults(df$weights)
GetIVTestResults(df$weights + runif(nrow(df)))

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
print(dse_mat_diag_dw_num - iv_se_list$dse_mat_diag_dw)
# Check the relative error for this one
AssertNearlyZero((dse_mat_diag_dw_num - iv_se_list$dse_mat_diag_dw) /
                     (iv_se_list$dse_mat_diag_dw + 1e-3), tol=1e-6)

dse_dw_num <- numDeriv::jacobian(function(w) { GetIVTestResults(w)$se }, w0)
#plot(dse_dw_num, iv_se_list$dse_dw); abline(0, 1)
# Check the relative error for this one
AssertNearlyZero((dse_dw_num - iv_se_list$dse_dw) / (iv_se_list$dse_dw + 1e-3), tol=1e-6)







######################
# Combine to get grads





#######################################
#######################################
#######################################


#######################################
# New code for regression

lm_res <- reg_fit

# This is duplicated in python_lib
GetLMWeights <- function(lm_res) {
    n_obs <- nrow(lm_res$x)
    if (is.null(lm_res$weights)) {
        weights <- rep(1.0, n_obs)
    } else {
        weights <- lm_res$weights
    }
    return(weights)
}


# Provide our own implementation so you can be sure
# it matches our gradient computations
GetStandardErrorMatrix <- function(betahat, y, x, w, se_group=NULL) {
    
}


# get_regression_w_grads
GetBetaGrads <- function(lm_res, se_group=NULL) {
    n_obs <- nrow(lm_res$x)
    w0 <- GetLMWeights(lm_res)
    x <- as.matrix(lm_res$x)
    y <- lm_res$y
    betahat <- as.numeric(lm_res$coefficients)

    x_w <- x * w0
    xtwx <- t(x) %*% x_w
    yhat <- as.numeric(x_w %*% betahat)
    eps <- y - yhat
    x_times_w_eps <- eps * x_w
    betahat_grad <- solve(xtwx, t(x_times_w_eps))
    
    # Compute standard errors
    if (is.null(se_group)) {
        # The derivative of sigma2hat wrt beta is zero.
        # But there is explicit w dependence in xtx_bar.
        xtx_bar <- xtwx / num_obs
        sigma2hat <- sum(w0 * (eps^2)) / (num_obs - length(betahat))
        
        
        xtwx_inv <- solve(xtx_bar)
        se2 <- sigma2hat * xtwx_inv / num_obs
    } else {
        stop("Not implemented yet")
    }
    
    se <- sqrt(se2)
    se_grad <- 
    
    return(list(betahat_grad=betahat_grad, se_grad=se_grad))
}

beta_grad <- GetBetaGrad(reg_fit)
infl_list <- ComputeRegressionInfluence(reg_fit)
max(abs(infl_list$beta_grad - beta_grad))
max(abs(infl_list$se_grad - se_grad))
max(abs(infl_list$se_grad))
max(abs(infl_list$beta_grad))

#plot(beta_grad, infl_list$beta_grad); abline(0, 1)

