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


AssertNearlyZero <- function(x, tol=1e-9) {
    x_norm <- max(abs(x))
    info_str <- sprintf("%e > %e", x_norm, tol)
    expect_true(x_norm < tol, info=info_str)
}


###############

ComputeRegressionInfluencePython <- function(lm_result, se_group=NULL) {
    py_main <- SetPythonRegressionVariables(lm_result, se_group=se_group)
    reg <- broom::tidy(lm_result)
    reticulate::py_run_string("
betahat = regsens_rgiordandev.reg(y, x, w=w0)
se, betahat_grad, se_grad = regsens_rgiordandev.get_regression_w_grads(
    betahat, y, x, w0, se_group=se_group)
")
    if (max(abs(py_main$betahat - reg$estimate)) > 1e-8) {
        warning("Regression coefficients do not match.")
    }
    
    # Note that the standard errors may not match lm_result when using se_group.
    return(list(model_fit=lm_result,
                n_obs=nrow(lm_result$x),
                regressor_names=colnames(lm_result$x),
                grad_fun="get_regression_w_grads",
                
                betahat=py_main$betahat,
                se=py_main$se,
                weights=py_main$w0,
                
                beta_grad=py_main$betahat_grad,
                se_grad=py_main$se_grad)
    )
}


ComputeIVRegressionInfluencePython <- function(iv_res, se_group=NULL) {
    py_main <- SetPythonIVRegressionVariables(iv_res, se_group=se_group)
    reg <- broom::tidy(iv_res)
    reticulate::py_run_string("
betahat = iv_lib.iv_reg(y, x, z, w=w0)
se, betahat_grad, se_grad = iv_lib.get_iv_regression_w_grads(
    betahat, y, x, z, w0, se_group=se_group)
")
    if (max(abs(py_main$betahat - reg$estimate)) > 1e-8) {
        warning("Regression coefficients do not match.")
    }
    
    # Note that the standard errors may not match iv_res when using se_group.
    return(list(model_fit=iv_res,
                n_obs=length(iv_res$y),
                regressor_names=colnames(iv_res$x$regressors),
                grad_fun="get_iv_regression_w_grads",
                
                betahat=py_main$betahat,
                se=py_main$se,
                weights=py_main$w0,
                
                beta_grad=py_main$betahat_grad,
                se_grad=py_main$se_grad)
    )
}


###################

do_iv <- FALSE

x_dim <- 3
beta_true <- runif(x_dim)
num_groups <- NULL

if (do_iv) {
    df <- GenerateIVRegressionData(20, beta_true, num_groups=num_groups)
} else {
    df <- GenerateRegressionData(20, beta_true, num_groups=num_groups)
}
df$weights <- runif(nrow(df)) + 1
df$se_group <- 1:nrow(df) - 1
df$se_group <- NULL

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
}


#######################
# Let's do some checks

reg_infl <- ComputeModelInfluence(reg_fit)
grad_df <- GetTargetRegressorGrads(iv_infl, "x1")
influence_dfs <- SortAndAccumulate(grad_df)

if (do_iv) {
    reg_infl_python <- ComputeIVRegressionInfluencePython(reg_fit)
} else {
    reg_infl_python <- ComputeRegressionInfluencePython(reg_fit)
}

AssertNearlyZero(reg_infl_python$betahat - reg_infl$betahat)
AssertNearlyZero(reg_infl_python$beta_grad - reg_infl$beta_grad)
AssertNearlyZero(reg_infl_python$se_grad - reg_infl$se_grad)


