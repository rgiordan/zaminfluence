#base_dir <- Sys.getenv("REPO")
base_dir  <- "/home/rgiordan/Documents/git_repos/zaminfluence"
setwd(base_dir)

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


######################
# New code

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

