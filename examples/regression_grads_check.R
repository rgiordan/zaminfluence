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


#######################
# Let's do some checks

# Generate data.
x_dim <- 3
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


iv_res <- iv_fit



##################
# New code

AssertNearlyZero <- function(x, tol=1e-15) {
    stopifnot(max(abs(x)) < tol)
}


x <- iv_res$x$regressors
num_obs <- nrow(x)
z <- iv_res$x$instruments
y <- as.numeric(iv_res$y)
if (is.null(iv_res$weights)) {
    w0 <- rep(1.0, num_obs)
} else {
    w0 <- iv_res$weights
}
AssertNearlyZero(w0 - df$weights)





GetIVSEMat <- function(x, z, y, beta, w0) {
    z_w <- z * w0
    zwz <- t(z_w) %*% z
    zwx <- t(z_w) %*% x
    
    beta <- as.numeric(iv_res$coefficients)
    
    eps <- as.numeric(y - x %*% beta)
    sig2_hat <- sum(w0 * eps^2) / (num_obs - length(beta))
    
    zwx_inv_zwz <- solve(zwx, zwz)
    sand_mat <- solve(zwx, t(zwx_inv_zwz))
    se_mat <- sig2_hat * sand_mat 
    
    # Derivatives
    dsig2_hat_dw <- eps^2 / (num_obs - length(beta))
    
    return(list(se_mat=se_mat, sig2_hat=sig2_hat,
                dsig2_hat_dw=dsig2_hat_dw,
                sand_mat=sand_mat))
}


iv_se_list <- GetIVSEMat(x=x, z=z, y=y, beta=beta, w0=w0)

AssertNearlyZero(iv_se_list$sig2_hat - iv_res$sigma^2)
AssertNearlyZero(iv_se_list$se_mat - vcov(iv_res), tol=1e-14)

########################
# Test the derivatives

library(numDeriv)

dsig2_hat_dw_num <-
    numDeriv::jacobian(function(w) { GetIVSEMat(x=x, z=z, y=y, beta=beta, w0=w)$sig2_hat }, w0) %>%
    as.numeric()

AssertNearlyZero(dsig2_hat_dw_num - iv_se_list$dsig2_hat_dw, tol=1e-8)

#### copied...

z_w <- z * w0
zwz <- t(z_w) %*% z
zwx <- t(z_w) %*% x

beta <- as.numeric(iv_res$coefficients)

eps <- as.numeric(y - x %*% beta)
sig2_hat <- sum(w0 * eps^2) / (num_obs - length(beta))

zwx_inv_zwz <- solve(zwx, zwz)
sand_mat <- solve(zwx, t(zwx_inv_zwz))
se_mat <- sig2_hat * sand_mat 

#### copied ^


# See notes for the definition of these terms

AR_x <- zwx_inv_zwz %*% solve(t(zwx), t(x))
R_z <- solve(zwx, t(z))

dsand_mat_dw_num <-
    numDeriv::jacobian(function(w) {
        GetIVSEMat(x=x, z=z, y=y, beta=beta, w0=w)$sand_mat %>% diag() }, w0)

dsand_mat_dw <- R_z * R_z - 2 * AR_x * R_z

dsand_mat_dw - dsand_mat_dw_num


######################
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

