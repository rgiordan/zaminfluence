library(torch)
library(tidyverse)
library(zaminfluence)
library(AER)

library(devtools)
devtools::load_all("/home/rgiordan/Documents/git_repos/zaminfluence/zaminfluence")
zam_dir <- "/home/rgiordan/Documents/git_repos/zaminfluence"
#source(file.path(zam_dir, "zaminfluence/R/ols_iv_grads_lib.R"))
set.seed(42)

AssertNearlyZero <- function(v, tol=1e-8) {
    if (max(abs(v)) > tol) {
        stop(sprintf("Not approximately zero: %e > %e", max(abs(v)), tol))
    }
}



# Hmm
# https://discuss.pytorch.org/t/groupby-aggregate-mean-in-pytorch/45335

inds <- torch_tensor(c(1, 1, 3, 2, 2) %>% matrix(5, 1), dtype=torch_int64())
src <- torch_tensor(c(0.1) %>% matrix(5, 2))
agg_mat <- torch_tensor(matrix(0, 3, 2))

agg_mat$scatter_add(dim=1, index=inds[["repeat"]](list(1, 2)), src=src)





######################

num_obs <- 20
x_dim <- 3
param_true <- 0.1 * runif(x_dim)


# Generate data.
set.seed(42)
x_dim <- 3
param_true <- 0.1 * runif(x_dim)
df <- GenerateIVRegressionData(num_obs, param_true, num_groups=5)

# Fit an IV model.
x_names <- sprintf("x%d", 1:x_dim)
z_names <- sprintf("z%d", 1:x_dim)
iv_form <- formula(sprintf("y ~ %s - 1 | %s - 1",
                           paste(x_names, collapse=" + "),
                           paste(z_names, collapse=" + ")))
iv_res <- ivreg(data=df, formula = iv_form, x=TRUE, y=TRUE)



DefineTorchVars <- function(iv_vars) {
    stopifnot(all(c("x", "y", "z", "w0") %in% names(iv_vars)))
    
    num_obs <- nrow(iv_vars$x)
    num_cols <- ncol(iv_vars$x)
    
    stopifnot(length(iv_vars$y) == num_obs)
    stopifnot(nrow(iv_vars$x) == num_obs)
    stopifnot(nrow(iv_vars$z) == num_obs)
    stopifnot(ncol(iv_vars$x) == num_cols)
    stopifnot(ncol(iv_vars$z) == num_cols)
    
    x <- torch_tensor(iv_vars$x, requires_grad=FALSE, dtype=torch_double())
    z <- torch_tensor(iv_vars$z, requires_grad=FALSE, dtype=torch_double())
    y <- torch_tensor(matrix(iv_vars$y, ncol=1), requires_grad=FALSE, dtype=torch_double())
    w <- torch_tensor(matrix(iv_vars$w0, ncol=1), requires_grad=TRUE, dtype=torch_double())

    z_w <- z * w
    z_w_t <- z_w$transpose(2, 1)
    zwx <- torch_matmul(z_w_t, x)
    #t(iv_vars$z) %*% iv_vars$z - zwz # Why 1e-6?
    betahat <- torch::linalg_solve(zwx, torch_matmul(z_w_t, y))
    eps <- y - torch_matmul(x, betahat)

    return(list(
        num_obs=num_obs,
        num_cols=num_cols,
        x=x,
        z=z,
        y=y,
        w=w,
        zwx=zwx,
        eps=eps,
        betahat=betahat
    ))    
}


iv_vars <- GetIVVariables(iv_res)
tv <- DefineTorchVars(iv_vars)

tv$score_mat <- tv$z * tv$eps * tv$w
tv$se_group <- torch_tensor(as.integer(factor(df$se_group)) - 1, dtype=torch_int64())
# Problem is bincount does not accept vectors of weights
tv$score_agg <- torch_bincount(tv$se_group, tv$score_mat)
GroupedSum(iv_vars$z * as.numeric(tv$eps) * iv_vars$w, as.integer(factor(df$se_group)) - 1)
tv$score_agg

s_mat <- s_mat - rep(colMeans(s_mat), each=nrow(s_mat))





#if (!is.null(se_group)) {
    # Grouped standard errors
    # Enforces that the grouping variable is sequential and zero-indexed.
    
    score_mat <- GroupedSum(z * eps * w0, se_group)
    
    # colMeans(s_mat) is zero at the weights used for regression, but include
    # it so we can test partial derivatives.
    s_mat <- s_mat - rep(colMeans(s_mat), each=nrow(s_mat))
    
    num_groups <- nrow(s_mat)
    
    # v is for variance, so I take the average, but we then need to multiply
    # again by num_groups to get the standard error variance.
    v_mat <- t(s_mat) %*% s_mat / num_groups
    
    # Covariance matrix
    zwx_inv_vmat <- solve(zwx_qr, v_mat)
    se_mat <- solve(zwx_qr, t(zwx_inv_vmat)) * num_groups
# } else {
#     
# }
 