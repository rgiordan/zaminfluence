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

se_group <- df$se_group
#se_group <- NULL

model_grads <- 
    zaminfluence::ComputeModelInfluence(iv_res, se_group=se_group) %>%
    zaminfluence::AppendTargetRegressorInfluence("x1")

keep_inds <- 1:x_dim


# Torch stuff

TorchGroupedAggregate <- function(src_mat, inds) {
    # https://discuss.pytorch.org/t/groupby-aggregate-mean-in-pytorch/45335
    max_ind <- inds$max() %>% as.integer()
    zero_mat <- torch_tensor(
        matrix(0, nrow=max_ind, ncol=src_mat$shape[2]),
        dtype=src_mat$dtype)
    inds_rep <- inds[["repeat"]](list(1, src_mat$shape[2]))
    agg_mat <- zero_mat$scatter_add(dim=1, index=inds_rep, src=src_mat)
    return(agg_mat)    
}

# inds <- torch_tensor(c(1, 1, 3, 2, 2) %>% matrix(5, 1), dtype=torch_int64())
# src <- torch_tensor(c(rep(0.1, 5), rep(0.2, 5)) %>% matrix(5, 2))
# TorchGroupedAggregate(src, inds)


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
        zw_t=z_w_t,
        eps=eps,
        betahat=betahat
    ))    
}

GetIVRegressionSEDerivs <- function(iv_vars, se_group=NULL, keep_inds=NULL) {
    keep_inds <- ValidateBetaInds(keep_inds, iv_vars$x)
    tv <- DefineTorchVars(iv_vars)
    
    if (!is.null(se_group)) {
        tv$score_mat <- tv$z * tv$eps * tv$w
        tv$se_group <-
            torch_tensor(
                as.integer(factor(df$se_group)) %>% matrix(ncol=1), 
                dtype=torch_int64())
        tv$score_sum <- TorchGroupedAggregate(src_mat=tv$score_mat, inds=tv$se_group)
        tv$s_mat <- tv$score_sum - tv$score_sum$mean(dim=1, keepdim=TRUE)
        
        num_groups <- tv$s_mat$shape[1]
        tv$v_mat <- torch_matmul(tv$s_mat$transpose(2, 1), tv$s_mat) / num_groups
        
        tv$zwx_inv_vmat <- linalg_solve(tv$zwx, tv$v_mat)
        tv$se_cov_mat <- linalg_solve(tv$zwx, torch_transpose(tv$zwx_inv_vmat, 1, 2)) * num_groups
    } else {
        tv$sig2_hat <- torch_sum(tv$w * (tv$eps ** 2)) / (tv$num_obs - tv$num_cols)
        
        if (iv_vars$z_equals_x) {
            # Regression is when Z == X and we can save some computation
            tv$se_cov_mat <- tv$sig2_hat * torch_inverse(tv$zwx)
        } else {
            # IV is when Z != X
            tv$zwz <- torch_matmul(tv$zw_t, tv$z)
            tv$zwx_inv_zwz <- linalg_solve(tv$zwx, tv$zwz)
            tv$se_cov_mat <- tv$sig2_hat * linalg_solve(tv$zwx, tv$zwx_inv_zwz$transpose(2, 1))
        }
    }
    
    
    # betahat
    betahat_infl_mat <- matrix(NA, nrow=length(keep_inds), ncol=num_obs)
    for (di in 1:length(keep_inds)) {
        d <- keep_inds[di]
        betahat_infl_mat[di, ] <- 
            torch::autograd_grad(
                tv$betahat[d], tv$w, retain_graph=TRUE)[[1]] %>% as.numeric()
    }
    
    # standard errors
    tv$betahat_se <- torch_sqrt(torch_diag(tv$se_cov_mat))
    se_infl_mat <- matrix(NA, nrow=length(keep_inds), ncol=num_obs)
    for (di in 1:length(keep_inds)) {
        d <- keep_inds[di]
        se_infl_mat[di, ] <- 
            torch::autograd_grad(
                tv$betahat_se[d],
                tv$w, retain_graph=TRUE)[[1]] %>% as.numeric()
    }
    
    return(list(
        tv=tv, 
        betahat=as.numeric(tv$betahat),
        betahat_se=as.numeric(tv$betahat_se),
        betahat_infl_mat=betahat_infl_mat, 
        se_infl_mat=se_infl_mat))
}



# Compute

z_equals_x <- FALSE
iv_vars <- GetIVVariables(iv_res)
iv_vars$z_equals_x <- FALSE

torch_grads <- GetIVRegressionSEDerivs(iv_vars, se_group=se_group, keep_inds=keep_inds)

AssertNearlyZero(model_grads$param_grad - torch_grads$betahat_infl_mat)
AssertNearlyZero(model_grads$se_grad - torch_grads$se_infl_mat)


#plot(model_grads$se_grad, se_infl_mat); abline(0, 1)

 