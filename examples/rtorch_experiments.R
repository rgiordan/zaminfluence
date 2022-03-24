library(torch)
library(tidyverse)
library(zaminfluence)
set.seed(42)

AssertNearlyZero <- function(v, tol=1e-8) {
    if (max(abs(v)) > tol) {
        stop(sprintf("Not approximately zero: %e > %e", max(abs(v)), tol))
    }
}


num_obs <- 1000
beta_dim <- 10
param_true <- 0.1 * runif(beta_dim)

use_fe <- TRUE

if (use_fe) {
    df <- GenerateRegressionData(num_obs, param_true, num_groups=100) %>%
        rename(z=se_group) %>%
        mutate(z=factor(z))
    
    # Fit a regression model.
    x_names <- sprintf("x%d", 1:beta_dim)
    reg_form <- formula(sprintf("y ~ %s + z - 1", paste(x_names, collapse=" + ")))
    fit_object <- lm(data = df, formula=reg_form, x=TRUE, y=TRUE)
    
} else {
    x_names <- sprintf("x%d", 1:beta_dim)
    reg_form <- formula(sprintf("y ~ %s - 1", paste(x_names, collapse=" + ")))
    df <- GenerateRegressionData(num_obs, param_true, num_groups=NULL)
}

fit_object <- lm(data=df, formula=reg_form, x=TRUE, y=TRUE)

model_grads <- 
    zaminfluence::ComputeModelInfluence(fit_object) %>%
    zaminfluence::AppendTargetRegressorInfluence("x1")


zam_time <- Sys.time()
grouped_model_grads <- 
    zaminfluence::ComputeModelInfluence(fit_object, se_group=1:num_obs) %>%
    zaminfluence::AppendTargetRegressorInfluence("x1")
zam_time <- Sys.time() - zam_time



##############
# Torch

DefineTorchVars <- function(x, y, w) {
    num_obs <- nrow(x)
    num_cols <- ncol(x)
    x_tens <- torch_tensor(x, requires_grad=FALSE, dtype=torch_double())
    y_tens <- torch_tensor(matrix(y, ncol=1), requires_grad=FALSE, dtype=torch_double())
    w_tens <- torch_tensor(matrix(w, ncol=1), requires_grad=TRUE, dtype=torch_double())
    
    xw <- x_tens * w_tens
    xwx <- torch_matmul(torch_transpose(xw, 1, 2), x_tens)
    xwy <- torch_matmul(torch_transpose(xw, 1, 2), y_tens)
    betahat <- linalg_solve(xwx, xwy)
    eps <- y_tens - torch_matmul(x_tens, betahat)
    return(list(
        n_obs=num_obs,
        x_dim=num_cols,
        x=x_tens,
        y=y_tens,
        w=w_tens,
        xwx=xwx,
        eps=eps,
        betahat=betahat
    ))    
}

w_base <- rep(1, num_obs)


##############
# Grads

torch_time <- Sys.time()
tv <- DefineTorchVars(x=fit_object$x, y=fit_object$y, w=w_base)

keep_dims <- 1:min(10, tv$x_dim)
#keep_dims <- 1:tv$x_dim

# betahat
torch_infl_mat <- matrix(NA, nrow=length(keep_dims), ncol=num_obs)
for (di in 1:length(keep_dims)) {
    d <- keep_dims[di]
    torch_infl_mat[di, ] <- 
        torch::autograd_grad(
            tv$betahat[d], tv$w, retain_graph=TRUE)[[1]] %>% as.numeric()
}

# grouped ses
tv$score_mat  <- tv$eps * tv$x
tv$score_cov <- torch_matmul(
    torch_transpose(tv$w * tv$score_mat, 2, 1), tv$w * tv$score_mat) / num_obs
tv$xwx_inv_score_cov <- linalg_solve(tv$xwx / num_obs, tv$score_cov)
tv$betahat_grouped_cov <- 
    linalg_solve(tv$xwx / num_obs, torch_transpose(tv$xwx_inv_score_cov, 1, 2)) / num_obs
tv$betahat_grouped_se <- torch_sqrt(torch_diag(tv$betahat_grouped_cov))

torch_grouped_se_infl_mat <- matrix(NA, nrow=length(keep_dims), ncol=num_obs)
for (di in 1:length(keep_dims)) {
    d <- keep_dims[di]
    torch_grouped_se_infl_mat[di, ] <- 
        torch::autograd_grad(
            tv$betahat_grouped_se[d],
            tv$w, retain_graph=TRUE)[[1]] %>% as.numeric()
}

torch_time <- Sys.time() - torch_time

# ses (don't time this one)
tv$sig2hat <- torch_sum(tv$w * (tv$eps ** 2))  / (num_obs - tv$x_dim)
tv$betahat_cov <- torch_inverse(tv$xwx) * tv$sig2hat
tv$betahat_se <- torch_sqrt(torch_diag(tv$betahat_cov))

torch_se_infl_mat <- matrix(NA, nrow=length(keep_dims), ncol=num_obs)
for (di in 1:length(keep_dims)) {
    d <- keep_dims[di]
    torch_se_infl_mat[di, ] <- 
        torch::autograd_grad(
            tv$betahat_se[d], tv$w, retain_graph=TRUE)[[1]] %>% as.numeric()
}

#########################
# Compare times

AssertNearlyZero(torch_infl_mat - model_grads$param_grad[keep_dims, ])
AssertNearlyZero(torch_infl_mat - grouped_model_grads$param_grad[keep_dims, ])
AssertNearlyZero(torch_grouped_se_infl_mat - grouped_model_grads$se_grad[keep_dims, ])
AssertNearlyZero(torch_se_infl_mat - model_grads$se_grad[keep_dims, ])


torch_time
zam_time



##################################################
# Can matrix inverses be cached?  I don't see how.

# TODO: try
# https://torch.mlverse.org/docs/reference/torch_cholesky_solve.html

mat_dim <- 500
a_mat <- runif(mat_dim * mat_dim) %>% matrix(mat_dim, mat_dim)
a_mat <- a_mat + t(a_mat) + diag(a_mat) 
b_vec <- runif(mat_dim)
ainv_b <- solve(a_mat, b_vec)

# Base solve

a_tens <- torch_tensor(a_mat, dtype=torch_double())
b_tens <- torch_tensor(b_vec %>% matrix(ncol=1), dtype=torch_double())

solve_time <- Sys.time()
ainv_b_torch <- torch::linalg_solve(a_tens, b_tens) %>% as.numeric()
solve_time <- Sys.time() - solve_time
AssertNearlyZero(ainv_b_torch -  ainv_b)

# QR
a_tens <- torch_tensor(a_mat, dtype=torch_double())
b_tens <- torch_tensor(b_vec %>% matrix(ncol=1), dtype=torch_double())

qr_time <- Sys.time()
qr_tens <- linalg_qr(a_tens)
qr_solve_time <- Sys.time()
qtb_tens <- torch::torch_matmul(torch_transpose(qr_tens[[1]], 2, 1), b_tens)
ainv_b_qr <- 
    torch::linalg_solve(qr_tens[[2]], qtb_tens) %>% as.numeric()
qr_time <- Sys.time() - qr_time
qr_solve_time <- Sys.time() - qr_solve_time

AssertNearlyZero(ainv_b_qr -  ainv_b)
AssertNearlyZero(ainv_b_qr -  ainv_b_torch)
print(solve_time)
print(qr_solve_time)
print(qr_time)



##################################################
# Be careful about cyclic dependencies.


cat("-----------\n")
x <- torch_tensor(5, requires_grad=TRUE)
print(sprintf("x = %f", x))
x <- 2 * x
print(sprintf("x = %f", x))
y <- 2 * x 
print(sprintf("x = %f", x))
print(sprintf("y = %f", y))
print(sprintf("dy/dx = %f", torch::autograd_grad(y, x)[[1]]))



##################################################
# Can you force it to compute a Hessian?  Yes.


par_dim <- 5

x_r <- runif(par_dim) %>% matrix(ncol=1)
a_r <- runif(par_dim * par_dim) %>% matrix(par_dim, par_dim)
a_r <- a_r + t(a_r)

x <- torch_tensor(x_r, dtype=torch_float64(), requires_grad=TRUE)
a <- torch_tensor(a_r, dtype=torch_float64(), requires_grad=FALSE)

y <- 0.5 * torch_einsum("ik,ij,jk", list(x, a, x))
AssertNearlyZero(as.numeric(y) - 0.5 * t(x_r) %*% a_r %*% x_r, tol=1e-7) # why?

y_grad <- torch::autograd_grad(y, x, retain_graph=TRUE, create_graph=TRUE)

y_hess <- matrix(NA, nrow=par_dim, ncol=par_dim)
for (d in 1:par_dim) {
    y_hess[d, ] <- 
        torch::autograd_grad(y_grad[[1]][d], x, retain_graph=TRUE)[[1]] %>% as.numeric()
}

AssertNearlyZero(y_hess - a_r)
