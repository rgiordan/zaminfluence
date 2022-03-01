library(torch)
library(tidyverse)
library(zaminfluence)

num_obs <- 10000
set.seed(42)
x_dim <- 10
param_true <- 0.1 * runif(x_dim)
df <- GenerateRegressionData(num_obs, param_true, num_groups=NULL)

x_names <- sprintf("x%d", 1:x_dim)
reg_form <- formula(sprintf("y ~ %s - 1", paste(x_names, collapse=" + ")))
fit_object <- lm(data=df, formula=reg_form, x=TRUE, y=TRUE)

model_grads <- 
    zaminfluence::ComputeModelInfluence(fit_object) %>%
    zaminfluence::AppendTargetRegressorInfluence("x1")


zam_time <- Sys.time()
grouped_model_grads <- 
    zaminfluence::ComputeModelInfluence(fit_object, se_group=1:num_obs) %>%
    zaminfluence::AppendTargetRegressorInfluence("x1")
zam_time <- Sys.time() - zam_time

AssertNearlyZero <- function(v, tol=1e-8) {
    if (max(abs(v)) > tol) {
        stop(sprintf("Not approximately zero: %f > %f", max(abs(v)), tol))
    }
}




##############
# Torch

DefineTorchVars <- function(x, y, w) {
    num_obs <- nrow(x)
    x_tens <- torch_tensor(x, requires_grad=FALSE)
    y_tens <- torch_tensor(matrix(y, ncol=1), requires_grad=FALSE)
    w_tens <- torch_tensor(matrix(w, ncol=1), requires_grad=TRUE)
    
    xw <- x_tens * w_tens
    xwx <- torch_matmul(torch_transpose(xw, 1, 2), x_tens)
    xwy <- torch_matmul(torch_transpose(xw, 1, 2), y_tens)
    betahat <- linalg_solve(xwx, xwy)
    eps <- y_tens - torch_matmul(x_tens, betahat)
    return(list(
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

# betahat
tv <- DefineTorchVars(x=fit_object$x, y=fit_object$y, w=w_base)
torch_infl_mat <- matrix(NA, nrow=x_dim, ncol=num_obs)
for (d in 1:x_dim) {
    torch_infl_mat[d, ] <- 
        torch::autograd_grad(
            tv$betahat[d], tv$w, retain_graph=TRUE)[[1]] %>% as.numeric()
}
AssertNearlyZero(torch_infl_mat - model_grads$param_grad)
AssertNearlyZero(torch_infl_mat - grouped_model_grads$param_grad)


# ses
tv$sig2hat <- torch_sum(tv$w * (tv$eps ** 2))  / (num_obs - x_dim)
tv$betahat_cov <- torch_inverse(tv$xwx) * tv$sig2hat
tv$betahat_se <- torch_sqrt(torch_diag(tv$betahat_cov))

torch_se_infl_mat <- matrix(NA, nrow=x_dim, ncol=num_obs)
for (d in 1:x_dim) {
    torch_se_infl_mat[d, ] <- 
        torch::autograd_grad(
            tv$betahat_se[d], tv$w, retain_graph=TRUE)[[1]] %>% as.numeric()
}
AssertNearlyZero(torch_se_infl_mat - model_grads$se_grad)

# grouped ses
tv$score_mat  <- tv$eps * tv$x
tv$score_cov <- torch_matmul(
    torch_transpose(tv$w * tv$score_mat, 2, 1), tv$w * tv$score_mat) / num_obs
tv$xwx_inv_score_cov <- linalg_solve(tv$xwx / num_obs, tv$score_cov)
tv$betahat_grouped_cov <- 
    linalg_solve(tv$xwx / num_obs, torch_transpose(tv$xwx_inv_score_cov, 1, 2)) / num_obs
tv$betahat_grouped_se <- torch_sqrt(torch_diag(tv$betahat_grouped_cov))

torch_grouped_se_infl_mat <- matrix(NA, nrow=x_dim, ncol=num_obs)
for (d in 1:x_dim) {
    torch_grouped_se_infl_mat[d, ] <- 
        torch::autograd_grad(
            tv$betahat_grouped_se[d], tv$w, retain_graph=TRUE)[[1]] %>% as.numeric()
}
AssertNearlyZero(torch_grouped_se_infl_mat - grouped_model_grads$se_grad)


torch_time <- Sys.time() - torch_time


torch_time
zam_time
