library(torch)
library(tidyverse)
library(zaminfluence)

if (FALSE) {
    x <- torch_tensor(1, requires_grad = TRUE)
    w <- torch_tensor(2, requires_grad = TRUE)
    b <- torch_tensor(3, requires_grad = TRUE)
    y <- w * x + b
    y$backward()
    x$grad
    w$grad
    b$grad
    
    # Fine
    
    x <- torch_tensor(c(1, 2, 3))
    x$repeat_interleave(0)
}




n_obs <- 10000
x_mat <- runif(n_obs * 2) %>% matrix(n_obs, 2)
beta_true <- matrix(c(1, 2), nrow=2, ncol=1)
y_vec <- x_mat %*% beta_true + 0.5 * rnorm(n_obs)
w_vec <- rep(1, n_obs) %>% matrix(ncol=1)



data_df <- data.frame(y=y_vec, x1=x_mat[,1], x2=x_mat[,2])
reg_fit <- lm(y ~ x1 + x2 - 1, data_df, x=TRUE, y=TRUE)
model_grads <- 
    zaminfluence::ComputeModelInfluence(reg_fit) %>%
    zaminfluence::AppendTargetRegressorInfluence("x1")


zam_time <- Sys.time()
grouped_model_grads <- 
    zaminfluence::ComputeModelInfluence(reg_fit, se_group=1:n_obs) %>%
    zaminfluence::AppendTargetRegressorInfluence("x1")
zam_time <- Sys.time() - zam_time


##############
# Torch

x <- torch_tensor(x_mat, requires_grad=FALSE)
y <- torch_tensor(y_vec, requires_grad=FALSE)
w <- torch_tensor(w_vec, requires_grad=TRUE)

xw <- x * w
xwx <- torch_matmul(torch_transpose(xw, 1, 2), x)
xwy <- torch_matmul(torch_transpose(xw, 1, 2), y)
betahat <- linalg_solve(xwx, xwy)

# COMPUTE GRAD
betahat_1 <- betahat[1]
betahat_1$backward()
betahat_1_grad <- as.numeric(w$grad)

# eps_vec <- y_vec - x_mat %*% as.numeric(betahat)
# xwx_mat <- t(x_mat * cbind(w_vec, w_vec)) %*% x_mat
# w_grad_closed <- solve(xwx_mat, t(x_mat * cbind(eps_vec, eps_vec))) %>% t()

#max(abs(betahat_1_grad - w_grad_closed[,1]))
max(abs(betahat_1_grad - model_grads$param_infls$x1$param$infl))



##############
# Ses

x <- torch_tensor(x_mat, requires_grad=FALSE)
y <- torch_tensor(y_vec, requires_grad=FALSE)
w <- torch_tensor(w_vec, requires_grad=TRUE)

xw <- x * w
xwx <- torch_matmul(torch_transpose(xw, 1, 2), x)
xwy <- torch_matmul(torch_transpose(xw, 1, 2), y)
betahat <- linalg_solve(xwx, xwy)

eps <- y - torch_matmul(x, betahat)
sig2hat <- torch_mean(w * (eps ** 2)) * (n_obs / (n_obs - 2))

betahat_cov <- torch_inverse(xwx) * sig2hat

# COMPUTE GRAD
betahat1_se <- torch_sqrt(betahat_cov[1, 1])
betahat1_se$backward()
betahat1_se_grad <- as.numeric(w$grad)

#max(abs(vcov(reg_fit) - as.numeric(betahat_cov))) # ok

max(abs(betahat1_se_grad - model_grads$param_infls$x1$se$infl))
#plot(betahat1_se_grad, model_grads$param_infls$x1$se$infl); abline(0, 1)




##############
# Grouped Ses

torch_time <- Sys.time()

x <- torch_tensor(x_mat, requires_grad=FALSE)
y <- torch_tensor(y_vec, requires_grad=FALSE)
w <- torch_tensor(w_vec, requires_grad=TRUE)

xw <- x * w
xwx <- torch_matmul(torch_transpose(xw, 1, 2), x)
xwy <- torch_matmul(torch_transpose(xw, 1, 2), y)
betahat <- linalg_solve(xwx, xwy)

eps <- y - torch_matmul(x, betahat)
score_mat  <- eps * x
score_cov <- torch_matmul(torch_transpose(w * score_mat, 2, 1), w * score_mat) / n_obs
xwx_inv_score_cov <- linalg_solve(xwx / n_obs, score_cov)
betahat_grouped_cov <- 
    linalg_solve(xwx / n_obs, torch_transpose(xwx_inv_score_cov, 1, 2)) / n_obs

betahat1_grouped_se <- torch_sqrt(betahat_grouped_cov[1, 1])

abs(betahat1_grouped_se - grouped_model_grads$param_infls$x1$se$base_value)


# COMPUTE GRAD
betahat1_grouped_se$backward()
betahat1_grouped_se_grad <- as.numeric(w$grad)

torch_time <- Sys.time() - torch_time


#max(abs(vcov(reg_fit) - as.numeric(betahat_cov))) # ok

max(abs(betahat1_grouped_se_grad - grouped_model_grads$param_infls$x1$se$infl))
if (FALSE) {
    plot(betahat1_grouped_se_grad, grouped_model_grads$param_infls$x1$se$infl); abline(0, 1)
}


torch_time
zam_time
