library(torch)
library(tidyverse)
library(zaminfluence)
library(AER)

library(devtools)
#devtools::load_all("/home/rgiordan/Documents/git_repos/zaminfluence/zaminfluence")
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

model_grads <- 
    zaminfluence::ComputeModelInfluence(iv_res, se_group=se_group, keep_pars=c("x2", "x1")) %>%
    zaminfluence::AppendTargetRegressorInfluence("x1")

GetParameterIndex(model_grads$model_fit, "x1")
GetParameterIndex(model_grads, "x1")

rownames(model_grads$param_grad)

