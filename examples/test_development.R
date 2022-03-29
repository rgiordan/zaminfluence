library(zaminfluence)
library(testthat)

set.seed(42)
x_dim <- 3
num_obs <- 10

param_true <- 1:x_dim / x_dim
df <- GenerateRegressionData(num_obs, param_true, num_groups=NULL)

# Fit a regression model.
x_names <- sprintf("x%d", 1:x_dim)
reg_form <- formula(sprintf("y ~ %s - 1", paste(x_names, collapse=" + ")))
fit_object <- lm(data=df, formula=reg_form, x=TRUE, y=TRUE)

# Check the behavior of keep_pars.
model_grads_full <- ComputeModelInfluence(fit_object, keep_pars=c("x1", "x2", "x3"))
model_grads_null <- ComputeModelInfluence(fit_object)
model_grads_2 <- ComputeModelInfluence(fit_object, keep_pars="x2")
model_grads_23 <- ComputeModelInfluence(fit_object, keep_pars=c("x2", "x3"))
model_grads_32 <- ComputeModelInfluence(fit_object, keep_pars=c("x3", "x2"))

TestKeepPars <- function(m, grad_pars, full_m=model_grads_full) {
    k <- length(grad_pars)
    
    # Make sure the kept parameter names are correct
    expect_equivalent(m$parameter_names, grad_pars)

    # The ModelFit should have the whole parameter vector
    expect_equivalent(m$model_fit$parameter_names, c("x1", "x2", "x3"))
    
    # The gradients should only be computed for the kept parameters
    # and should match the corresponding gradients in the full model.
    expect_equivalent(dim(m$param_grad), c(k, num_obs))
    expect_equivalent(dim(m$se_grad), c(k, num_obs))
    for (par in grad_pars) {
        i <- GetParameterIndex(m, par)
        full_i <- GetParameterIndex(full_m, par)
        expect_equivalent(m$param_grad[i, ], full_m$param_grad[full_i, ])
        expect_equivalent(m$se_grad[i, ], full_m$se_grad[full_i, ])
    }
}

TestKeepPars(model_grads_full, c("x1", "x2", "x3"))
TestKeepPars(model_grads_null, c("x1", "x2", "x3"))
TestKeepPars(model_grads_2, c("x2"))
TestKeepPars(model_grads_23, c("x2", "x3"))
TestKeepPars(model_grads_32, c("x3", "x2"))

# Check that it fails if you request a parameter that's not present
expect_error(AppendTargetRegressorInfluence(model_grads_23, "x4"), "x4 not found")
