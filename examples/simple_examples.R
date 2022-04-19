###################################################################
#
# These simple examples illustrate the use of zaminfluence.
# https://github.com/rgiordan/zaminfluence
#
# See the README.md file for installation instructions.
#
# These examples show only how to run zaminfluence in different configurations.
# For more examples of how to interpret and analyze the output, see
# the script interpreting_output.R in this directory.

library(tidyverse)
library(gridExtra)
library(zaminfluence)
library(AER)

compare <- function(x, y) { return(max(abs(x - y))) }
check_equivalent  <- function(x, y) { stopifnot(compare(x, y) < 1e-8) }

num_obs <- 10000

set.seed(42)


#############################
# Oridinary regression.

# Generate data.
set.seed(42)
x_dim <- 3
param_true <- 0.1 * runif(x_dim)
df <- GenerateRegressionData(num_obs, param_true, num_groups=NULL)

# Fit a regression model.
x_names <- sprintf("x%d", 1:x_dim)
reg_form <- formula(sprintf("y ~ %s - 1", paste(x_names, collapse=" + ")))
fit_object <- lm(data = df, formula=reg_form, x=TRUE, y=TRUE)

# Get influence and reruns.
# Derivatives are only computed for keep_pars, which can be in any order.
model_grads <-
    ComputeModelInfluence(fit_object, keep_pars=c("x2", "x1")) %>%
    AppendTargetRegressorInfluence("x1") %>%
    AppendTargetRegressorInfluence("x2")

signals <- GetInferenceSignals(model_grads)
reruns <- RerunForSignals(signals, model_grads)
preds <- PredictForSignals(signals, model_grads)



#############################
# Instrumental variables.

# Generate data.
set.seed(42)
x_dim <- 3
param_true <- 0.1 * runif(x_dim)
df <- GenerateIVRegressionData(num_obs, param_true, num_groups=NULL)

# Fit an IV model.
x_names <- sprintf("x%d", 1:x_dim)
z_names <- sprintf("z%d", 1:x_dim)
iv_form <- formula(sprintf("y ~ %s - 1 | %s - 1",
                           paste(x_names, collapse=" + "),
                           paste(z_names, collapse=" + ")))
fit_object <- ivreg(data = df, formula = iv_form, x=TRUE, y=TRUE)

# Get influence and reruns.
model_grads <-
    ComputeModelInfluence(fit_object, keep_pars=c("x2", "x1")) %>%
    AppendTargetRegressorInfluence("x1")

signals <- GetInferenceSignals(model_grads)
reruns <- RerunForSignals(signals, model_grads)
preds <- PredictForSignals(signals, model_grads)


#############################
# Grouped standard errors.

# Generate data.
set.seed(42)
x_dim <- 3
param_true <- 0.1 * runif(x_dim)
num_groups <- 50
df <- GenerateRegressionData(num_obs, param_true, num_groups=num_groups)

# se_group is zero-indexed group indicator with no missing entries.
table(df$se_group)

# Fit a regression model.
x_names <- sprintf("x%d", 1:x_dim)
reg_form <- formula(sprintf("y ~ %s - 1", paste(x_names, collapse=" + ")))
fit_object <- lm(data=df, formula=reg_form, x=TRUE, y=TRUE)

# Get influence and reruns.  Pass the grouping indicator to the `se_group`` argument
# of `ComputeModelInfluence`.
model_grads <-
    ComputeModelInfluence(fit_object, se_group=df$se_group) %>%
    AppendTargetRegressorInfluence("x1")


# The grouped standard error which zaminfluence computes...
cat("Zaminfluence SE:\t", model_grads$param_infls[["x1"]]$se$base_value, "\n")

# ...is equivalent to the that computed by the following standard command:
cat("vcovCL se:\t\t", 
    vcovCL(fit_object, cluster=df$se_group, type="HC0", cadjust=FALSE)["x1", "x1"] %>% sqrt(), 
    "\n")

signals <- GetInferenceSignals(model_grads)
reruns <- RerunForSignals(signals, model_grads)
preds <- PredictForSignals(signals, model_grads)


################################################
# Customizing some of what zaminfluence does


# Generate data.
set.seed(42)
x_dim <- 3
param_true <- 0.1 * runif(x_dim)
df <- GenerateRegressionData(num_obs, param_true, num_groups=NULL)

# Fit a regression model.
x_names <- sprintf("x%d", 1:x_dim)
reg_form <- formula(sprintf("y ~ %s - 1", paste(x_names, collapse=" + ")))
fit_object <- lm(data = df, formula=reg_form, x=TRUE, y=TRUE)

# Get influence and reruns.
model_grads <-
    ComputeModelInfluence(fit_object) %>%
    AppendTargetRegressorInfluence("x1") %>%
    AppendTargetRegressorInfluence("x2")

# By default, RerunForSignals re-runs the model for all parameters and 
# all quantities of interest.  You can also manually pick out a single
# signal to rerun.

# signals is a nested list of parameters and quantities of interest.
signals <- GetInferenceSignals(model_grads)
signal <- signals[["x1"]][["sign"]]

if (signal$apip$success) {
    cat("Rerunning for ", signal$description, ".\n", sep="")
    weights <- GetWeightVector(drop_inds=signal$apip$inds, 
                               orig_weights=model_grads$model_fit$weights)
    rerun <- model_grads$RerunFun(weights)
    pred <- PredictModelFit(model_grads, weights)
} else {
    cat("The linear approximation cannot reverse the signal  ", 
         signal$description, "; skipping rerun.\n", sep="")
}

cbind(coefficients(fit_object), rerun$param , pred$param)


# By default, RerunFun uses a weighted regression where some weights
# are set to zero.  You can also write your own rerun function that
# manually drops the rows.  This can help when dropping causes colinearity
# due to, say, eliminating some levels of a fixed effect indicator.

# A rerun function must take model weights and return a ModelFit object.
# The default is in model_grads$RerunFun, which you can use as a template.
CustomRerunFun <- function(weights) {
    keep_rows <- abs(weights) > 1e-8
    df_drop <- df[keep_rows, ]
    # We don't need x=TRUE and y=TRUE because we won't compute gradients from the
    # fit object with dropped rows.
    fit_object_drop <- lm(data=df_drop, formula=reg_form)
    
    model_fit_drop <- ModelFit(
        # The fit object isn't that important for a rerun.
        fit_object=fit_object_drop,
        
        # In the default RerunFun, the num_obs is the original number
        # not the number after dropping.
        num_obs=length(weights),
        
        param=coefficients(fit_object_drop),
        se=vcov(fit_object_drop) %>% diag() %>% sqrt(), 
        parameter_names=names(coefficients(fit_object_drop)), 
        weights=weights,
        se_group=model_grads$model_fit$se_group)
    return(model_fit_drop)
}

rerun_v2 <- CustomRerunFun(weights)
cbind(coefficients(fit_object), rerun$param , rerun_v2$param, pred$param)
