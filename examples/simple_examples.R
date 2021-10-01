###################################################################
#
# These simple examples illustrate the use of zaminfluence.
# https://github.com/rgiordan/zaminfluence
# See the README.md file for installation instructions.

library(ggplot2)
library(dplyr)
library(reshape2)
library(gridExtra)
library(sandwich)
library(zaminfluence)
library(purrr)
library(AER)


library(devtools)
load_all("/home/rgiordan/Documents/git_repos/zaminfluence/zaminfluence")

compare <- function(x, y) { return(max(abs(x - y))) }
check_equivalent  <- function(x, y) { stopifnot(compare(x, y) < 1e-8) }

n_obs <- 10000

set.seed(42)

SummarizeReruns <- function(reruns, preds) {
    reruns_df <- GetSignalsAndRerunsDataframe(signals, reruns, model_grads)
    preds_df <- GetSignalsAndRerunsDataframe(signals, preds, model_grads)
    
    summary_df <-
        rbind(reruns_df %>% mutate(method="rerun"),
              preds_df %>% mutate(method="prediction")) %>%
        pivot_wider(-method, names_from=method, values_from=value)
    return(summary_df)
}


#############################
# Oridinary regression.

load_all("/home/rgiordan/Documents/git_repos/zaminfluence/zaminfluence")

# Generate data.
set.seed(42)
x_dim <- 3
param_true <- 0.1 * runif(x_dim)
df <- GenerateRegressionData(n_obs, param_true, num_groups=NULL)

# Fit a regression model.
x_names <- sprintf("x%d", 1:x_dim)
reg_form <- formula(sprintf("y ~ %s - 1", paste(x_names, collapse=" + ")))
fit_object <- lm(data = df, formula=reg_form, x=TRUE, y=TRUE)

# Get influence and reruns.
model_grads <-
    ComputeModelInfluence(fit_object) %>%
    AppendTargetRegressorInfluence("x1") %>%
    AppendTargetRegressorInfluence("x2")

signals <- GetInferenceSignals(model_grads)
reruns <- RerunForSignals(signals, model_grads)
preds <- PredictForSignals(signals, model_grads)
summary_df <- SummarizeReruns(reruns, preds)

ggplot(summary_df) +
    geom_point(aes(x=prediction, y=rerun, color=param_name, shape=metric)) +
    geom_abline(aes(slope=1, intercept=0))

PlotSignal(model_grads, signals, "x1", "sign",
          reruns=reruns, apip_max=0.03)


# Visualize which points are being dropped
signal <- signals[["x1"]][["both"]] 
df$drop <- GetWeightVector(signal$apip$inds, nrow(df), bool=TRUE, invert=TRUE)
df$infl <- signal$qoi$infl

grid.arrange(
    ggplot(df) +
        geom_point(aes(x=x1, y=infl, color=drop)),
    ggplot(df) +
        geom_point(aes(x=x1, y=y, color=drop)),
    ggplot(df) +
        geom_point(aes(x=x1, y=x2, color=drop)),
    ncol=3
)






#############################
# Instrumental variables.

# Generate data.
set.seed(42)
x_dim <- 3
param_true <- 0.1 * runif(x_dim)
df <- GenerateIVRegressionData(n_obs, param_true, num_groups=NULL)

# Fit an IV model.
x_names <- sprintf("x%d", 1:x_dim)
z_names <- sprintf("z%d", 1:x_dim)
iv_form <- formula(sprintf("y ~ %s - 1 | %s - 1",
                           paste(x_names, collapse=" + "),
                           paste(z_names, collapse=" + ")))
fit_object <- ivreg(data = df, formula = iv_form, x=TRUE, y=TRUE)

# Get influence and reruns.
model_grads <-
    ComputeModelInfluence(fit_object) %>%
    AppendTargetRegressorInfluence("x1")

signals <- GetInferenceSignals(model_grads)
reruns <- RerunForSignals(signals, model_grads)
preds <- PredictForSignals(signals, model_grads)

summary_df <- SummarizeReruns(reruns, preds)
ggplot(summary_df) +
    geom_point(aes(x=prediction, y=rerun, color=param_name, shape=metric)) +
    geom_abline(aes(slope=1, intercept=0))

PlotSignal(model_grads, signals, "x1", "sign",
           reruns=reruns, apip_max=0.03)


#############################
# Grouped standard errors.

# Generate data.
set.seed(42)
x_dim <- 3
param_true <- 0.1 * runif(x_dim)
num_groups <- 50
df <- GenerateRegressionData(n_obs, param_true, num_groups=num_groups)

# se_group is zero-indexed group indicator with no missing entries.
table(df$se_group)

# Fit a regression model.
x_names <- sprintf("x%d", 1:x_dim)
reg_form <- formula(sprintf("y ~ %s - 1", paste(x_names, collapse=" + ")))
fit_object <- lm(data=df, formula=reg_form, x=TRUE, y=TRUE)

# Get influence and reruns.
model_grads <-
    ComputeModelInfluence(fit_object) %>%
    AppendTargetRegressorInfluence("x1")


signals <- GetInferenceSignals(model_grads)
reruns <- RerunForSignals(signals, model_grads)
preds <- PredictForSignals(signals, model_grads)
summary_df <- SummarizeReruns(reruns, preds)

ggplot(summary_df) +
    geom_point(aes(x=prediction, y=rerun, color=param_name, shape=metric)) +
    geom_abline(aes(slope=1, intercept=0))

# Summaries comparing reruns and predictions for each signal.
plots <-
    map(c("sign", "sig", "both"),
        ~ PlotSignal(model_grads, signals, "x1", ., reruns=reruns, apip_max=0.03))
do.call(function(...) { grid.arrange(..., ncol=1) }, plots)


#############################
# Pairing

# In the current version of the refactor, pairing is not yet supported.
