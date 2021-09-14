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

compare <- function(x, y) { return(max(abs(x - y))) }
check_equivalent  <- function(x, y) { stopifnot(compare(x, y) < 1e-8) }

n_obs <- 10000

set.seed(42)

RerunSummaryDf <- function(signals) {
    # A summary comparing reruns and predictions for each signal.
    rerun_df <-
        signals[c("sign", "sig", "both")] %>%
        map_dfr(~ .$rerun_df) %>%
        mutate(summary_orig=
                   sprintf("%f (%f, %f)", 
                           betahat_orig, beta_mzse_orig, beta_pzse_orig),
               summary_refit=
                   sprintf("%f (%f, %f)",
                           betahat_refit, beta_mzse_refit, beta_pzse_refit)) %>%
        select(description, num_removed, prop_removed, summary_orig, summary_refit)
    return(rerun_df)
}


#############################
# Oridinary regression.

# Generate data.
set.seed(42)
x_dim <- 3
beta_true <- 0.1 * runif(x_dim)
df <- GenerateRegressionData(n_obs, beta_true, num_groups=NULL)

# Fit a regression model.
x_names <- sprintf("x%d", 1:x_dim)
reg_form <- formula(sprintf("y ~ %s - 1", paste(x_names, collapse=" + ")))
model_fit <- lm(data = df, formula=reg_form, x=TRUE, y=TRUE)

# Get influence and reruns.
model_grads <- 
    ComputeModelInfluence(model_fit) %>%
    AppendTargetRegressorInfluence("x1")
signals <-
    GetInferenceSignals(model_grads$param_infl_list[["x1"]]) %>%
    RerunForTargetChanges(model_grads)

# Summaries comparing reruns and predictions for each signal.
RerunSummaryDf(signals)
PlotSignal(model_grads$param_infl_list[["x1"]], signals[["both"]], apip_max=0.03)

# Visualize which points are being dropped

df$drop <- FALSE
df$drop[signals[["both"]]$apip$inds] <- TRUE
df$infl <- model_grads$param_infl_list[["x1"]][[signals[["both"]]$qoi_name]]$infl

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
beta_true <- 0.1 * runif(x_dim)
df <- GenerateIVRegressionData(n_obs, beta_true, num_groups=NULL)

# Fit an IV model.
x_names <- sprintf("x%d", 1:x_dim)
z_names <- sprintf("z%d", 1:x_dim)
iv_form <- formula(sprintf("y ~ %s - 1 | %s - 1",
                           paste(x_names, collapse=" + "),
                           paste(z_names, collapse=" + ")))
model_fit <- ivreg(data = df, formula = iv_form, x=TRUE, y=TRUE)

# Get influence and reruns.
model_grads <- 
    ComputeModelInfluence(model_fit) %>%
    AppendTargetRegressorInfluence("x1")
signals <-
    GetInferenceSignals(model_grads$param_infl_list[["x1"]]) %>%
    RerunForTargetChanges(model_grads)

# Summaries comparing reruns and predictions for each signal.
RerunSummaryDf(signals)
PlotSignal(model_grads$param_infl_list[["x1"]], signals[["both"]], apip_max=0.03)



testthat::expect_equivalent(
    model_grads$parameter_names, model_fit$x %>% colnames(),
    info="column names")

model_fit$x$regressors


#############################
# Grouped standard errors.

# Generate data.
set.seed(42)
x_dim <- 3
beta_true <- 0.1 * runif(x_dim)
num_groups <- 50
df <- GenerateRegressionData(n_obs, beta_true, num_groups=num_groups)

# se_group is zero-indexed group indicator with no missing entries.
table(df$se_group)

# Fit a regression model.
x_names <- sprintf("x%d", 1:x_dim)
reg_form <- formula(sprintf("y ~ %s - 1", paste(x_names, collapse=" + ")))
model_fit <- lm(data=df, formula=reg_form, x=TRUE, y=TRUE)

# Get influence and reruns.
model_grads <- 
    ComputeModelInfluence(model_fit) %>%
    AppendTargetRegressorInfluence("x1")

signals <-
    GetInferenceSignals(model_grads$param_infl_list[["x1"]]) %>%
    RerunForTargetChanges(model_grads)

# Summaries comparing reruns and predictions for each signal.
RerunSummaryDf(signals)
grid.arrange(
    PlotSignal(model_grads$param_infl_list[["x1"]], signals[["sign"]], apip_max=0.03),
    PlotSignal(model_grads$param_infl_list[["x1"]], signals[["sig"]], apip_max=0.03),
    PlotSignal(model_grads$param_infl_list[["x1"]], signals[["both"]], apip_max=0.03),
    ncol=3
)


#############################
# Pairing

# In the current version of the refactor, pairing is not yet supported.
