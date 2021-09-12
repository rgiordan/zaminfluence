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
#library(zaminfluence)
library(purrr)
library(AER)

compare <- function(x, y) { return(max(abs(x - y))) }
check_equivalent  <- function(x, y) { stopifnot(compare(x, y) < 1e-8) }

n_obs <- 10000

set.seed(42)


library(devtools)
load_all("/home/rgiordan/Documents/git_repos/zaminfluence/zaminfluence")

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
        select(change, num_removed, summary_orig, summary_refit)
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
reg_fit <- lm(data = df, formula = reg_form, x=TRUE, y=TRUE)

# Get influence and reruns.
model_grads <- 
    ComputeModelInfluence(reg_fit) %>%
    AppendTargetRegressorInfluence("x1")
signals <-
    GetInferenceSignals(model_grads$targets[["x1"]]) %>%
    RerunForTargetChanges(model_grads)

# Summaries comparing reruns and predictions for each signal.
RerunSummaryDf(signals)
PlotSignal(model_grads$targets[["x1"]], signals[["both"]], apip_max=0.03)

# Visualize which points are being dropped

df$drop <- FALSE
df$drop[signals[["both"]]$apip$inds] <- TRUE
df$infl <- model_grads$targets[["x1"]][[signals[["both"]]$metric]]$infl

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
iv_fit <- ivreg(data = df, formula = iv_form, x=TRUE, y=TRUE)

# Get influence and reruns.
model_grads <- 
    ComputeModelInfluence(iv_fit) %>%
    AppendTargetRegressorInfluence("x1")
signals <-
    GetInferenceSignals(model_grads$targets[["x1"]]) %>%
    RerunForTargetChanges(model_grads)

# Summaries comparing reruns and predictions for each signal.
RerunSummaryDf(signals)
PlotSignal(model_grads$targets[["x1"]], signals[["both"]], apip_max=0.03)

