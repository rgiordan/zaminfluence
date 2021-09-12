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


#############################
# Oridinary regression.

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
reg_signals <-
    GetRegressionSignals(model_grads$targets[["x1"]]) %>%
    RerunForTargetChanges(model_grads)

# A summary comparing reruns and predictions for each signal.
rerun_df <-
    reg_signals[c("sign", "sig", "both")] %>%
    map_dfr(~ .$rerun_df)
print(rerun_df)

PlotSignal(model_grads$targets[["x1"]], reg_signals[["both"]], apip_max=0.03)

# Visualize which points are being dropped

df$drop <- FALSE
df$drop[reg_signals[["both"]]$apip$inds] <- TRUE
df$infl <- model_grads$targets[["x1"]][[reg_signals[["both"]]$metric]]$infl

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

# Get influence.
iv_infl <- ComputeModelInfluence(iv_fit)
grad_df <- GetTargetRegressorGrads(iv_infl, "x1")
influence_dfs <- SortAndAccumulate(grad_df)

target_change <- GetRegressionTargetChange(influence_dfs, "prop_removed")
if (FALSE) {
    PlotInfluence(influence_dfs$sign, "prop_removed", 0.01, target_change)
}

rerun_df <- RerunForTargetChanges(influence_dfs, target_change, iv_fit)
select(rerun_df, change, beta, beta_pzse, beta_mzse, prop_removed)

