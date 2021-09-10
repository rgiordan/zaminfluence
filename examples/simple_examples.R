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
    RerunForTargetChanges(model_grads, RerunFun=RerunRegression)

rerun_df <- 
    lapply(c("sign", "sig", "both"), 
           function(x) { reg_signals[[x]]$rerun_df }) %>%
    do.call(bind_rows, .)

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
