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

git_repo_dir <- "/home/rgiordan/Documents/git_repos/zaminfluence"

#source(file.path(git_repo_dir, "zaminfluence/R/influence_lib.R"))


#############################
# Oridinary regression.

x_dim <- 3
beta_true <- 0.1 * runif(x_dim)
df <- GenerateRegressionData(n_obs, beta_true, num_groups=NULL)

# Fit a regression model.
x_names <- sprintf("x%d", 1:x_dim)
reg_form <- formula(sprintf("y ~ %s - 1", paste(x_names, collapse=" + ")))
reg_fit <- lm(data = df, formula = reg_form, x=TRUE, y=TRUE)

# Get influence.
reg_infl <- ComputeModelInfluence(reg_fit)
reg_infl <- AppendTargetRegressorInfluence(reg_infl, "x1")
reg_signals <-
    GetRegressionSignals(reg_infl$targets[["x1"]]) %>%
    RerunForTargetChanges(reg_infl)

# Make the reruns.
rerun_df <- 
    lapply(c("sign", "sig", "both"), 
           function(x) { reg_signals[[x]]$rerun_df }) %>%
    do.call(bind_rows, .)


parameter_infl <- reg_infl$targets[["x1"]]
names(parameter_infl)
sorting_qoi_name <- "beta"
qoi_infl <- parameter_infl[[sorting_qoi_name]]
qoi_infl %>% names()
qoi_infl$neg %>% names()
plot_num_dropped <- FALSE

##################

names(reg_signals[["sign"]])

influence_df <- GetSortedInfluenceDf(parameter_infl, "beta")
PlotInfluence(influence_df, reg_signals=list(reg_signals[["sign"]]), apip_max=0.05)

reg_signal <- reg_signals[["sign"]]
names(reg_signal$apip)

PlotSignal(parameter_infl, reg_signals[["both"]], apip_max=0.03)

