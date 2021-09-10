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


reg_infl$targets[["x1"]] %>% names()

target_change_df <- 
    lapply(reg_signals, GetSignalDataFrame) %>%
    do.call(bind_rows, .)


rerun_df <- data.frame()
for (reg_signal in reg_signals) {
    rerun_df <- bind_rows(
        rerun_df,
    )
}

rerun_df
