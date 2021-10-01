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
library(testthat)

#library(zaminfluence)
library(AER)

compare <- function(x, y) { return(max(abs(x - y))) }
check_equivalent  <- function(x, y) { stopifnot(compare(x, y) < 1e-8) }

num_obs <- 10000

set.seed(42)

git_repo_dir <- "/home/rgiordan/Documents/git_repos/zaminfluence"

library(devtools)
load_all("/home/rgiordan/Documents/git_repos/zaminfluence/zaminfluence")


source(file.path(git_repo_dir, "zaminfluence/tests/testthat/test_derivs.R"))



##########################################################
##########################################################
##########################################################
##########################################################
# Work in progress:



############################################
# Demonstration of unnesting with indices

foo <- list()
foo[["a"]] <- list()
foo[["b"]] <- list()

foo[["a"]][["x"]] <- data.frame(x=runif(1), z=1)
foo[["b"]][["x"]] <- data.frame(x=runif(1), z=2)
foo[["a"]][["y"]] <- data.frame(x=runif(1), z=3)
foo[["b"]][["y"]] <- data.frame(x=runif(1), z=4)
foo[["a"]][["z"]] <- data.frame(x=runif(1), z=5, no=10)

# This appears to work
tibble(list=foo) %>%
    mutate(letter=names(list)) %>%
    unnest_longer(list) %>%
    mutate(other_letter=names(list)) %>%
    unnest(list)


# This also appears to work
tibble(list=foo) %>%
    mutate(letter=names(list)) %>%
    unnest_longer(list) %>%
    mutate(other_letter=names(list)) %>%
    mutate(new_df=map(list, ~ data.frame(foo=.x$z))) %>%
    unnest(new_df) %>%
    select(-list)





#############################
# Pairing

x_dim <- 3
param_true <- 0.1 * runif(x_dim)
df <- GenerateRegressionData(num_obs, param_true, num_groups=NULL)
df$arm <- as.integer(runif(nrow(df)) < 0.5)  # Control or treatment
df$group <- sample.int(3, num_obs, replace=TRUE)  # Some subgrouping
df$row <- 1:nrow(df)

# Fit a regression model.
x_names <- sprintf("x%d", 1:x_dim)
reg_form <- formula(sprintf("y ~ %s + arm", paste(x_names, collapse=" + ")))
reg_fit <- lm(data = df, formula = reg_form, x=TRUE, y=TRUE)

# Get influence, where you remove arm pairs within groups.
reg_infl <- ComputeModelInfluence(reg_fit)
grad_df <-
    GetTargetRegressorGrads(reg_infl, "x1") %>%
    bind_cols(df[c("arm", "group")])

influence_dfs <-
    PairInfluenceScores(grad_df,
                        assignment_col="arm",
                        group_cols="group",
                        level0=0,
                        level1=1)

target_change <- GetRegressionTargetChange(influence_dfs, "prop_removed")

if (FALSE) {
    PlotInfluence(influence_dfs$sign, "prop_removed", 0.01, target_change)
}

# Rerun
rerun_df <- RerunForTargetChanges(influence_dfs, target_change, reg_fit)
select(rerun_df, change, param, param_pzse, param_mzse, prop_removed)



#############################
# Aggregated pairing

x_dim <- 3
param_true <- 0.1 * runif(x_dim)
df <- GenerateRegressionData(num_obs, param_true, num_groups=NULL)
df$arm <- as.integer(runif(nrow(df)) < 0.5)  # Control or treatment
df$group <- sample.int(floor(num_obs / 10), num_obs, replace=TRUE)  # Some subgrouping
df$row <- 1:nrow(df)

# Fit a regression model.
x_names <- sprintf("x%d", 1:x_dim)
reg_form <- formula(sprintf("y ~ %s + arm", paste(x_names, collapse=" + ")))
reg_fit <- lm(data = df, formula = reg_form, x=TRUE, y=TRUE)

# Get influence, where you aggregate arm pairs aggregated within a group.
reg_infl <- ComputeModelInfluence(reg_fit)
base_grad_df <- GetTargetRegressorGrads(reg_infl, "x1")
grad_df <-
    base_grad_df %>%
    bind_cols(df[c("arm", "group")]) %>%
    group_by(arm, group) %>%
    summarize(se_grad=sum(se_grad),
              param_grad=sum(param_grad),
              param_pzse_grad=sum(param_pzse_grad),
              param_mzse_grad=sum(param_mzse_grad),
              obs_per_row=sum(obs_per_row),
              .groups="drop") %>%
    mutate(grouped_row=1:n()) %>%
    CopyGradAttributes(base_grad_df)

# There is no longer a column that matches rows to the original data.
attr(grad_df, "data_row_cols") <- "grouped_row"

# Note that I still consider one of the original rows to be an observation,
# so the num_obs attribute remains the same.  If a single group were to be
# an "observation", then you'd want ``obs_per_row`` to contain all ones,
# and the ``num_obs`` attribute to be set to the number of unique groups.
print(attr(grad_df, "num_obs"))

# However, we still need to record how many distinct gradients are being
# sorted so we can index into grad_df using our results.
attr(grad_df, "n_grad_rows") <- nrow(grad_df)

head(grad_df)

influence_dfs <-
    PairInfluenceScores(grad_df,
                        assignment_col="arm",
                        group_cols="group",
                        level0=0,
                        level1=1)

target_change <- GetRegressionTargetChange(influence_dfs, "prop_removed")

if (FALSE) {
    PlotInfluence(influence_dfs$sign, "prop_removed", 0.01, target_change)
}

# Note that rerun does not work with aggregated influence functions, because
# it does not know which rows of the original data correspond to the rows of
# the influence function.

# However, we can still get weights into the grad_df data frame.
this_change <- filter(target_change, change == "sign")
infl_df <- influence_dfs[["sign"]][[this_change[["direction"]]]]
w_bool <- GetWeightForAlpha(infl_df, "prop_removed",
                            this_change[["prop_removed"]], rows_to_keep=FALSE)
grad_df[w_bool, ] %>% select(arm, group, param_grad, grouped_row)
table(grad_df[w_bool, "arm"])

