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

# Get influence.
reg_infl <- ComputeModelInfluence(reg_fit)
grad_df <- GetTargetRegressorGrads(reg_infl, "x1")
influence_dfs <- SortAndAccumulate(grad_df)

# Compute the AMIP and friends for changes of sign and significance.
target_change <- GetRegressionTargetChange(influence_dfs, "prop_removed")

if (FALSE) {
    PlotInfluence(influence_dfs$sign, "prop_removed", 0.01, target_change)
}

# Rerun.
rerun_df <- RerunForTargetChanges(influence_dfs, target_change, reg_fit)
select(rerun_df, change, beta, beta_pzse, beta_mzse, prop_removed)


# See which points were left out for, e.g., a sign change.
target_change_row <- filter(target_change, change == "sign")
w_drop <- GetWeightForTargetChangeRow(influence_dfs, target_change_row, rows_to_keep=FALSE)

drop_df <- data.frame(
    leverage=hatvalues(reg_fit),
    residual=reg_fit$residuals,
    x1=reg_fit$x[, 1],
    x2=reg_fit$x[, 2],
    drop=w_drop)

if (FALSE) {
    # The beta gradient is not precisely the leverage score times the residual,
    # but it's closely related.  See the paper for more details.
    ggplot(drop_df) +
        geom_point(aes(x=leverage, y=residual)) +
        geom_point(aes(x=leverage, y=residual), data=filter(drop_df, drop),
                   color="red", size=3)
}

# To compare an influence score to the original data, use the row column of
# an influence dataframe.
infl_df <- zaminfluence::GetInflDfForTargetChangeRow(influence_dfs, target_change_row)

drop_infl_df <- 
    drop_df %>%
    mutate(row=1:nrow(drop_df)) %>%
    inner_join(select(infl_df, row, beta_grad), by="row")

if (FALSE) {
    # Sanity check that we're dropping the points with the largest gradients.
    ggplot(drop_infl_df) +
        geom_histogram(aes(x=beta_grad, y=..density.., fill=drop, group=drop), alpha=0.5)
    
    # You can then examine the resulation between the influence scores and the
    # values of the regressors, for example.
    grid.arrange(
        ggplot(drop_infl_df) +
            geom_point(aes(x=x1, y=beta_grad)),
        ggplot(drop_infl_df) +
            geom_point(aes(x=x2, y=beta_grad)),
        ncol=2
    )
}


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


#############################
# Pairing

x_dim <- 3
beta_true <- 0.1 * runif(x_dim)
df <- GenerateRegressionData(n_obs, beta_true, num_groups=NULL)
df$arm <- as.integer(runif(nrow(df)) < 0.5)  # Control or treatment
df$group <- sample.int(3, n_obs, replace=TRUE)  # Some subgrouping
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
select(rerun_df, change, beta, beta_pzse, beta_mzse, prop_removed)



#############################
# Aggregated pairing

x_dim <- 3
beta_true <- 0.1 * runif(x_dim)
df <- GenerateRegressionData(n_obs, beta_true, num_groups=NULL)
df$arm <- as.integer(runif(nrow(df)) < 0.5)  # Control or treatment
df$group <- sample.int(floor(n_obs / 10), n_obs, replace=TRUE)  # Some subgrouping
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
              beta_grad=sum(beta_grad),
              beta_pzse_grad=sum(beta_pzse_grad),
              beta_mzse_grad=sum(beta_mzse_grad),
              obs_per_row=sum(obs_per_row),
              .groups="drop") %>%
    mutate(grouped_row=1:n()) %>%
    CopyGradAttributes(base_grad_df)

# There is no longer a column that matches rows to the original data.
attr(grad_df, "data_row_cols") <- "grouped_row"

# Note that I still consider one of the original rows to be an observation,
# so the n_obs attribute remains the same.  If a single group were to be
# an "observation", then you'd want ``obs_per_row`` to contain all ones,
# and the ``n_obs`` attribute to be set to the number of unique groups.
print(attr(grad_df, "n_obs"))

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
grad_df[w_bool, ] %>% select(arm, group, beta_grad, grouped_row)
table(grad_df[w_bool, "arm"])


#############################
# Grouped standard errors.

x_dim <- 3
beta_true <- runif(x_dim)
num_groups <- 50
df <- GenerateRegressionData(n_obs, beta_true, num_groups=num_groups)

# se_group is zero-indexed group indicator with no missing entries.
table(df$se_group)

# Fit a regression model.
x_names <- sprintf("x%d", 1:x_dim)
reg_form <- formula(sprintf("y ~ %s - 1", paste(x_names, collapse=" + ")))
reg_fit <- lm(data = df, formula = reg_form, x=TRUE, y=TRUE)

# Get influence.
reg_infl <- ComputeModelInfluence(reg_fit, se_group=df$se_group)
grad_df <- GetTargetRegressorGrads(reg_infl, "x1")
influence_dfs <- SortAndAccumulate(grad_df)

# Here's what our grouped standard errors match.
se_cov <- vcovCL(reg_fit, cluster=df$se_group, type="HC0", cadjust=FALSE)
base_vals <- attr(grad_df, "base_vals")
target_index <- attr(grad_df, "target_index")
cat(
    sprintf("vcovCL:\t\t%f", sqrt(se_cov[target_index, target_index])), "\n",
    sprintf("zaminfluence:\t%f", base_vals["se"]), "\n", sep="")

target_change <- GetRegressionTargetChange(influence_dfs, "prop_removed")
if (FALSE) {
    PlotInfluence(influence_dfs$sign, "prop_removed", 0.01, target_change)
}

rerun_df <- RerunForTargetChanges(influence_dfs, target_change, reg_fit, se_group=df$se_group)
select(rerun_df, change, beta, beta_pzse, beta_mzse, prop_removed)


