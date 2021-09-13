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

n_obs <- 10000

set.seed(42)

git_repo_dir <- "/home/rgiordan/Documents/git_repos/zaminfluence"

library(devtools)
load_all("/home/rgiordan/Documents/git_repos/zaminfluence/zaminfluence")

#source(file.path(git_repo_dir, "zaminfluence/R/influence_lib.R"))


GenerateTestInstance <- function(do_iv, do_grouping) {
    x_dim <- 3
    beta_true <- runif(x_dim)
    num_obs <- 500
    
    GenerateFun <- if (do_iv) GenerateIVRegressionData else GenerateRegressionData
    if (do_grouping) {
        df <- GenerateFun(num_obs, beta_true, num_groups=5)
    } else {
        df <- GenerateFun(num_obs, beta_true)
    }
    
    df$weights <- runif(nrow(df)) + 1
    
    # Fit a model.
    if (do_iv) {
        # IV:
        x_names <- sprintf("x%d", 1:x_dim)
        z_names <- sprintf("z%d", 1:x_dim)
        reg_form <- formula(sprintf("y ~ %s - 1 | %s - 1",
                                    paste(x_names, collapse=" + "),
                                    paste(z_names, collapse=" + ")))
        model_fit <- ivreg(data=df, formula = reg_form,
                           x=TRUE, y=TRUE, weights=weights)
    } else {
        # Regression:
        x_names <- sprintf("x%d", 1:x_dim)
        reg_form <- formula(sprintf("y ~ %s - 1",
                                    paste(x_names, collapse=" + ")))
        model_fit <- lm(data=df, formula=reg_form,
                        x=TRUE, y=TRUE, weights=weights)
    }
    
    se_group <- if (do_grouping) df$se_group else NULL
    
    model_grads <-
        ComputeModelInfluence(model_fit) %>%
        AppendTargetRegressorInfluence("x1")
    signals <-
        GetInferenceSignals(model_grads$param_infl_list[["x1"]]) %>%
        RerunForTargetChanges(model_grads)
    
    return(list(
        model_grads=model_grads,
        signals=signals,
        model_fit=model_fit,
        se_group=se_group,
        df=df
    ))
}


test_instance <- GenerateTestInstance(FALSE, FALSE)


model_fit <- test_instance$model_fit
model_grads <- test_instance$model_grads
param_infl <- model_grads$param_infl_list[["x1"]]


# Sanity check an APIP.
# - Influence scores should match the sign
# - Cumulative influence scores should decreas or increase according to the sign
# - Lengths should match
TestAPIP <- function(qoi,  sign) {
    apip <- qoi[[sign]]
    apip_infl <- qoi$infl[apip$infl_inds]
    if (sign == "pos") {
        expect_true(all(apip_infl > 0), info="positive influence is positive")
        expect_true(all(diff(apip$infl_cumsum) > 0), info="positive influence is increasing")
    } else if (sign == "neg") {
        expect_true(all(apip_infl < 0), info="negative influence is negative")
        expect_true(all(diff(apip$infl_cumsum) < 0), info="negative influence is decreasing")
    } else {
        stop("Bad sign passed to test")
    }
    expect_true(length(apip$infl_inds) == length(apip$infl_cumsum), "Inds len == cumsum len")    
}


# Check the validity of a prediction when leaving out a small number of influential points.
# We check that the relative error in the difference is less than 100 * tol %.
TestPredictions <- function(param_infl, qoi_name, sign, num_leave_out=2, tol=0.03) {
    # Get the indices to drop.
    qoi <- param_infl[[qoi_name]]
    apip <- qoi[[sign]]
    drop_inds <- apip$infl_inds[1:num_leave_out]
    w_bool <- GetWeightVector(drop_inds, num_obs=model_grads$n_obs, bool=TRUE)
    
    # The original values
    base_values <- GetBaseValues(param_infl)

    # Rerun
    rerun <- model_grads$RerunFun(model_grads$model_fit, w_bool)
    rerun_base_values <- GetRerunBaseValues(rerun, param_infl)
    diff_rerun <- rerun_base_values - base_values
    
    # Prediction
    diff_pred <- 
        map_dbl(names(base_values), 
                ~ PredictChange(param_infl[[.]], drop_inds))
    
    # Check the maximum relative error amongst beta, beta_mzse, and beta_pzse
    max_rel_err <- max(abs((diff_pred - diff_rerun) / diff_rerun))
    expect_true(max_rel_err < tol, 
                info=sprintf("%s %s prediction error: %f", 
                             qoi_name, sign, max_rel_err))
}


TestSignalPrediction <- function(param_infl, signal_name) {
    signal <- signals[[signal_name]]
    qoi_name <- signal$qoi_name
    
    # Form the prediction
    drop_inds <- signal$apip$inds
    base_value <- GetBaseValues(param_infl)[qoi_name]
    pred_diff <- PredictChange(param_infl[[qoi_name]], drop_inds)
    pred_value <- base_value + pred_diff
    
    # Assert that a sign change took place.
    # Everything we look for is a sign change.
    expect_true(
        sign(base_value) != sign(pred_value),
        sprintf("%s predicted to acheive sign change",
                signal_name))
    
    if (length(drop_inds) > 1) {
        pred_diff <- PredictChange(param_infl[[qoi_name]], 
                                   drop_inds[1:(length(drop_inds) - 1)])
        pred_value <- base_value + pred_diff
        expect_true(
            sign(base_value) == sign(base_value + pred_diff),
            sprintf("%s predicted to fail to sign change with one fewer point",
                    signal_name))
        
    }
}



# Check the validity of the influence scores.
qoi_names <- c("beta", "beta_mzse", "beta_pzse")
for (qoi_name in qoi_names) {
    for (sign in c("pos", "neg")) {
        TestAPIP(param_infl[[qoi_name]], sign)
        TestPredictions(param_infl, qoi_name, sign)
    }
}

# Check that the APIP predicts the appropriate change, and that
# one fewer point does not.
signals <- GetInferenceSignals(param_infl)
for (signal_name in c("sign", "sig", "both")) {
    TestSignalPrediction(param_infl, signal_name)

    # Just test that this runs.
    GetSignalDataFrame(signals[[signal_name]])
}



load_all("/home/rgiordan/Documents/git_repos/zaminfluence/zaminfluence")

# Test GetAMIS
qoi <- param_infl$beta_mzse
n_drops <- c(0, 1, 10, 10000)
for (sign in c("pos", "neg")) {
    for (n_drop in n_drops) {
        suppressWarnings(amis <- GetAMIS(qoi, sign=sign, n_drop=n_drop))
        if (n_drop == 0) {
            expect_equivalent(amis, NULL)        
        } else if (n_drop == 1) {
            expect_equivalent(amis, qoi[[sign]]$infl_inds[1])
        } else if (n_drop == 10) {
            expect_equivalent(amis, qoi[[sign]]$infl_inds[1:10])
        } else if (n_drop == 10000) {
            expect_equivalent(amis, qoi[[sign]]$infl_inds)
        } else {
            stop(sprintf("Bad number of drops: %d", n_drop))
        }
    }
}


##########################################################
##########################################################
##########################################################
##########################################################
# Work in progress:





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

