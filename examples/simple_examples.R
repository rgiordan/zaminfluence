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
fit_object <- lm(data = df, formula=reg_form, x=TRUE, y=TRUE)

# Get influence and reruns.
model_grads <-
    ComputeModelInfluence(fit_object) %>%
    AppendTargetRegressorInfluence("x1") %>%
    AppendTargetRegressorInfluence("x2")

signals <- GetInferenceSignals(model_grads)




#############################################

RerunForSignals <- function(signals, model_grads, RerunFun=NULL, verbose=FALSE) {
    stopifnot(class(model_grads) == "ModelGrads")
    stopifnot(setequal(names(signals), names(model_grads$param_infls)))
    for (param_name in names(signals)) {
        for (signal in signals[[param_name]]) {
            stopifnot(class(signal) == "QOISignal")
        }
    }

    if (is.null(RerunFun)) {
        RerunFun <- model_grads$RerunFun
    }
    
    verbosePrint <- function(...) {
        if (verbose) cat(..., sep="")
    }
    reruns <- list()
    for (param_name in names(signals)) {
        verbosePrint("Re-running for ", param_name, " signals: ")
        reruns[[param_name]] <- list()
        param_signals <- signals[[param_name]]
        #pram_infl <- model_grads$param_infls[[param_name]]
        for (signal_name in names(param_signals)) {
            verbosePrint(signal_name, "...")
            signal <- param_signals[[signal_name]]
            w_bool <- GetWeightVector(
                drop_inds=signal$apip$inds,
                num_obs=model_grads$model_fit$n_obs,
                bool=TRUE)
            reruns[[param_name]][[signal_name]] <- RerunFun(w_bool)
        }
        verbosePrint("done.\n")
    }
    return(reruns)
}

reruns <- RerunForSignals(signals, model_grads, verbose=TRUE)


library(tibble)
map_depth(signals, 2, ~ data.frame(signal=.$signal, row.names=NULL))


# Alternative to RerunForSignals
RerunFun <- model_grads$RerunFun
num_obs <- model_grads$model_fit$n_obs
RerunSignal <- function(signal) {
    cat("Rerunning ", signal$description, "\n", sep="")
    w_bool <- GetWeightVector(
        drop_inds=signal$apip$inds,
        num_obs=num_obs,
        bool=TRUE)
    return(RerunFun(w_bool))
}
reruns_v2 <- map_depth(signals, 2, RerunSignal)

names(reruns_v2)
names(reruns)

reruns_v2[[1]][[2]]$betahat
reruns[[1]][[2]]$betahat

reruns_v2[[2]][[2]]$betahat
reruns[[2]][[2]]$betahat

# A compact way to make dataframes
reruns_dfs <- map_depth(reruns, 2, ~ GetModelFitInferenceDataframe(., model_grads$param_infls))

rerun_df <- 
    tibble(list=reruns_dfs) %>%
    mutate(target_param_name=names(list)) %>%
    unnest_longer(col=list, indices_to="signal") %>%
    unnest(list)

rerun_df

signal_dfs <-
    map_depth(signals, 2, ~ data.frame(
        description=.$description, n_drop=.$apip$n, prop_drop=.$apip$prop,
        target_qoi=.$qoi$name))
signal_df <-
    tibble(list=signal_dfs) %>%
    mutate(target_param_name=names(list)) %>%
    unnest_longer(col=list, indices_to="signal") %>%
    unnest(list)

signal_df


inner_join(rerun_df, signal_df, by=c("target_param_name", "signal"))

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










# Summarize the values of each QOI for each parameter for a given model_fit.
GetModelFitInferenceDataframe <- function(model_fit, param_infls) {
    stopifnot(class(model_fit) == "ModelFit")
    stopifnot(all(names(param_infls) %in% model_fit$parameter_names))
    
    GetParameterInferenceDataframe <- function(model_fit, target_index, sig_num_ses) {
        GetInferenceQOIs(beta=model_fit$betahat[target_index],
                         se=model_fit$se[target_index],
                         sig_num_ses=sig_num_ses) %>%
            purrr::imap_dfr(~ data.frame(metric=.y, value=.x))
    }
    
    summary_df <- data.frame()
    AppendRow <- function(row) bind_rows(summary_df, row)
    for (param_name in names(param_infls)) {
        param_infl <- param_infls[[param_name]]
        # We checked above that each parameter name is found.
        target_index <- which(model_fit$parameter_names == param_name)
        summary_df <-
            GetParameterInferenceDataframe(
                model_fit=model_fit,
                target_index=target_index,
                sig_num_ses=param_infl$sig_num_ses) %>%
            mutate(param_name=param_name) %>%
            AppendRow()
    }
    return(summary_df)
}



ConstructDf <- function(model_fit, signal) {
    GetModelFitInferenceDataframe(model_fit, param_infls=model_grads$param_infls) %>%
        mutate(fit="refit",
               n_drop=signal$apip$n,
               prop_drop=signal$apip$prop,
               signal_description=signal$description,
               target_metric=signal$qoi$name)
}


GetSignalsAndRefitsDataframe <- function(reruns, signals, ConstructDf) {
    stopifnot(setequal(names(reruns), names(signals)))
    summary_df <- data.frame()
    AppendRow <- function(row) bind_rows(summary_df, row)
    for (target_param_name in names(reruns)) {
        param_reruns  <- reruns[[target_param_name]]
        param_signals  <- signals[[target_param_name]]
        stopifnot(names(param_reruns) == names(param_signals))
        for (signal_name in names(param_reruns)) {
            rerun <- param_reruns[[signal_name]]
            signal <- param_signals[[signal_name]]
            stopifnot(class(rerun) == "ModelFit")
            stopifnot(class(signal) == "QOISignal")
            summary_df <-
                ConstructDf(rerun, signal) %>%
                mutate(target_param_name=target_param_name,
                       signal_name=signal_name) %>%
                AppendRow()
        }
    }
    return(summary_df)    
}


####################
# More compact

TraverseParamSignalList <- function(signals, SignalFun) {
    for (param_name in names(signals)) {
        param_signal <- signals[[param_name]]
        for (signal_name in names(param_signal)) {
            signal <- param_signal[[signal_name]]
            SignalFun(param_name, signal_name, signal)
        }
    }
}

summary_df <- data.frame()
SignalFun <- function(param_name, signal_name, signal) {
    stopifnot(param_name %in% names(reruns))
    rerun  <- reruns[[param_name]][[signal_name]]
    stopifnot(class(rerun) == "ModelFit")
    stopifnot(class(signal) == "QOISignal")
    summary_df <<- bind_rows(
        summary_df,
        ConstructDf(rerun, signal)
    )
}


summary_df <- data.frame()
TraverseParamSignalList(signals, SignalFun)
refit_df <-
    summary_df %>% mutate(fit="refit")
base_df <-
    GetModelFitInferenceDataframe(model_grads$model_fit, model_grads$param_infls) %>%
    mutate(fit="initial")






# What is the right way to do this?
library(tidyr)
summary_df %>%
    filter(metric == "beta", param_name == "x1") %>%
    select(fit, param_name, target_param_name, signal_name, value)

#############################################

# %>%
#     RerunForTargetChanges(model_grads)

# Summaries comparing reruns and predictions for each signal.
RerunSummaryDf(signals)
PlotSignal(model_grads$param_infls[["x1"]], signals[["both"]], apip_max=0.03)

# Visualize which points are being dropped

df$drop <- FALSE
df$drop[signals[["both"]]$apip$inds] <- TRUE
df$infl <- model_grads$param_infls[["x1"]][[signals[["both"]]$qoi_name]]$infl

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
fit_object <- ivreg(data = df, formula = iv_form, x=TRUE, y=TRUE)

# Get influence and reruns.
model_grads <-
    ComputeModelInfluence(fit_object) %>%
    AppendTargetRegressorInfluence("x1")
signals <-
    GetInferenceSignalsForParameter(model_grads$param_infls[["x1"]]) %>%
    RerunForTargetChanges(model_grads)

# Summaries comparing reruns and predictions for each signal.
RerunSummaryDf(signals)
PlotSignal(model_grads$param_infls[["x1"]], signals[["both"]], apip_max=0.03)



testthat::expect_equivalent(
    model_grads$model_fit$parameter_names, fit_object$x %>% colnames(),
    info="column names")

fit_object$x$regressors


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
fit_object <- lm(data=df, formula=reg_form, x=TRUE, y=TRUE)

# Get influence and reruns.
model_grads <-
    ComputeModelInfluence(fit_object) %>%
    AppendTargetRegressorInfluence("x1")

signals <-
    GetInferenceSignalsForParameter(model_grads$param_infls[["x1"]]) %>%
    RerunForTargetChanges(model_grads)

# Summaries comparing reruns and predictions for each signal.
RerunSummaryDf(signals)
grid.arrange(
    PlotSignal(model_grads$param_infls[["x1"]], signals[["sign"]], apip_max=0.03),
    PlotSignal(model_grads$param_infls[["x1"]], signals[["sig"]], apip_max=0.03),
    PlotSignal(model_grads$param_infls[["x1"]], signals[["both"]], apip_max=0.03),
    ncol=3
)


#############################
# Pairing

# In the current version of the refactor, pairing is not yet supported.
