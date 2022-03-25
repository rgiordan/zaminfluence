###################################################################
#
# These simple examples illustrate the use of zaminfluence.
# https://github.com/rgiordan/zaminfluence
# See the README.md file for installation instructions.

library(tidyverse)
library(gridExtra)
library(zaminfluence)
library(AER)

compare <- function(x, y) { return(max(abs(x - y))) }
check_equivalent  <- function(x, y) { stopifnot(compare(x, y) < 1e-8) }

num_obs <- 10000

set.seed(42)

SummarizeReruns <- function(reruns, preds) {
    reruns_df <- GetSignalsAndRerunsDataframe(signals, reruns, model_grads)
    preds_df <- GetSignalsAndRerunsDataframe(signals, preds, model_grads)
    
    summary_df <-
        rbind(reruns_df %>% mutate(method="rerun"),
              preds_df %>% mutate(method="prediction")) %>%
        pivot_wider(names_from=method, values_from=value)
    return(summary_df)
}


#############################
# Oridinary regression.

# Generate data.
set.seed(42)
x_dim <- 3
param_true <- 0.1 * runif(x_dim)
df <- GenerateRegressionData(num_obs, param_true, num_groups=NULL)

# Fit a regression model.
x_names <- sprintf("x%d", 1:x_dim)
reg_form <- formula(sprintf("y ~ %s - 1", paste(x_names, collapse=" + ")))
fit_object <- lm(data = df, formula=reg_form, x=TRUE, y=TRUE)

# Get influence and reruns.
model_grads <-
    ComputeModelInfluence(fit_object, keep_pars=c("x2", "x1")) %>%
    AppendTargetRegressorInfluence("x1") %>%
    AppendTargetRegressorInfluence("x2")

signals <- GetInferenceSignals(model_grads)
reruns <- RerunForSignals(signals, model_grads)



preds <- PredictForSignals(signals, model_grads)
base_df <- GetModelFitInferenceDataframe(model_grads$model_fit, model_grads$param_infls)

summary_df <- SummarizeReruns(reruns, preds)
View(summary_df)

ggplot(summary_df) +
    geom_point(aes(x=prediction, y=rerun, color=param_name, shape=metric)) +
    geom_abline(aes(slope=1, intercept=0))

PlotSignal(model_grads, signals, "x1", "sign",
          reruns=reruns, apip_max=0.03)


# Visualize which points are being dropped
signal <- signals[["x1"]][["both"]] 
df$drop <- GetWeightVector(signal$apip$inds, nrow(df), bool=TRUE, invert=TRUE)
df$infl <- signal$qoi$infl

grid.arrange(
    ggplot(df) +
        geom_point(aes(x=x1, y=infl, color=drop)),
    ggplot(df) +
        geom_point(aes(x=x1, y=y, color=drop)),
    ggplot(df) +
        geom_point(aes(x=x1, y=x2, color=drop)),
    ncol=3
)


# Can we structure the signals better?  Yes.

param_name <- "x1"
this_signals <- GetInferenceSignalsForParameter(model_grads$param_infls[[param_name]])
names(this_signals)

signals_df <- tibble()
for (qoi_name in names(this_signals)) {
    signals_df <- bind_rows(
        signals_df,
        tibble(param=param_name, qoi_name=qoi_name) %>%
            mutate(signal=list(this_signals[[qoi_name]])))
}

this_signals[["sign"]] %>% pluck("apip") %>% pluck("prop")
this_signals[["sig"]] %>% pluck("apip") %>% pluck("prop")
this_signals[["both"]] %>% pluck("apip") %>% pluck("prop")

signals_df %>% rowwise() %>% mutate(apip=pluck(signal, "apip") %>% pluck("prop"))
signals_df$signal[[1]]$apip$prop
signals_df$signal[[2]]$apip$prop



foo <- list(a=1)
bar <- list(a=2, c=3)
baz <- list(a=3.5, c=2,8, d=6)


tb <- tibble(ind=1:3) %>%
    mutate(l=list(foo, bar, baz))


tb %>% mutate(len=length(l))
tb %>%  rowwise() %>% mutate(len=length(l))
tb %>%  rowwise() %>% mutate(myval=pluck(l, "a"))


#############################
# Instrumental variables.

# Generate data.
set.seed(42)
x_dim <- 3
param_true <- 0.1 * runif(x_dim)
df <- GenerateIVRegressionData(num_obs, param_true, num_groups=NULL)

# Fit an IV model.
x_names <- sprintf("x%d", 1:x_dim)
z_names <- sprintf("z%d", 1:x_dim)
iv_form <- formula(sprintf("y ~ %s - 1 | %s - 1",
                           paste(x_names, collapse=" + "),
                           paste(z_names, collapse=" + ")))
fit_object <- ivreg(data = df, formula = iv_form, x=TRUE, y=TRUE)

# Get influence and reruns.
model_grads <-
    ComputeModelInfluence(fit_object, keep_pars=c("x2", "x1")) %>%
    AppendTargetRegressorInfluence("x1")

signals <- GetInferenceSignals(model_grads)
reruns <- RerunForSignals(signals, model_grads)
preds <- PredictForSignals(signals, model_grads)

summary_df <- SummarizeReruns(reruns, preds)
ggplot(summary_df) +
    geom_point(aes(x=prediction, y=rerun, color=param_name, shape=metric)) +
    geom_abline(aes(slope=1, intercept=0))

PlotSignal(model_grads, signals, "x1", "sign",
           reruns=reruns, apip_max=0.03)


#############################
# Grouped standard errors.

# Generate data.
set.seed(42)
x_dim <- 3
param_true <- 0.1 * runif(x_dim)
num_groups <- 50
df <- GenerateRegressionData(num_obs, param_true, num_groups=num_groups)

# se_group is zero-indexed group indicator with no missing entries.
table(df$se_group)

# Fit a regression model.
x_names <- sprintf("x%d", 1:x_dim)
reg_form <- formula(sprintf("y ~ %s - 1", paste(x_names, collapse=" + ")))
fit_object <- lm(data=df, formula=reg_form, x=TRUE, y=TRUE)

# Get influence and reruns.  Pass the grouping indicator to the `se_group`` argument
# of `ComputeModelInfluence`.
model_grads <-
    ComputeModelInfluence(fit_object, se_group=df$se_group) %>%
    AppendTargetRegressorInfluence("x1")


# The grouped standard error which zaminfluence computes...
cat("Zaminfluence SE:\t", model_grads$param_infls[["x1"]]$se$base_value, "\n")

# ...is equivalent to the that computed by the following standard command:
cat("vcovCL se:\t\t", 
    vcovCL(fit_object, cluster=df$se_group, type="HC0", cadjust=FALSE)["x1", "x1"] %>% sqrt(), 
    "\n")

signals <- GetInferenceSignals(model_grads)
reruns <- RerunForSignals(signals, model_grads)
preds <- PredictForSignals(signals, model_grads)
summary_df <- SummarizeReruns(reruns, preds)

ggplot(summary_df) +
    geom_point(aes(x=prediction, y=rerun, color=param_name, shape=metric)) +
    geom_abline(aes(slope=1, intercept=0))

# Summaries comparing reruns and predictions for each signal.
plots <-
    map(c("sign", "sig", "both"),
        ~ PlotSignal(model_grads, signals, "x1", ., reruns=reruns, apip_max=0.03))
do.call(function(...) { grid.arrange(..., ncol=1) }, plots)


#############################
# Pairing

# In the current version of zaminfluence, pairing is not yet supported.  If you're particularly
# interested in paried analysis, please reach out to the package authors.
