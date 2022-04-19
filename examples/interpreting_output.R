###################################################################
#
# These simple examples illustrate the use of zaminfluence.
# https://github.com/rgiordan/zaminfluence
# See the README.md file for installation instructions.

library(tidyverse)
library(gridExtra)
library(zaminfluence)
#library(AER)

set.seed(42)


#############################
# Oridinary regression.

# Generate data.
num_obs <- 10000

set.seed(42)
x_dim <- 2
param_true <- c(0.1, 100)
df <- GenerateRegressionData(num_obs, param_true, num_groups=NULL)

# Fit a regression model.  We'll investigate the coefficient in front of x1.
fit_object <- lm(data=df, formula="y ~ x1 + x2 + 1", x=TRUE, y=TRUE)

# First compute the influence function itself and store the output in model_grads.
# There's no need for the user to interact with model_grads directly --- it
# just stores the relatively expensive influence function computations.
model_grads <-
    ComputeModelInfluence(fit_object, keep_pars=c("x1", "x2")) %>%
    AppendTargetRegressorInfluence("x1") %>%
    AppendTargetRegressorInfluence("x2")

# The signals object tries to identify points to change certain
# quantities of interest: the sign of an estimator, its significance,
# or to produce a significant change of the opposite sign.
signals <- GetInferenceSignals(model_grads)


###################################
# Interpreting signals

# Signals is a list with one entry for every parameter for which you
# ran AppendTargetRegressorInfluence.
signals %>% names()

# Each parameter has an analysis for each quantity of interest (QOI):
# sign, significance ("sig"), or a significant result of the opposite sign ("both")
signals[["x1"]] %>% names()

signal <- signals[["x1"]][["sign"]] 

# For a particular QOI, the "apip" tells you the "approximate perturbation-inducing proportion" ---
# that is, the proportion of points you need to drop to change the QOI.
signal$apip

# If success is TRUE, then the linear approximation is capable of finding
# points that reverse the QOI.  (See below for an example where sucess is FALSE)
stopifnot(signal$apip$success)

sprintf("Proportion of points needed to change x1 sign: %f",
        signal$apip$prop)

sprintf("Number of points needed to change x1 sign: %d",
        signal$apip$n)

# inds gives which points you need to drop from the original dataset.
# You can use this, for exmample, to visualize which points are being dropped.

# You can do this manually like so:
df$drop <- FALSE
df$drop[signal$apip$inds] <- TRUE

# Or with the GetWeightVector helper:
df$drop_v2 <- GetWeightVector(signal$apip$inds, nrow(df), bool=TRUE, invert=TRUE)
stopifnot(all(df$drop_v2 == df$drop))

# Show that points that are dropped have large residuals and large |x1|
df$resid <- fit_object$residuals
ggplot(df) + geom_point(aes(x=x1, y=resid, color=drop))


# Because the x2 coefficient was so large, the linear approximation
# cannot find a way to change it.  So "success" is FALSE:
signals[["x2"]][["sign"]]$apip$success

# ... and all the information is NA.
signals[["x2"]][["sign"]]$apip$prop
signals[["x2"]][["sign"]]$apip$n
signals[["x2"]][["sign"]]$apip$inds

# Just because the linear approximation cannot produce a change doesn't
# mean it can't be produced, though it is evidence that it cannot be produced
# by dropping a /small/ proportion of the data.  See Section 3.3.2 of 
# the paper https://arxiv.org/pdf/2011.14999.pdf for more discussion.




###################################
# Reruns and predictions

# You can (and should!) check the linear approximation by rerunning
# without the data points selected for the APIP.  There are convenient
# wrapper functions that do this for every signal:

reruns <- RerunForSignals(signals, model_grads)

# PredictForSignals doesn't re-run, but it uses the linear approximation to
# predict the changes, and returns a list in the same format as RerunForSignals.
# This can be convenient for comparing the reruns and predictions side-by-side.
preds <- PredictForSignals(signals, model_grads)

# The structure of the output is a list of lists, just like signals:
reruns %>% names()
reruns[["x1"]] %>% names()

preds %>% names()
preds[["x1"]] %>% names()

# For easy summarization, you can convert the lists into dataframes.
base_df <- GetModelFitInferenceDataframe(model_grads$model_fit, model_grads$param_infls)
reruns_df <- GetSignalsAndRerunsDataframe(signals, reruns, model_grads)
preds_df <- GetSignalsAndRerunsDataframe(signals, preds, model_grads)

# For example: when we drop n_drop points targeting a sign change in x1 and re-run
# the regression, we can see how the x2 coefficient ("param"), standard error ("se"),
# and confidence interval ("param_mzse", "param_pzse") change:
reruns_df %>%
    filter(target_param_name == "x1", target_signal == "sign",
           param_name == "x2") %>%
    select(param_name, metric, value, n_drop)

# It's easy then to summarize things however you like.
summary_df <-
    bind_rows(reruns_df %>% mutate(method="rerun"),
              preds_df %>% mutate(method="prediction")) %>%
    pivot_wider(names_from=method, values_from=value) %>%
    inner_join(base_df %>% rename(original=value),
               by=c("metric", "param_name"))

# For example, you can graph the reruns and predictions versus one another like so:
ggplot(summary_df) +
    geom_point(aes(x=prediction, y=rerun, color=param_name, shape=metric)) +
    geom_abline(aes(slope=1, intercept=0))

# View the effects on the target parameter of reruns and predictions for different signals.
summary_df %>%
    filter(param_name == target_param_name, metric == "param") %>%
    select(param_name, description, n_drop, prop_drop, original, rerun, prediction)


# An earlier version of the paper visualized signals like so, and the
# code is still there.  The plot shows the cumulative change in certain
# quantities of interest as more and more points are dropped, as well as the
# actual effect of dropping points given by the APIP.
PlotSignal(model_grads, signals, "x1", "sign",
           reruns=reruns, apip_max=0.03)

# Summaries comparing reruns and predictions for each signal.
plots <-
    map(c("sign", "sig", "both"),
        ~ PlotSignal(model_grads, signals, "x1", ., reruns=reruns, apip_max=0.03))
do.call(function(...) { grid.arrange(..., ncol=1) }, plots)

