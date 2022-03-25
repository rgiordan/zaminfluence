###################################################################
#
# These simple examples illustrate the use of zaminfluence.
# https://github.com/rgiordan/zaminfluence
# See the README.md file for installation instructions.

# Consider doing autodiff in C++
# http://www.autodiff.org/?module=Tools&language=C%2FC%2B%2B
# Adept?

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
        pivot_wider(-method, names_from=method, values_from=value)
    return(summary_df)
}


#############################
# Oridinary regression with a large number of fixed effects.

# Generate data.
set.seed(42)
x_dim <- 3
param_true <- 0.1 * runif(x_dim)
df <- GenerateRegressionData(num_obs, param_true, num_groups=100) %>%
    rename(z=se_group) %>%
    mutate(z=factor(z))
table(df$z)

# Fit a regression model.
x_names <- sprintf("x%d", 1:x_dim)
reg_form <- formula(sprintf("y ~ %s + z - 1", paste(x_names, collapse=" + ")))
fit_object <- lm(data = df, formula=reg_form, x=TRUE, y=TRUE)

# Get influence and reruns.
model_grads <-
    ComputeModelInfluence(fit_object, se_group=df$z, keep_pars="x1") %>%
    AppendTargetRegressorInfluence("x1")

if (FALSE) {
    # This is really slow because it computes grads for all the fixed effects
    model_grads <-
        ComputeModelInfluence(fit_object, se_group=df$z) %>%
        AppendTargetRegressorInfluence("x1")
}



signals <- GetInferenceSignals(model_grads)
reruns <- RerunForSignals(signals, model_grads)
preds <- PredictForSignals(signals, model_grads)
base_df <- GetModelFitInferenceDataframe(model_grads$model_fit, model_grads$param_infls)

summary_df <- SummarizeReruns(reruns, preds)

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



