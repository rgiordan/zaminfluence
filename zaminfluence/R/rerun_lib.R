#' Rerun the model at the AMIS for a set of signals.
#' @param signals `A list of signal objects, "sign", "sig", "both",
#' as produced by [GetInferenceSignalsForParameter].
#' @param model_grads `r docs$model_grads`
#' @param RerunFun (Optional) A function taking as inputs
#' a vector of weights, and returning a ModelFit object.  If unspecified,
#' `model_grads$RerunFun` is used.
#' @param verbose (Optional) If TRUE, print status updates.
#'
#' @return The original `signals`, with two new fields:
#' - rerun_df: A dataframe summarizing the rerun.
#' - rerun: The complete rerun result.
#' @export
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

    num_obs <- model_grads$model_fit$n_obs
    RerunSignal <- function(signal) {
        verbosePrint("Rerunning ", signal$description, "\n")
        w_bool <- GetWeightVector(
            drop_inds=signal$apip$inds,
            num_obs=num_obs,
            bool=TRUE)
        return(RerunFun(w_bool))
    }
    reruns <- map_depth(signals, 2, RerunSignal)
    return(reruns)
}


#' Predict the model at the AMIS for a set of signals.
#' @param signals `A list of signal objects, "sign", "sig", "both",
#' as produced by [GetInferenceSignalsForParameter].
#' @param model_grads `r docs$model_grads`
#' @param verbose (Optional) If TRUE, print status updates.
#'
#' @return The original `signals`, with two new fields:
#' - rerun_df: A dataframe summarizing the rerun.
#' - rerun: The complete rerun result.
#' @export
PredictForSignals <- function(signals, model_grads, verbose=FALSE) {
  PredictFun <- function(weights) {
    PredictModelFit(model_grads, weights)
  }
  RerunForSignals(signals, model_grads, verbose=verbose, RerunFun=PredictFun)
}

# Here's the old format:
#       signals[[target]]$rerun_df <-
#           as.data.frame(signal) %>%
#               mutate(
#                   betahat_refit=rerun_base_values["beta"],
#                   se_refit=rerun_base_values["se"],
#                   beta_mzse_refit=rerun_base_values["beta_mzse"],
#                   beta_pzse_refit=rerun_base_values["beta_pzse"],
#                   betahat_orig=orig_base_values["beta"],
#                   se_orig=se_orig,
#                   beta_mzse_orig=orig_base_values["beta_mzse"],
#                   beta_pzse_orig=orig_base_values["beta_pzse"])
