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
        weights <- GetWeightVector(
            drop_inds=signal$apip$inds,
            orig_weights=model_grads$model_fit$weights)
        return(RerunFun(weights))
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
#                   param_refit=rerun_base_values["param"],
#                   se_refit=rerun_base_values["se"],
#                   param_mzse_refit=rerun_base_values["param_mzse"],
#                   param_pzse_refit=rerun_base_values["param_pzse"],
#                   param_orig=orig_base_values["param"],
#                   se_orig=se_orig,
#                   param_mzse_orig=orig_base_values["param_mzse"],
#                   param_pzse_orig=orig_base_values["param_pzse"])
