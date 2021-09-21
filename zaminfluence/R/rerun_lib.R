#' Compute the "base values" (beta and the CI) for a parameter from a rerun.
#' @param rerun `r docs$rerun`
#' @param param_infl `r docs$param_infl`
#'
#' @return A named numeric vector of base values.
#' @export
GetRerunBaseValues <- function(rerun, param_infl) {
    target_index <- param_infl$target_index
    sig_num_ses <- param_infl$sig_num_ses
    beta <- rerun$betahat[target_index]
    se <- rerun$se[target_index]
    base_values <-
        c(beta,
          beta - sig_num_ses * se,
          beta + sig_num_ses * se)
    names(base_values) <- c("beta", "beta_mzse", "beta_pzse")
    return(base_values)
}



#' Rerun the model at the AMIS for a set of signals.
#' @param signals `A list of signal objects, "sign", "sig", "both",
#' as produced by [GetInferenceSignals].
#' @param model_grads `r docs$model_grads`
#' @param RerunFun (Optional) A function taking as inputs
#' `model_grads$model_fit`
#' and a vector of weights, and return a rerun object.  If unspecified,
#' `model_grads$RerunFun` is used.
#'
#' @return The original `signals`, with two new fields:
#' - rerun_df: A dataframe summarizing the rerun.
#' - rerun: The complete rerun result.
#' @export
RerunForTargetChanges <- function(signals, model_grads, RerunFun=NULL) {

  if (is.null(RerunFun)) {
    RerunFun <- model_grads$RerunFun
  }
  for (target in c("sign", "sig", "both")) {
      signal <- signals[[target]]
      w_bool <- GetWeightVector(
          drop_inds=signal$apip$inds,
          num_obs=model_grads$n_obs,
          bool=TRUE)

      rerun <- RerunFun(model_grads$model_fit, w_bool=w_bool)

      # Save the whole rerun
      signals[[target]]$rerun <- rerun

      param_infl <- model_grads$param_infl_list[[signals$target_parameter]]

      # Make a nice dataframe with the targeted regressor
      orig_base_values <- GetBaseValues(param_infl)
      rerun_base_values <- GetRerunBaseValues(rerun, param_infl)

      # The standard error is not a QOI, but it is convenient to
      # report it.
      se_orig <- model_grads$se[param_infl$target_index]
      se_refit <- rerun$se[param_infl$target_index]

      signals[[target]]$rerun_df <-
          GetSignalDataFrame(signal) %>%
              mutate(
                  betahat_refit=rerun_base_values["beta"],
                  se_refit=rerun_base_values["se"],
                  beta_mzse_refit=rerun_base_values["beta_mzse"],
                  beta_pzse_refit=rerun_base_values["beta_pzse"],
                  betahat_orig=orig_base_values["beta"],
                  se_orig=se_orig,
                  beta_mzse_orig=orig_base_values["beta_mzse"],
                  beta_pzse_orig=orig_base_values["beta_pzse"])
  }
  return(signals)
}
