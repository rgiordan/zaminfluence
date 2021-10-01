library(purrr)


# Define an ParameterInferenceInfluence S3 object.
new_ParameterInferenceInfluence <- function(
    target_index, target_parameter, sig_num_ses,
    se_qoi, param_qoi, param_mzse_qoi, param_pzse_qoi) {
  return(structure(
    list(
        target_index=as.integer(target_index),
        target_parameter=as.character(target_parameter),
        se=se_qoi,
        param=param_qoi,
        param_mzse=param_mzse_qoi,
        param_pzse=param_pzse_qoi,
        sig_num_ses=sig_num_ses,
        qoi_names=c("param", "param_mzse", "param_pzse", "se")),
    class="ParameterInferenceInfluence"))
}



validate_ParameterInferenceInfluence <- function(param_infl) {
  stopifnot(class(param_infl) == "ParameterInferenceInfluence")
  stopifnot(setequal(
    param_infl$qoi_names,
    c("param", "param_mzse", "param_pzse", "se")))

  for (qoi_name in param_infl$qoi_names) {
    stopifnot(class(param_infl[[qoi_name]]) == "QOIInfluence")
  }

  sig_num_ses <- param_infl$sig_num_ses
  stopifnot(is.numeric(sig_num_ses))
  stopifnot(length(sig_num_ses) == 1)
  stopifnot(sig_num_ses >= 0)

  target_index <- param_infl$target_index
  stopifnot(is.numeric(target_index))
  stopifnot(length(target_index) == 1)
  stopifnot(target_index > 0)

  return(invisible(param_infl))
}


#'@export
GetParameterInferenceQOIs <- function(model_fit, target_parameter,
                                      sig_num_ses=qnorm(0.975)) {
  stopifnot(class(model_fit) == "ModelFit")
  target_index <- which(model_fit$parameter_names == target_parameter)
  if (length(target_index) != 1) {
      stop(paste0("Target parameter ",
         target_parameter, " not found in the parameter names (",
         paste(model_fit$parameter_names, collapse=", ")
       ), ")\n")
  }

  param <- model_fit$param[target_index]
  se <- model_fit$se[target_index]

  values <- GetInferenceQOIs(param=param, se=se, sig_num_ses=sig_num_ses)
  values$target_index <- target_index
  return(values)
}


GetInferenceQOIs <- function(param, se, sig_num_ses) {
  # Remove names so we can use unlist() and get expected names.
  param <- unname(param)
  se <- unname(se)
  sig_num_ses <- unname(sig_num_ses)
  return(list(
    param=param,
    se=se,
    param_mzse=param - sig_num_ses * se,
    param_pzse=param + sig_num_ses * se
  ))
}


ParameterInferenceInfluence <- function(model_grads, target_parameter,
                                        sig_num_ses=qnorm(0.975)) {
    stopifnot(class(model_grads) == "ModelGrads")

    weights <- model_grads$model_fit$weights

    qoi_base_values <- GetParameterInferenceQOIs(
      model_grads$model_fit, target_parameter=target_parameter,
      sig_num_ses=sig_num_ses)
    target_index <- qoi_base_values$target_index

    # se_grad and param_grad are the gradients for a parameter along the path
    # taking a weight
    # from its current value to zero.  That way the "gradient" measures the
    # effect of removing a datapoint (taking its value from the current weight
    # to zero.) So we multiply the raw weight
    # derivatives by the actual base weights.
    se_grad <- weights * model_grads$se_grad[target_index,]
    param_grad <- weights * model_grads$param_grad[target_index, ]
    n_obs <- model_grads$model_fit$n_obs
    qoi_gradients <- GetInferenceQOIs(
      param=param_grad,
      se=se_grad,
      sig_num_ses=sig_num_ses)

    param_infl <- new_ParameterInferenceInfluence(
          target_index=target_index,
          target_parameter=target_parameter,
          sig_num_ses=sig_num_ses,
          se_qoi=QOIInfluence(
              name="se",
              infl=qoi_gradients$se,
              base_value=qoi_base_values$se,
              num_obs=n_obs),
          param_qoi=QOIInfluence(
              name="param",
              infl=qoi_gradients$param,
              base_value=qoi_base_values$param,
              num_obs=n_obs),
          param_mzse_qoi=QOIInfluence(
              name="param_mzse",
              infl=qoi_gradients$param_mzse,
              base_value=qoi_base_values$param_mzse,
              num_obs=n_obs),
          param_pzse_qoi=QOIInfluence(
              name="param_pzse",
              infl=qoi_gradients$param_pzse,
              base_value=qoi_base_values$param_pzse,
              num_obs=n_obs))

    validate_ParameterInferenceInfluence(param_infl)

    return(param_infl)
}


#' Compute the influence scores for a particular parameter.
#' @param model_grads `r docs$model_grads`
#' @param target_parameter The string naming a regressor (must be an
#' entry in [model_grads$parameter_names]).
#' @param sig_num_ses `r docs$sig_num_ses`
#'
#' @return The original `model_grads`, with an entry
#' `model_grads$param_infls[[target_parameter]]` containing a
#' parameter influence object.
#'
#'@export
AppendTargetRegressorInfluence <- function(model_grads, target_parameter,
                                           sig_num_ses=qnorm(0.975)) {
    stopifnot(class(model_grads) == "ModelGrads")
    if (is.null(model_grads[["param_infls"]])) {
        model_grads$param_infls <- list()
    }

    param_infl <- ParameterInferenceInfluence(
      model_grads, target_parameter, sig_num_ses=sig_num_ses)

    model_grads$param_infls[[target_parameter]] <- param_infl
    return(invisible(model_grads))
}


#' @export
GetBaseValues <- function(param_infl) {
  stopifnot(class(param_infl) == "ParameterInferenceInfluence")
  return(map_dbl(param_infl[param_infl$qoi_names], ~ .$base_value))
}


# Define an QOISignal S3 object.
new_QOISignal <- function(qoi, signal, description, apip) {
  return(structure(
    list(qoi=qoi, signal=signal, description=description, apip=apip),
    class="QOISignal"
  ))
}


validate_QOISignal <- function(signal) {
  stopifnot(class(signal) == "QOISignal")
  stopifnot(class(signal$qoi) == "QOIInfluence")
  stopifnot(class(signal$apip) == "APIP")
  StopIfNotNumericScalar(signal$signal)
  return(invisible(signal))
}


QOISignal <- function(qoi, signal, description) {
  # Note that the qoi data is not copied.
  return(validate_QOISignal(new_QOISignal(
    qoi=qoi,
    signal=signal,
    description=description,
    apip=GetAPIPForQOI(qoi=qoi, signal=signal)
  )))
}


#' Compute the signals for changes to sign, significance, and both.
#' @param model_grads `r docs$model_grads`
#'
#' @return A list lists of of signals, one for each parameter in
#' named model_grads$param_infls.  See the output of
#' `GetInferenceSignalsForParameter`.
#'
#' @export
GetInferenceSignals <- function(model_grads) {
  stopifnot(class(model_grads) == "ModelGrads")
  signals <- list()
  for (param_name in names(model_grads$param_infls)) {
    signals[[param_name]] <- GetInferenceSignalsForParameter(
      model_grads$param_infls[[param_name]]
    )
  }
  return(signals)
}


#' Compute the signals for changes to sign, significance, and both.
#' @param param_infl `r docs$param_infl`
#'
#' @return A list of signals, named "sign", "sig", and "both".  Each
#' entry is a `QOISignal` object.
#' @export
GetInferenceSignalsForParameter <- function(param_infl) {
    stopifnot(class(param_infl) == "ParameterInferenceInfluence")
    base_values <- GetBaseValues(param_infl)
    param <- base_values["param"]
    param_mzse <- base_values["param_mzse"]
    param_pzse <- base_values["param_pzse"]

    sign_label <- "sign"
    sig_label <- "significance"
    both_label <- "sign and significance"

    signals <- list()
    #signals$target_parameter <- param_infl$target_parameter
    signals$sign <- QOISignal(
      qoi=param_infl[["param"]],
      signal=-1 * param,
      description=sign_label)

    is_significant <- sign(param_mzse) == sign(param_pzse)
    if (is_significant) {
        if (param_mzse >= 0) { # then param_pzse > 0 too because significant
            signals$sig <- QOISignal(
              qoi=param_infl[["param_mzse"]],
              signal=-1 * param_mzse,
              description=sig_label)
            signals$both  <- QOISignal(
              qoi=param_infl[["param_pzse"]],
              signal=-1 * param_pzse,
              description=both_label)
        } else if (param_pzse < 0) { # then param_mzse < 0 too because significant
            signals$sig <- QOISignal(
              qoi=param_infl[["param_pzse"]],
              signal=-1 * param_pzse,
              description=sig_label)
            signals$both <- QOISignal(
                qoi=param_infl[["param_mzse"]],
                signal=-1 * param_mzse,
                description=both_label)
        } else {
            stop("Impossible for a significant result")
        }
    } else { # Not significant.  Choose to change the interval endpoint which
             # is closer.
        if (abs(param_mzse) >= abs(param_pzse)) {
            signals$sig <- QOISignal(
              qoi=param_infl[["param_pzse"]],
              signal=-1 * param_pzse,
              description=sig_label)
        } else  {
            signals$sig <- QOISignal(
              qoi=param_infl[["param_mzse"]],
              signal=-1 * param_mzse,
              description=sig_label)
        }

        if (param >= 0) {
            signals$both <- QOISignal(
                qoi=param_infl[["param_mzse"]],
                signal=-1 * param_mzse,
                description=both_label)
        } else {
            signals$both <- QOISignal(
                qoi=param_infl[["param_mzse"]],
                signal=-1 * param_mzse,
                description=both_label)
        }
    }

    return(signals)
}


#' Produce a tidy dataframe summarizing a signal.
#' @param signal `r docs$signal`
#'
#'@export
as.data.frame.QOISignal <- function(signal) {
  data.frame(
      qoi_name=signal$qoi$name,
      description=signal$description,
      signal=signal$signal,
      num_removed=signal$apip$n,
      prop_removed=signal$apip$prop)
}
