library(ggplot2)
library(purrr)
library(latex2exp)


# Define an ParameterInferenceInfluence S3 object.

new_ParameterInferenceInfluence <- function(
    target_index, target_parameter, sig_num_ses,
    se_qoi, beta_qoi, beta_mzse_qoi, beta_pzse_qoi) {
  return(structure(
    list(
        target_index=as.integer(target_index),
        target_parameter=as.character(target_parameter),
        se=se_qoi,
        beta=beta_qoi,
        beta_mzse=beta_mzse_qoi,
        beta_pzse=beta_pzse_qoi,
        sig_num_ses=sig_num_ses,
        qoi_names=c("beta", "beta_mzse", "beta_pzse", "se")),
    class="ParameterInferenceInfluence"))
}



validate_ParameterInferenceInfluence <- function(param_infl) {
  stopifnot(class(param_infl) == "ParameterInferenceInfluence")
  stopifnot(setequal(
    param_infl$qoi_names,
    c("beta", "beta_mzse", "beta_pzse", "se")))

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


ParameterInferenceInfluence <- function(model_grads, target_parameter,
                                        sig_num_ses=qnorm(0.975)) {
    stopifnot(class(model_grads) == "ModelGrads")
    target_index <- which(
      model_grads$model_fit$parameter_names == target_parameter)
    if (length(target_index) != 1) {
        stop("Error finding target regressor in the regression.")
    }

    weights <- model_grads$model_fit$weights
    se_grad <- weights * model_grads$se_grad[target_index,]
    beta_grad <- weights * model_grads$beta_grad[target_index, ]
    betahat <- model_grads$model_fit$betahat[target_index]
    sehat <- model_grads$model_fit$se[target_index]
    n_obs <- model_grads$model_fit$n_obs

    param_infl <- new_ParameterInferenceInfluence(
          target_index=target_index,
          target_parameter=target_parameter,
          sig_num_ses=sig_num_ses,
          se_qoi=QOIInfluence(
              name="se",
              infl=se_grad,
              base_value=sehat,
              num_obs=n_obs),
          beta_qoi=QOIInfluence(
              name="beta",
              infl=beta_grad,
              base_value=betahat,
              num_obs=n_obs),
          beta_mzse_qoi=QOIInfluence(
              name="beta_mzse",
              infl=beta_grad - sig_num_ses * se_grad,
              base_value=betahat - sig_num_ses * sehat,
              num_obs=n_obs),
          beta_pzse_qoi=QOIInfluence(
              name="beta_pzse",
              infl=beta_grad + sig_num_ses * se_grad,
              base_value=betahat + sig_num_ses * sehat,
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
#' `model_grads$param_infl_list[[target_parameter]]` containing a
#' parameter influence object.
#'
#'@export
AppendTargetRegressorInfluence <- function(model_grads, target_parameter,
                                           sig_num_ses=qnorm(0.975)) {
    stopifnot(class(model_grads) == "ModelGrads")
    if (is.null(model_grads[["param_infl_list"]])) {
        model_grads$param_infl_list <- list()
    }

    param_infl <- ParameterInferenceInfluence(
      model_grads, target_parameter, sig_num_ses=sig_num_ses)

    model_grads$param_infl_list[[target_parameter]] <- param_infl
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
  return(validate_QOISignal(new_QOISignal(
    qoi=qoi,
    signal=signal,
    description=description,
    apip=GetAPIPForQOI(qoi=qoi, signal=signal)
  )))
}


#' Compute the signals for changes to sign, significance, and both.
#' @param param_infl `r docs$param_infl`
#'
#' @return A list of signals, named "sign", "sig", and "both".  Each
#' entry is a `signal` object.
#' @export
GetInferenceSignals <- function(param_infl) {
    stopifnot(class(param_infl) == "ParameterInferenceInfluence")
    base_values <- GetBaseValues(param_infl)
    betahat <- base_values["beta"]
    beta_mzse <- base_values["beta_mzse"]
    beta_pzse <- base_values["beta_pzse"]

    sign_label <- "sign"
    sig_label <- "significance"
    both_label <- "sign and significance"

    signals <- list()
    signals$target_parameter <- param_infl$target_parameter
    signals$sign <- QOISignal(
      qoi=param_infl[["beta"]],
      signal=-1 * betahat,
      description=sign_label)

    is_significant <- sign(beta_mzse) == sign(beta_pzse)
    if (is_significant) {
        if (beta_mzse >= 0) { # then beta_pzse > 0 too because significant
            signals$sig <- QOISignal(
              qoi=param_infl[["beta_mzse"]],
              signal=-1 * beta_mzse,
              description=sig_label)
            signals$both  <- QOISignal(
              qoi=param_infl[["beta_pzse"]],
              signal=-1 * beta_pzse,
              description=both_label)
        } else if (beta_pzse < 0) { # then beta_mzse < 0 too because significant
            signals$sig <- QOISignal(
              qoi=param_infl[["beta_pzse"]],
              signal=-1 * beta_pzse,
              description=sig_label)
            signals$both <- QOISignal(
                qoi=param_infl[["beta_mzse"]],
                signal=-1 * beta_mzse,
                description=both_label)
        } else {
            stop("Impossible for a significant result")
        }
    } else { # Not significant.  Choose to change the interval endpoint which
             # is closer.
        if (abs(beta_mzse) >= abs(beta_pzse)) {
            signals$sig <- QOISignal(
              qoi=param_infl[["beta_pzse"]],
              signal=-1 * beta_pzse,
              description=sig_label)
        } else  {
            signals$sig <- QOISignal(
              qoi=param_infl[["beta_mzse"]],
              signal=-1 * beta_mzse,
              description=sig_label)
        }

        if (betahat >= 0) {
            signals$both <- QOISignal(
                qoi=param_infl[["beta_mzse"]],
                signal=-1 * beta_mzse,
                description=both_label)
        } else {
            signals$both <- QOISignal(
                qoi=param_infl[["beta_mzse"]],
                signal=-1 * beta_mzse,
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


################################################################################
# Plotting and visualization functions


#' Produce an influence dataframe suitable for visualization.
#' @param param_infl `r docs$param_infl`
#' @param sorting_qoi_name The name of a QOI ("beta", "beta_mzse", or
#' "beta_pzse") whose influence is used to sort the dataframe.
#' @param max_num_obs (Optional)  Include at most the `max_num_obs`
#' most influential observations for the sorting QOI.
#' Default is to include all observations.
#'
#' @return A dataframe with predictions, leaving out cumulatively more
#' points according to the sorting QOI's influence scores.
#' @export
GetSortedInfluenceDf <- function(param_infl, sorting_qoi_name,
                                 max_num_obs=Inf) {
    stopifnot(class(param_infl) == "ParameterInferenceInfluence")
    stopifnot(sorting_qoi_name %in% param_infl$qoi_names)
    qoi_for_sorting <- param_infl[[sorting_qoi_name]]

    GetQOIDf <- function(infl_sign) {
        ordered_inds <- qoi_for_sorting[[infl_sign]]$infl_inds
        if (max_num_obs < length(ordered_inds)) {
            ordered_inds <- ordered_inds[1:max_num_obs]
        }
        qoi_df <- data.frame(
            num_dropped=c(0, 1:length(ordered_inds))) %>%
            mutate(prop_dropped=num_dropped /
                       qoi_for_sorting[[infl_sign]]$num_obs,
                   sign=infl_sign)
        for (qoi_name in c("beta", "beta_mzse", "beta_pzse")) {
            base_value <- param_infl[[qoi_name]]$base_value
            infl_sorted <- param_infl[[qoi_name]]$infl[ordered_inds]
            qoi_df[[qoi_name]] <- base_value + cumsum(c(0, infl_sorted))
        }
        return(qoi_df)
    }

    qoi_df <-
        bind_rows(GetQOIDf("pos"), GetQOIDf("neg")) %>%
        mutate(sorted_by=sorting_qoi_name)

    return(qoi_df)
}


#' Plot influence scores, signals, and reruns.
#' @param influence_df The output of [GetSortedInfluenceDf]
#' @param plot_num_dropped If TRUE, plot the number dropped on the x-axis.
#' If FALSE (the default), plot the proportion dropped.
#' @param apip_max The maximum value for the x-axis (as a number or proportion
#' according to the value of `plot_num_dropped`).
#' @param signals (Optional) A list of signals to plot.
#' @param include_y_zero (Optional) If TRUE (the default), force the y-axis
#' to include zero and plot a horizontal line.
#'
#' @return A plot.
#' @export
PlotInfluence <- function(influence_df,
                          plot_num_dropped=FALSE,
                          apip_max=NULL,
                          signals=NULL,
                          include_y_zero=TRUE) {

    PlotRegSignal <- function(plot, signal) {
        alpha_type <- if (plot_num_dropped) "n" else "prop"
        alpha <- signal$apip[[alpha_type]]
        if (is.null(apip_max) || (!is.null(apip_max) && alpha <= apip_max)) {
            plot <- plot + geom_vline(aes(xintercept=!!alpha,
                                          linetype=!!signal$description))
            if (!is.null(signal$rerun_df)) {
                rerun_df <- signal$rerun_df
                stopifnot(nrow(rerun_df) == 1)
                plot <-
                    plot +
                    geom_errorbar(aes(
                        x=!!alpha,
                        ymin=beta_mzse_refit,
                        ymax=beta_pzse_refit),
                        data=rerun_df,
                        width=errorbar_width,
                        lwd=1.5) +
                    geom_point(aes(x=!!alpha, y=betahat_refit),
                               data=rerun_df,
                               shape=8)
            }
        }
        return(plot)
    }

    influence_df$alpha <-
        if (plot_num_dropped) influence_df$num_dropped else
            influence_df$prop_dropped

    if (!is.null(apip_max)) {
      influence_df <- filter(influence_df, alpha  <= apip_max)
    }
    errorbar_width <- diff(range(influence_df$alpha)) / 50

    plot <- ggplot(influence_df, aes(x=alpha))
    if (include_y_zero) {
        plot <-
            plot +
            geom_line(aes(y=0.0), col="gray50")
    }
    if (!is.null(signals)) {
        for (signal in signals) {
            plot <- PlotRegSignal(plot,  signal)
        }
        plot <- plot + guides(linetype=guide_legend(title="Change type"))
    }

    base_beta <- filter(influence_df, alpha == 0) %>% pull("beta") %>% unique()
    stopifnot(length(base_beta) == 1)
    plot <-
        plot +
        geom_line(aes(y=!!base_beta), col="blue", lwd=2) +
        geom_ribbon(aes(
            ymin=beta_mzse,
            ymax=beta_pzse,
            group=sign),
            fill="blue", color=NA, alpha=0.1) +
        geom_line(aes(y=beta, group=sign), lwd=2)

    xlab_name <- if (plot_num_dropped)
        "Number of points removed" else "Proportion of points removed"
    plot <- plot + guides(color="none") + xlab(xlab_name)

    return(plot)
}


#' Plot influence scores, signals, and reruns.
#' @param param_infl `r docs$param_infl`
#' @param signal `r docs$signal`
#' @return A plot for the specified signal.
#'@export
PlotSignal <- function(param_infl, signal, ...) {
    stopifnot(class(param_infl) == "ParameterInferenceInfluence")
    influence_df <- GetSortedInfluenceDf(param_infl, signal$qoi$name)
    PlotInfluence(influence_df, signals=list(signal), ...)
}
