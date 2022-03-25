library(tidyr)
library(purrr)
library(latex2exp)
library(ggplot2)
library(tibble)


################################################################################
# Plotting and visualization functions


#' Summarize the values of each QOI for each parameter for a given model_fit.
#'@param model_fit `r docs$model_fit`
#'@param param_infls A list of ParameterInferenceInfluence objects.
#'@return A dataframe summarizing the values of all quantities of interest
#' in `param_infls` for `model_fit`.
#'@export
GetModelFitInferenceDataframe <- function(model_fit, param_infls) {
    if (is.null(model_fit)) {
      return(data.frame())
    }
    stopifnot(class(model_fit) == "ModelFit")
    stopifnot(all(names(param_infls) %in% model_fit$parameter_names))

    GetParameterInferenceDataframe <-
      function(model_fit, target_index, sig_num_ses) {
        GetInferenceQOIs(param=model_fit$param[target_index],
                         se=model_fit$se[target_index],
                         sig_num_ses=sig_num_ses) %>%
            purrr::imap_dfr(~ data.frame(metric=.y, value=.x))
    }

    summary_df <- data.frame()
    AppendRow <- function(row) bind_rows(summary_df, row)
    for (param_name in names(param_infls)) {
        param_infl <- param_infls[[param_name]]
        # We checked above that each parameter name is found.
        target_index <- GetParameterIndex(model_fit, param_name)
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


# The signals and reruns are expected to have a matching list structure,
# which we enforce with this function.
ValidateSignalsAndReruns <- function(signals, reruns) {
  stopifnot(setequal(names(reruns), names(signals)))
  for (target_param_name in names(reruns)) {
      param_reruns  <- reruns[[target_param_name]]
      param_signals  <- signals[[target_param_name]]
      stopifnot(names(param_reruns) == names(param_signals))
      for (signal_name in names(param_reruns)) {
          rerun <- param_reruns[[signal_name]]
          signal <- param_signals[[signal_name]]
          stopifnot(class(signal) == "QOISignal")
          if (signal$apip$success) {
            stopifnot(class(rerun) == "ModelFit")
          } else {
            stopifnot(is.null(rerun))
          }
      }
  }
}

#' For signals and reruns lists, produce a dataframe summarizing the two.
#'@export
GetSignalsAndRerunsDataframe <- function(signals, reruns, model_grads) {
  ValidateSignalsAndReruns(signals, reruns)

  reruns_dfs <- map_depth(
    reruns, 2, ~ GetModelFitInferenceDataframe(., model_grads$param_infls))

  rerun_df <-
      tibble(list=reruns_dfs) %>%
      mutate(target_param_name=names(list)) %>%
      unnest_longer(col=list, indices_to="target_signal") %>%
      unnest(list)

  signal_dfs <-
      map_depth(signals, 2, ~ data.frame(
          description=.$description, n_drop=.$apip$n, prop_drop=.$apip$prop,
          target_qoi=.$qoi$name))
  signal_df <-
      tibble(list=signal_dfs) %>%
      mutate(target_param_name=names(list)) %>%
      unnest_longer(col=list, indices_to="target_signal") %>%
      unnest(list)

  return(inner_join(
    rerun_df, signal_df, by=c("target_param_name", "target_signal")))
}


#' Produce an influence dataframe suitable for visualization.
#' @param param_infl `r docs$param_infl`
#' @param sorting_qoi_name The name of a QOI ("param", "param_mzse", or
#' "param_pzse") whose influence is used to sort the dataframe.
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
        for (qoi_name in c("param", "param_mzse", "param_pzse")) {
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
PlotInfluenceDf <- function(influence_df, signal, rerun_vals=NULL,
                            plot_num_dropped=FALSE,
                            apip_max=NULL,
                            include_y_zero=TRUE) {
    influence_df$alpha <-
        if (plot_num_dropped) influence_df$num_dropped else
            influence_df$prop_dropped

    if (!is.null(apip_max)) {
      influence_df <- filter(influence_df, alpha  <= apip_max)
    }

    plot <- ggplot(influence_df, aes(x=alpha))
    if (include_y_zero) {
        plot <-
            plot +
            geom_line(aes(y=0.0), col="gray50")
    }

    base_param <- filter(influence_df, alpha == 0) %>% pull("param") %>% unique()
    stopifnot(length(base_param) == 1)
    plot <-
        plot +
        geom_line(aes(y=!!base_param), col="blue", lwd=2) +
        geom_ribbon(aes(
            ymin=param_mzse,
            ymax=param_pzse,
            group=sign),
            fill="blue", color=NA, alpha=0.1) +
        geom_line(aes(y=param, group=sign), lwd=2)

    xlab_name <- if (plot_num_dropped)
        "Number of points removed" else "Proportion of points removed"
    plot <- plot + guides(color="none") + xlab(xlab_name)

    # Plot the signal
    if (signal$apip$success) {
      alpha_type <- if (plot_num_dropped) "n" else "prop"
      alpha <- signal$apip[[alpha_type]]
      if (is.null(apip_max) || (!is.null(apip_max) && alpha <= apip_max)) {
          plot <- plot + geom_vline(aes(xintercept=!!alpha,
                                        linetype=!!signal$description)) +
                  guides(linetype=guide_legend(title="Change type"))
      }
    }

    if (!is.null(rerun_vals)) {
      errorbar_width <- diff(range(influence_df$alpha)) / 50
      plot <-
          plot +
          geom_errorbar(aes(
              x=!!alpha,
              ymin=rerun_vals$param_mzse,
              ymax=rerun_vals$param_pzse),
              data=NULL,
              width=errorbar_width,
              lwd=1.5) +
          geom_point(aes(x=!!alpha, y=rerun_vals$param),
                     data=NULL,
                     shape=8)
    }

    return(plot)
}


#' Plot influence scores, signals, and reruns.
#' @param param_infl `r docs$param_infl`
#' @param signal `r docs$signal`
#' @return A plot for the specified signal.
#'@export
PlotSignal <- function(model_grads, signals, parameter_name, target_signal,
                       reruns=NULL, ...) {
    stopifnot(class(model_grads) == "ModelGrads")
    stopifnot(parameter_name %in% names(model_grads$param_infls))
    param_infl <- model_grads$param_infls[[parameter_name]]

    stopifnot(all(names(signals) %in% names(model_grads$param_infls)))

    stopifnot(parameter_name %in% names(signals))
    stopifnot(target_signal %in%  names(signals[[parameter_name]]))
    signal <- signals[[parameter_name]][[target_signal]]
    stopifnot(class(signal) == "QOISignal")

    rerun_vals <- NULL
    if (!is.null(reruns)) {
      ValidateSignalsAndReruns(signals, reruns)
      rerun <- reruns[[parameter_name]][[target_signal]]
      if (!is.null(rerun)) {
        rerun_vals <- GetParameterInferenceQOIs(
          rerun, parameter_name, sig_num_ses=param_infl$sig_num_ses)
      }
    }

    influence_df <- GetSortedInfluenceDf(param_infl, signal$qoi$name)
    plot <- PlotInfluenceDf(influence_df, signal, rerun_vals=rerun_vals, ...)
    return(plot)
}
