library(ggplot2)
library(latex2exp)


#'@export
AppendTargetRegressorInfluence <- function(reg_infl, target_regressor,
                                           sig_num_ses=qnorm(0.975)) {
    if (is.null(reg_infl[["targets"]])) {
        reg_infl$targets <- list()
    }
    target_index <- which(reg_infl$regressor_names == target_regressor)
    if (length(target_index) != 1) {
        stop("Error finding target regressor in the regression.")
    }

    se_grad <- reg_infl$weights * reg_infl$se_grad[target_index,]
    beta_grad <- reg_infl$weights * reg_infl$beta_grad[target_index, ]
    betahat <- reg_infl$betahat[target_index]
    sehat <- reg_infl$se[target_index]
    n_obs <- reg_infl$n_obs

    target_list <- list(
        target_index=target_index,
        target_regressor=target_regressor,
        beta=ProcessInfluenceVector(
            infl=beta_grad,
            base_value=betahat,
            num_obs=n_obs),
        beta_mzse=ProcessInfluenceVector(
            infl=beta_grad - sig_num_ses * se_grad,
            base_value=betahat - sig_num_ses * sehat,
            num_obs=n_obs),
        beta_pzse=ProcessInfluenceVector(
            infl=beta_grad + sig_num_ses * se_grad,
            base_value=betahat + sig_num_ses * sehat,
            num_obs=n_obs),
        sig_num_ses=sig_num_ses
    )

    reg_infl$targets[[target_regressor]] <- target_list
    return(reg_infl)
}


#' @export
GetRegressionSignals <- function(estimator_infl) {
    betahat <- estimator_infl$beta$base_value
    beta_mzse <- estimator_infl$beta_mzse$base_value
    beta_pzse <- estimator_infl$beta_pzse$base_value

    sign_label <- "sign"
    sig_label <- "significance"
    both_label <- "sign and significance"

    signals <- list()
    # signals$target_index <- estimator_infl$target_index
    # signals$sig_num_ses <- estimator_infl$sig_num_ses
    signals$target_regressor <- estimator_infl$target_regressor
    signals$sign <- list(metric="beta", signal=-1 * betahat, change=sign_label)

    is_significant <- sign(beta_mzse) == sign(beta_mzse)
    if (is_significant) {
        if (beta_mzse >= 0) { # then beta_pzse > 0 too because significant
            signals$sig <- list(
              metric="beta_mzse", signal=-1 * beta_mzse, change=sig_label)
            signals$both <- list(
              metric="beta_pzse", signal=-1 * beta_pzse, change=both_label)
        } else if (beta_pzse < 0) { # then beta_mzse < 0 too because significant
            signals$sig <- list(
              metric="beta_pzse", signal=-1 * beta_pzse, change=sig_label)
            signals$both <- list(
              metric="beta_mzse", signal=-1 * beta_mzse, change=both_label)
        } else {
            stop("Impossible for a significant result")
        }
    } else { # Not significant.  Choose to change the interval endpoint which
             # is closer.
        if (abs(beta_mzse) >= abs(beta_pzse)) {
            signals$sig <- list(
              metric="beta_pzse", signal=-1 * beta_pzse, change=sig_label)
        } else  {
            signals$sig <- list(
              metric="beta_mzse", signal=-1 * beta_mzse, change=sig_label)
        }

        if (betahat >= 0) {
            signals$both <- list(
              metric="beta_pzse", signal=-1 * beta_pzse, change=both_label)
        } else {
            signals$both <- list(
              metric="beta_mzse", signal=-1 * beta_mzse, change=both_label)
        }
    }

    # Get the APIP for all the signals
    for (target in c("sign", "sig", "both")) {
        reg_signal <- signals[[target]]
        apip <- GetAPIP(
            infl_lists=estimator_infl[[reg_signal$metric]],
            signal=reg_signal$signal)
        signals[[target]]$apip <- apip
    }

    return(signals)
}


#' @export
GetSignalDataFrame <- function(reg_signal) {
    data.frame(
        metric=reg_signal$metric,
        change=reg_signal$change,
        signal=reg_signal$signal,
        num_removed=reg_signal$apip$n,
        prop_removed=reg_signal$apip$prop)
}


#' @export
RerunForTargetChanges <- function(signals, reg_infl) {
  for (target in c("sign", "sig", "both")) {
      reg_signal <- signals[[target]]
      w_bool <- GetWeightVector(
          drop_inds=reg_signal$apip$inds,
          num_obs=reg_infl$n_obs,
          bool=TRUE)

      rerun <- RerunRegression(reg_infl$model_fit, w_bool=w_bool)

      # Save the whole rerun
      signals[[target]]$rerun <- rerun

      regressor_infl <- reg_infl$targets[[signals$target_regressor]]
      # Make a nice dataframe with the targeted regressor
      betahat <- rerun$betahat[[regressor_infl$target_index]]
      sehat <- rerun$se[[regressor_infl$target_index]]
      sig_num_ses <- regressor_infl$sig_num_ses

      signals[[target]]$rerun_df <-
          GetSignalDataFrame(reg_signal) %>%
              mutate(
                  betahat_refit=betahat,
                  beta_mzse_refit=betahat - sig_num_ses * sehat,
                  beta_pzse_refit=betahat + sig_num_ses * sehat,
                  betahat_orig=regressor_infl$beta$base_value,
                  beta_mzse_orig=regressor_infl$beta_mzse$base_value,
                  beta_pzse_orig=regressor_infl$beta_pzse$base_value)
  }
  return(signals)
}



#' @export
GetSortedInfluenceDf <- function(parameter_infl, sorting_qoi_name, max_num_obs=Inf) {
    qoi_for_sorting <- parameter_infl[[sorting_qoi_name]]

    GetQOIDf <- function(infl_sign) {
        ordered_inds <- qoi_for_sorting[[infl_sign]]$infl_inds
        if (max_num_obs < length(ordered_inds)) {
            ordered_inds <- ordered_inds[1:max_num_obs]
        }
        qoi_df <- data.frame(
            num_dropped=c(0, 1:length(ordered_inds))) %>%
            mutate(prop_dropped=num_dropped / qoi_for_sorting[[infl_sign]]$num_obs,
                   sign=infl_sign)
        for (qoi_name in c("beta", "beta_mzse", "beta_pzse")) {
            base_value <- parameter_infl[[qoi_name]]$base_value
            infl_sorted <- parameter_infl[[qoi_name]]$infl[ordered_inds]
            qoi_df[[qoi_name]] <- base_value + cumsum(c(0, infl_sorted))
        }
        return(qoi_df)
    }

    qoi_df <-
        bind_rows(GetQOIDf("pos"), GetQOIDf("neg")) %>%
        mutate(sorted_by=sorting_qoi_name)

    return(qoi_df)
}



#' @export
PlotInfluence <- function(influence_df,
                          plot_num_dropped=FALSE,
                          apip_max=NULL,
                          reg_signals=NULL,
                          include_y_zero=TRUE) {

    PlotRegSignal <- function(plot, reg_signal) {
        alpha_type <- if (plot_num_dropped) "n" else "prop"
        alpha <- reg_signal$apip[[alpha_type]]
        if (is.null(apip_max) || (!is.null(apip_max) && alpha <= apip_max)) {
            plot <- plot + geom_vline(aes(xintercept=!!alpha,
                                          linetype=!!reg_signal$change))
            if (!is.null(reg_signal$rerun_df)) {
                rerun_df <- reg_signal$rerun_df
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
    if (!is.null(reg_signals)) {
        for (reg_signal in reg_signals) {
            plot <- PlotRegSignal(plot,  reg_signal)
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


#'@export
PlotSignal <- function(parameter_infl, reg_signal, ...) {
    influence_df <- GetSortedInfluenceDf(parameter_infl, reg_signal$metric)
    PlotInfluence(influence_df, reg_signals=list(reg_signal), ...)
}
