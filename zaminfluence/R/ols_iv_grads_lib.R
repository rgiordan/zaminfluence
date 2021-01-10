



######################################################3
# Ordinary least squares


# Compute the estimate, standard errors, and their derivatives for
# ordinary least squares regression.
GetRegressionSEDerivs <- function(x, y, beta, w0, se_group=NULL, testing=FALSE) {
  num_obs <- length(y)

  x_w <- x * w0
  xwx <- t(x_w) %*% x

  eps <- as.numeric(y - x %*% beta)

  # TODO: you can make this faster by factorizing xwx.
  if (!is.null(se_group)) {
    # Grouped ("robust") standard errors
    stop("Not implemented.")
  } else {
    # Un-grouped ("non-robust") standard errors
    sig2_hat <- sum(w0 * eps^2) / (num_obs - length(beta))

    sand_mat <- solve(xwx)
    se_mat <- sig2_hat * sand_mat

    ##############################
    # Derivatives

    dbetahat_dw <- solve(xwx, t(x * eps))

    # standard error partial derivatives
    dsig2_hat_dw_partial <- eps^2 / (num_obs - length(beta))
    dsig2_hat_dbeta <- -2 * colSums(w0 * eps * x) / (num_obs - length(beta))

    # Derivative of the diagonal of the sandwich matrix, which does not
    # depend on beta.
    # See notes for the definition of these terms and the derivation.
    R_x <- solve(xwx, t(x))
    dsand_mat_diag_dw_partial <- -1 * R_x * R_x

    dse_mat_diag_dw_partial <-
        dsand_mat_diag_dw_partial * sig2_hat +
        outer(diag(sand_mat), dsig2_hat_dw_partial)

    # combine with the chain rule.
    # Note that dsig2_hat_dbeta gets broadcasted like a column vector.
    dsig2_hat_dw <- colSums(dsig2_hat_dbeta * dbetahat_dw) + dsig2_hat_dw_partial
    dse_mat_diag_dw <-
        dse_mat_diag_dw_partial +
        outer(diag(sand_mat), colSums(dsig2_hat_dbeta * dbetahat_dw))

    dse_dw <- 0.5 * dse_mat_diag_dw / sqrt(diag(se_mat))

    # Specify return values
    ret_list <- list(
        se_mat=se_mat,
        se=sqrt(diag(se_mat)),
        sig2_hat=sig2_hat,
        dbetahat_dw=dbetahat_dw,
        dsig2_hat_dw=dsig2_hat_dw,
        dse_mat_diag_dw=dse_mat_diag_dw,
        dse_dw=dse_dw)

    if (testing) {
        # For testing and debugging, it's useful to get the intermediate results.
        ret_list$betahat <- solve(xwx, t(x_w) %*% y) %>% as.numeric()
        ret_list$sand_mat <- sand_mat
        ret_list$dsig2_hat_dbeta <- dsig2_hat_dbeta
        ret_list$dsig2_hat_dw_partial <- dsig2_hat_dw_partial
        ret_list$dsand_mat_diag_dw_partial <- dsand_mat_diag_dw_partial
        ret_list$dse_mat_diag_dw_partial <- dse_mat_diag_dw_partial
    }

    return(ret_list)
  }
}


GetRegressionVariables <- function(lm_res) {
  if (!(("x" %in% names(lm_res)) &
        ("y" %in% names(lm_res)))) {
      stop("You must run lm with the arguments x=TRUE and y=TRUE.")
  }

  x <- lm_res$x
  num_obs <- nrow(x)
  y <- as.numeric(lm_res$y)
  betahat <- as.numeric(lm_res$coefficients)
  if (is.null(lm_res$weights)) {
      w0 <- rep(1.0, num_obs)
  } else {
      w0 <- lm_res$weights
  }
  return(list(x=x, y=y, num_obs=num_obs, w0=w0, betahat=betahat))
}


#' Compute all influence scores for a regression.
#' @param lm_result The output of a call to lm.
#' @param se_group Optional, a vector of integers defining a standard error grouping.
#' @return A list containing the regression and influence result.
#' @export
ComputeRegressionInfluence <- function(lm_result, se_group=NULL) {
  reg_vars <- GetRegressionVariables(lm_result)
  reg_grad_list <- GetRegressionSEDerivs(
    x=reg_vars$x, y=reg_vars$y, beta=reg_vars$betahat,
    w0=reg_vars$w0, se_group=se_group)

  # Note that the standard errors may not match iv_res when using se_group.
  return(list(model_fit=lm_result,
              n_obs=reg_vars$num_obs,
              regressor_names=colnames(lm_result$x),
              grad_fun="GetIVSEDerivs",

              betahat=reg_vars$betahat,
              se=reg_grad_list$se,
              weights=reg_vars$w0,

              beta_grad=reg_grad_list$dbetahat_dw,
              se_grad=reg_grad_list$dse_dw)
  )
}


######################################################3
# Instrumental variables


# Compute the estimate, standard errors, and their derivatives for
# instrumental variables regression.
GetIVSEDerivs <- function(x, z, y, beta, w0, se_group=NULL, testing=FALSE) {
    # TODO: you can make this faster by factorizing zwx.
    if (!is.null(se_group)) {
        stop("Not implemented.")
    }

    num_obs <- length(y)

    z_w <- z * w0
    zwz <- t(z_w) %*% z
    zwx <- t(z_w) %*% x

    eps <- as.numeric(y - x %*% beta)
    sig2_hat <- sum(w0 * eps^2) / (num_obs - length(beta))

    zwx_inv_zwz <- solve(zwx, zwz)
    sand_mat <- solve(zwx, t(zwx_inv_zwz))
    se_mat <- sig2_hat * sand_mat

    ##############################
    # Derivatives

    dbetahat_dw <- solve(zwx, t(z * eps))

    # standard error partial derivatives
    dsig2_hat_dw_partial <- eps^2 / (num_obs - length(beta))
    dsig2_hat_dbeta <- -2 * colSums(w0 * eps * x) / (num_obs - length(beta))

    # Derivative of the diagonal of the sandwich matrix, which does not
    # depend on beta.
    # See notes for the definition of these terms and the derivation.
    AR_x <- zwx_inv_zwz %*% solve(t(zwx), t(x))
    R_z <- solve(zwx, t(z))
    dsand_mat_diag_dw_partial <- R_z * R_z - 2 * AR_x * R_z

    dse_mat_diag_dw_partial <-
        dsand_mat_diag_dw_partial * sig2_hat +
        outer(diag(sand_mat), dsig2_hat_dw_partial)

    # combine with the chain rule.
    # Note that dsig2_hat_dbeta gets broadcasted like a column vector.
    dsig2_hat_dw <- colSums(dsig2_hat_dbeta * dbetahat_dw) + dsig2_hat_dw_partial
    dse_mat_diag_dw <-
        dse_mat_diag_dw_partial +
        outer(diag(sand_mat), colSums(dsig2_hat_dbeta * dbetahat_dw))

    dse_dw <- 0.5 * dse_mat_diag_dw / sqrt(diag(se_mat))

    # Specify return values
    ret_list <- list(
        se_mat=se_mat,
        se=sqrt(diag(se_mat)),
        sig2_hat=sig2_hat,
        dbetahat_dw=dbetahat_dw,
        dsig2_hat_dw=dsig2_hat_dw,
        dse_mat_diag_dw=dse_mat_diag_dw,
        dse_dw=dse_dw)

    if (testing) {
        # For testing and debugging, it's useful to get the intermediate results.
        ret_list$betahat <- solve(zwx, t(z_w) %*% y) %>% as.numeric()
        ret_list$sand_mat <- sand_mat
        ret_list$dsig2_hat_dbeta <- dsig2_hat_dbeta # tested
        ret_list$dsig2_hat_dw_partial <- dsig2_hat_dw_partial # tested
        ret_list$dsand_mat_diag_dw_partial <- dsand_mat_diag_dw_partial # tested
        ret_list$dse_mat_diag_dw_partial <- dse_mat_diag_dw_partial # tested
    }

    return(ret_list)
}


GetIVVariables <- function(iv_res) {
  if (!(("x" %in% names(iv_res)) &
        ("y" %in% names(iv_res)))) {
      stop("You must run ivreg with the arguments x=TRUE and y=TRUE.")
  }

  x <- iv_res$x$regressors
  num_obs <- nrow(x)
  z <- iv_res$x$instruments
  y <- as.numeric(iv_res$y)
  betahat <- as.numeric(iv_res$coefficients)
  if (is.null(iv_res$weights)) {
      w0 <- rep(1.0, num_obs)
  } else {
      w0 <- iv_res$weights
  }
  return(list(x=x, y=y, z=z, num_obs=num_obs, w0=w0, betahat=betahat))
}


#' Compute all influence scores for an IV regression.
#' @param iv_res The output of an IV regression computed with AER::ivreg.
#' @param se_group Optional, a vector of integers defining a standard error grouping.
#' @return A list containing the regression and influence result.
#' @export
ComputeIVRegressionInfluence <- function(iv_res, se_group=NULL) {
    iv_vars <- GetIVVariables(iv_res)
    iv_grad_list <- GetIVSEDerivs(
      x=iv_vars$x, z=iv_vars$z, y=iv_vars$y,
      beta=iv_vars$betahat, w0=iv_vars$w0, se_group=se_group)

    # Note that the standard errors may not match iv_res when using se_group.
    return(list(model_fit=iv_res,
                n_obs=iv_vars$num_obs,
                regressor_names=colnames(iv_res$x$regressors),
                grad_fun="GetIVSEDerivs",

                betahat=iv_vars$betahat,
                se=iv_grad_list$se,
                weights=iv_vars$w0,

                beta_grad=iv_grad_list$dbetahat_dw,
                se_grad=iv_grad_list$dse_dw)
    )
}


#' Compute all influence scores for an IV regression.
#' @param iv_res The output of an IV regression computed with AER::ivreg.
#' @param se_group Optional, a vector of integers defining a standard error grouping.
#' @return The standard error matrix.
#' @export
ComputeIVRegressionErrorCovariance <- function(iv_res, se_group=NULL) {
  iv_vars <- GetIVVariables(iv_res)
  iv_grad_list <- GetIVSEDerivs(
    x=iv_vars$x, z=iv_vars$z, y=iv_vars$y,
    beta=iv_vars$betahat, w0=iv_vars$w0, se_group=se_group)
  return(iv_grad_list$se_mat / iv_vars$num_obs)
}



####################################################################
# The remaining functions are common to ordinary and IV regression.

#' Compute the influence functions for all regressors given a model fit.
#' @param model_fit A model fit (currently from lm or AER::iv_reg).
#' @param se_group Optional, a vector of integers defining a standard error grouping.
#' @return A list containing the regression and influence result.
#' @export
ComputeModelInfluence <- function(model_fit, se_group=NULL) {
  valid_classes <- c("lm", "ivreg")
  model_class <- class(model_fit)
  if (!(model_class %in% valid_classes)) {
    stop(sprintf("The class of `model_fit` must be one of %s",
                 paste(valid_classes, collapse=", ")))
  }
  if (model_class == "lm") {
    return(ComputeRegressionInfluence(model_fit, se_group))
  } else if (model_class == "ivreg") {
    return(ComputeIVRegressionInfluence(model_fit, se_group))
  } else {
    # Redundant, so sue me.
    stop(sprint("Unknown model class %s", model_class))
  }
}


#' @export
AppendIntervalColumns <- function(grad_df, sig_num_ses) {
  grad_df[["beta_pzse_grad"]] <-
    grad_df[["beta_grad"]] + sig_num_ses * grad_df[["se_grad"]]
  grad_df[["beta_mzse_grad"]] <-
    grad_df[["beta_grad"]] - sig_num_ses * grad_df[["se_grad"]]
  base_vals <- attr(grad_df, "base_vals")
  base_vals_names <- names(base_vals)
  base_vals <- c(base_vals,
                 base_vals["beta"] + sig_num_ses * base_vals["se"],
                 base_vals["beta"] - sig_num_ses * base_vals["se"])
  names(base_vals) <- c(base_vals_names,
                        "beta_pzse",
                        "beta_mzse")
  attr(grad_df, "base_vals") <- base_vals
  attr(grad_df, "sig_num_ses") <- sig_num_ses
  return(grad_df)
}


#' @export
GetTargetRegressorGrads <- function(reg_infl, target_regressor,
                                    sig_num_ses=qnorm(0.975)) {
    target_index <- which(reg_infl$regressor_names == target_regressor)
    if (length(target_index) != 1) {
        stop("Error finding target regressor in the regression.")
    }

    # The reg_infl$*_grad columns are derivatives with respect to each
    # observation at `weights`.  They are converted to derivatives with respect
    # to a weight scaled to be one at inclusion and zero at exclusion by the
    # chain rule.
    grad_df <-
      data.frame(
        row=1:reg_infl$n_obs,
        weights=reg_infl$weights,
        se_grad=reg_infl$weights * reg_infl$se_grad[target_index,],
        beta_grad=reg_infl$weights * reg_infl$beta_grad[target_index, ],
        obs_per_row=1)

    attr(grad_df, "n_obs") <- reg_infl$n_obs
    attr(grad_df, "n_grad_rows") <- reg_infl$n_obs
    attr(grad_df, "obs_per_row_col") <- "obs_per_row"
    base_vals <- c(
        reg_infl$betahat[target_index],
        reg_infl$se[target_index])
    names(base_vals) <- c("beta", "se")
    attr(grad_df, "base_vals") <- base_vals
    attr(grad_df, "target_regressor") <- target_regressor
    attr(grad_df, "target_index") <- target_index
    attr(grad_df, "data_row_cols") <- "row"

    grad_df <- AppendIntervalColumns(grad_df, sig_num_ses=sig_num_ses)

    return(grad_df)
}


#' @export
CopyGradAttributes <- function(
    new_df, grad_df,
    attrs=c("n_obs", "obs_per_row_col", "base_vals",
            "target_regressor", "target_index",
            "sig_num_ses", "data_row_cols", "n_grad_rows")) {
    for (a in attrs) {
        attr(new_df, a) <- attr(grad_df, a)
    }
    return(new_df)
}
