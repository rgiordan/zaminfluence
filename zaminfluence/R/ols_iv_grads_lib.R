######################################################3
# Ordinary least squares

# Sum the rows of the matrix mat according to the 0-based integers group.
GroupedSum <- function(mat, group) {
  stopifnot(length(group) == nrow(mat))
  group_rows <- split(1:length(group), group)
  summed_mat <- lapply(group_rows,
         function(rows) { colSums(mat[rows, , drop=FALSE])}) %>%
      do.call(rbind, .)
  return(summed_mat)
}


# Repeat grouped sum rows to match the shape of the original data.
# s_mat might be the output of GroupedSum.
ExpandGroupedSum <- function(s_mat, group) {
  return(s_mat[group + 1, ])
}


# Compute the estimate, standard errors, and their derivatives for
# ordinary least squares regression.
# See the file inst/regression_derivatives.pdf for the derivation of
# the derivatives computed by this function.
GetRegressionSEDerivs <- function(x, y, beta, w0,
                                  se_group=NULL,
                                  testing=FALSE,
                                  compute_derivs=TRUE) {
  num_obs <- length(y)

  eps <- as.numeric(y - x %*% beta)
  x_w <- x * w0

  xwx_qr <- qr(t(x_w) %*% x)
  betahat <- solve(xwx_qr, t(x_w) %*% y) %>% as.numeric()

  if (compute_derivs) {
    # This derivative is common to both standard error methods.
    dbetahat_dw <- solve(xwx_qr, t(x * eps))
  } else {
    dbetahat_dw <- NA
  }

  if (!is.null(se_group)) {
    ##############################################
    # Derivatives for grouped standard errors:
    s_mat <- GroupedSum(x * w0 * eps, se_group)

    # colMeans(s_mat) is zero at the weights used for regression, but include
    # it so we can test partial derivatives.
    s_mat <- s_mat - rep(colMeans(s_mat), each=nrow(s_mat))

    # Check for well-formedness of the groups.
    all(as.numeric(row.names(s_mat)) ==
        unique(se_group) %>% sort()) %>%
      stopifnot()
    all(as.numeric(row.names(s_mat)) ==
        (min(se_group):max(se_group))) %>%
      stopifnot()
    num_groups <- nrow(s_mat)

    # v is for variance, so I take the average, but we then need to multiply
    # again by num_groups to get the standard error variance.
    v_mat <- t(s_mat) %*% s_mat / num_groups

    # Covariance matrix
    xwx_inv_vmat <- solve(xwx_qr, v_mat)
    se_mat <- solve(xwx_qr, t(xwx_inv_vmat)) * num_groups

    if (compute_derivs) {
      # Covariance matrix derivatives.

      # First the partials wrt the weights.
      # s_mat_expanded is s_mat with rows repeated to match the shape of z.
      xwx_inv_s_mat <- solve(xwx_qr, t(s_mat))
      s_mat_expanded <- ExpandGroupedSum(s_mat, se_group)
      xwx_inv_s_mat_expanded <-
        ExpandGroupedSum(t(xwx_inv_s_mat), se_group) %>% t()
      #xwx_inv_s_mat_expanded <- solve(xwx_qr, t(s_mat_expanded))

      # To understand the next line,
      # note that dbetahat_dw = solve(xwx_qr, t(x * eps))
      ddiag_semat_dw_partial <-
        2 * xwx_inv_s_mat_expanded * dbetahat_dw -
        2 * (solve(xwx_qr, t(x))) * (se_mat %*% t(x))

      # Second, the partials through the beta dependence.
      beta_dim <- ncol(x)
      # ddiag_semat_dbeta_partial will be beta_dim x beta_dim; the columns
      # correspond to the entries of betahat.
      ddiag_semat_dbeta_partial <- matrix(NA, nrow(se_mat), nrow(dbetahat_dw))
      for (d in 1:beta_dim) {
        xi_d <- GroupedSum(-1 * x_w * x[, d], se_group)
        xi_d <- xi_d - rep(colMeans(xi_d), each=nrow(xi_d))
        ddiag_semat_dbeta_partial[, d] <-
          2 * xwx_inv_s_mat %*%
              t(solve(xwx_qr, t(xi_d))) %>% diag()
      }

      ddiag_semat_dw <-
        ddiag_semat_dw_partial + ddiag_semat_dbeta_partial %*% dbetahat_dw
      dse_dw <- 0.5 * ddiag_semat_dw / sqrt(diag(se_mat))
    } else {
      # Don't compute any derivative information.
      dse_mat_diag_dw <- NA
      dse_dw <- NA
      ddiag_semat_dw_partial <- NA
      s_mat_expanded <- NA
      ddiag_semat_dbeta_partial <- NA
      ddiag_semat_dw <- NA
    }

    # Specify return values
    ret_list <- list(
        betahat=betahat,
        se_mat=se_mat,
        se=sqrt(diag(se_mat)),
        dbetahat_dw=dbetahat_dw,
        dse_mat_diag_dw=ddiag_semat_dw,
        dse_dw=dse_dw
      )

    if (testing) {
        # For testing and debugging, it's useful to get the
        # intermediate results.
        ret_list$v_mat <- v_mat
        ret_list$ddiag_semat_dw_partial <- ddiag_semat_dw_partial
        ret_list$s_mat <- s_mat
        ret_list$s_mat_expanded <- s_mat_expanded
        ret_list$ddiag_semat_dbeta_partial <- ddiag_semat_dbeta_partial
    }
    return(ret_list)
  } else {
    ##############################################
    # Derivatives for un-grouped standard errors:

    sig2_hat <- sum(w0 * eps^2) / (num_obs - length(beta))

    sand_mat <- solve(xwx_qr)
    se_mat <- sig2_hat * sand_mat

    if (compute_derivs) {
      ##############################
      # Derivatives

      # standard error partial derivatives
      dsig2_hat_dw_partial <- eps^2 / (num_obs - length(beta))
      dsig2_hat_dbeta <- -2 * colSums(w0 * eps * x) / (num_obs - length(beta))

      # Derivative of the diagonal of the sandwich matrix, which does not
      # depend on beta.
      # See notes for the definition of these terms and the derivation.
      R_x <- solve(xwx_qr, t(x))
      dsand_mat_diag_dw_partial <- -1 * R_x * R_x

      dse_mat_diag_dw_partial <-
          dsand_mat_diag_dw_partial * sig2_hat +
          outer(diag(sand_mat), dsig2_hat_dw_partial)

      # combine with the chain rule.
      # Note that dsig2_hat_dbeta gets broadcasted like a column vector.
      dsig2_hat_dw <- colSums(
        dsig2_hat_dbeta * dbetahat_dw) + dsig2_hat_dw_partial
      dse_mat_diag_dw <-
          dse_mat_diag_dw_partial +
          outer(diag(sand_mat), colSums(dsig2_hat_dbeta * dbetahat_dw))

      dse_dw <- 0.5 * dse_mat_diag_dw / sqrt(diag(se_mat))
    } else {
      dsig2_hat_dw <- NA
      dse_mat_diag_dw <- NA
      dse_dw <- NA
      dsig2_hat_dbeta <- NA
      dsig2_hat_dw_partial <- NA
      dsand_mat_diag_dw_partial <- NA
      dse_mat_diag_dw_partial <- NA
    }

    # Specify return values
    ret_list <- list(
      betahat=betahat,
      se_mat=se_mat,
      se=sqrt(diag(se_mat)),
      sig2_hat=sig2_hat,
      dbetahat_dw=dbetahat_dw,
      dsig2_hat_dw=dsig2_hat_dw,
      dse_mat_diag_dw=dse_mat_diag_dw,
      dse_dw=dse_dw)

    if (testing) {
        # For testing and debugging, it's useful to get the intermediate results.
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
# See the file inst/regression_derivatives.pdf for the derivation of
# the derivatives computed by this function.
GetIVSEDerivs <- function(x, z, y, beta, w0, se_group=NULL,
                          compute_derivs=TRUE, testing=FALSE) {
  num_obs <- length(y)

  eps <- as.numeric(y - x %*% beta)

  z_w <- z * w0
  zwz <- t(z_w) %*% z
  #zwx <- t(z_w) %*% x
  zwx_qr <- qr(t(z_w) %*% x)

  betahat <- solve(zwx_qr, t(z_w) %*% y) %>% as.numeric()

  if (compute_derivs) {
    # This derivative is the same for both SE methods.
    dbetahat_dw <- solve(zwx_qr, t(z * eps))
  } else {
    dbetahat_dw <- NA
  }

  if (!is.null(se_group)) {
    ##############################################
    # Derivatives for grouped standard errors:
    s_mat <- GroupedSum(z * eps * w0, se_group)

    # colMeans(s_mat) is zero at the weights used for regression, but include
    # it so we can test partial derivatives.
    s_mat <- s_mat - rep(colMeans(s_mat), each=nrow(s_mat))

    # Check for well-formedness of the groups.
    all(as.numeric(row.names(s_mat)) ==
        unique(se_group) %>% sort()) %>%
      stopifnot()
    all(as.numeric(row.names(s_mat)) ==
        (min(se_group):max(se_group))) %>%
      stopifnot()
    num_groups <- nrow(s_mat)

    # v is for variance, so I take the average, but we then need to multiply
    # again by num_groups to get the standard error variance.
    v_mat <- t(s_mat) %*% s_mat / num_groups

    # Covariance matrix
    zwx_inv_vmat <- solve(zwx_qr, v_mat)
    se_mat <- solve(zwx_qr, t(zwx_inv_vmat)) * num_groups

    if (compute_derivs) {
      # Covariance matrix derivatives.

      # First the partials wrt the weights.

      #s_mat_expanded <- ExpandGroupedSum(s_mat, se_group)

      zwx_inv_s_mat <- solve(zwx_qr, t(s_mat))
      #zwx_inv_s_mat_expanded <- solve(zwx_qr, t(s_mat_expanded))
      zwx_inv_s_mat_expanded <-
        ExpandGroupedSum(t(zwx_inv_s_mat), se_group) %>% t()

      # To understand the next line, note that
      # dbetahat_dw = solve(zwx_qr, t(z * eps))
      ddiag_semat_dw_partial <-
        2 * zwx_inv_s_mat_expanded * dbetahat_dw -
        2 * (solve(zwx_qr, t(z))) * (se_mat %*% t(x))

      # Second, the partials through the beta dependence.
      beta_dim <- ncol(x)
      # ddiag_semat_dbeta_partial will be beta_dim x beta_dim; the columns
      # correspond to the entries of betahat.
      ddiag_semat_dbeta_partial <- matrix(NA, nrow(se_mat), nrow(dbetahat_dw))
      for (d in 1:beta_dim) {
        xi_d <- GroupedSum(-1 * z_w * x[, d], se_group)
        xi_d <- xi_d - rep(colMeans(xi_d), each=nrow(xi_d))
        ddiag_semat_dbeta_partial[, d] <-
          2 * zwx_inv_s_mat %*%
              t(solve(zwx_qr, t(xi_d))) %>% diag()
      }

      ddiag_semat_dw <-
        ddiag_semat_dw_partial + ddiag_semat_dbeta_partial %*% dbetahat_dw
      dse_dw <- 0.5 * ddiag_semat_dw / sqrt(diag(se_mat))
    } else {
      ddiag_semat_dw <- NA
      dse_dw <- NA
      ddiag_semat_dw_partial <- NA
      ddiag_semat_dbeta_partial <- NA
    }

    # Specify return values
    ret_list <- list(
        betahat=betahat,
        se_mat=se_mat,
        se=sqrt(diag(se_mat)),
        dbetahat_dw=dbetahat_dw,
        dse_mat_diag_dw=ddiag_semat_dw,
        dse_dw=dse_dw
      )

    if (testing) {
        # For testing and debugging, it's useful to get the intermediate results.
        ret_list$v_mat <- v_mat
        ret_list$ddiag_semat_dw_partial <- ddiag_semat_dw_partial
        ret_list$s_mat <- s_mat
        ret_list$ddiag_semat_dbeta_partial <- ddiag_semat_dbeta_partial
    }
    return(ret_list)

  } else {
    ##############################################
    # Derivatives for un-grouped standard errors:

    sig2_hat <- sum(w0 * eps^2) / (num_obs - length(beta))

    zwx_inv_zwz <- solve(zwx_qr, zwz)
    sand_mat <- solve(zwx_qr, t(zwx_inv_zwz))
    se_mat <- sig2_hat * sand_mat

    if (compute_derivs) {
      ##############################
      # Derivatives

      # standard error partial derivatives
      dsig2_hat_dw_partial <- eps^2 / (num_obs - length(beta))
      dsig2_hat_dbeta <- -2 * colSums(w0 * eps * x) / (num_obs - length(beta))

      # Derivative of the diagonal of the sandwich matrix, which does not
      # depend on beta.
      # See notes for the definition of these terms and the derivation.
      #AR_x <- zwx_inv_zwz %*% solve(t(zwx), t(x))
      AR_x <- sand_mat %*% t(x)
      R_z <- solve(zwx_qr, t(z))
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
    } else {
      dsig2_hat_dw <- NA
      dse_mat_diag_dw <- NA
      dse_dw <- NA
      dsig2_hat_dbeta <- NA
      dsig2_hat_dw_partial <- NA
      dsand_mat_diag_dw_partial <- NA
      dse_mat_diag_dw_partial <- NA
    }

    # Specify return values
    ret_list <- list(
        betahat=betahat,
        se_mat=se_mat,
        se=sqrt(diag(se_mat)),
        sig2_hat=sig2_hat,
        dbetahat_dw=dbetahat_dw,
        dsig2_hat_dw=dsig2_hat_dw,
        dse_mat_diag_dw=dse_mat_diag_dw,
        dse_dw=dse_dw)

    if (testing) {
        # For testing and debugging, it's useful to get the intermediate results.
        ret_list$sand_mat <- sand_mat
        ret_list$dsig2_hat_dbeta <- dsig2_hat_dbeta
        ret_list$dsig2_hat_dw_partial <- dsig2_hat_dw_partial
        ret_list$dsand_mat_diag_dw_partial <- dsand_mat_diag_dw_partial
        ret_list$dse_mat_diag_dw_partial <- dse_mat_diag_dw_partial
    }
    return(ret_list)
  }
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


#' Compute the standard error matrix for an IV regression.
#' @param iv_res The output of an IV regression computed with AER::ivreg.
#' @param se_group Optional, a vector of integers defining a standard error grouping.
#' @return The standard error matrix.
#' @export
ComputeIVRegressionErrorCovariance <- function(iv_res, se_group=NULL) {
  iv_vars <- GetIVVariables(iv_res)
  iv_grad_list <- GetIVSEDerivs(
    x=iv_vars$x, z=iv_vars$z, y=iv_vars$y,
    beta=iv_vars$betahat, w0=iv_vars$w0, se_group=se_group)
  return(iv_grad_list$se_mat)
}



####################################################################
# The remaining functions are common to ordinary and IV regression.

#' Compute the influence functions for all regressors given a model fit.
#' @param model_fit A model fit (currently from lm or AER::iv_reg).
#' @param se_group Optional, a vector of integers defining a standard error
#  grouping.
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
