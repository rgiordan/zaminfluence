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


# Run regression using the QR decomposition.
GetRegressionCoefficients <- function(x, y, w0) {
  x_w <- x * w0
  xwx_qr <- qr(t(x_w) %*% x)
  betahat <- solve(xwx_qr, t(x_w) %*% y) %>% as.numeric()
  return(list(betahat=betahat, x_w=x_w, xwx_qr=xwx_qr))
}


# Compute the estimate, standard errors, and their derivatives for
# ordinary least squares regression.
# See the file inst/regression_derivatives.pdf for the derivation of
# the derivatives computed by this function.
GetRegressionSEDerivs <- function(x, y, beta, w0,
                                  se_group=NULL,
                                  testing=FALSE,
                                  compute_derivs=TRUE) {

  reg_coeff <- GetRegressionCoefficients(x, y, w0)
  betahat <- reg_coeff$betahat
  xwx_qr <- reg_coeff$xwx_qr
  x_w <- reg_coeff$x_w

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

    # Enforces that the grouping variable is sequential and zero-indexed.
    se_group <- as.numeric(factor(se_group)) - 1

    s_mat <- GroupedSum(x * w0 * eps, se_group)

    # colMeans(s_mat) is zero at the weights used for regression, but include
    # it so we can test partial derivatives.
    s_mat <- s_mat - rep(colMeans(s_mat), each=nrow(s_mat))

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


# Extract the relevant variables from the output of `lm()`
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




######################################################3
# Instrumental variables

# Get the IV estimates from regressors x, instruemnts z, responses y, and
# weights w0.
GetIVCoefficients <- function(x, z, y, w0) {
  z_w <- z * w0
  zwz <- t(z_w) %*% z
  zwx_qr <- qr(t(z_w) %*% x)

  betahat <- solve(zwx_qr, t(z_w) %*% y) %>% as.numeric()

  return(list(
    z_w=z_w,
    zwz=zwz,
    zwx_qr=zwx_qr,
    betahat=betahat))
}

# Compute the estimate, standard errors, and their derivatives for
# instrumental variables regression.
# See the file inst/regression_derivatives.pdf for the derivation of
# the derivatives computed by this function.
GetIVSEDerivs <- function(x, z, y, beta, w0, se_group=NULL,
                          compute_derivs=TRUE, testing=FALSE) {
  num_obs <- length(y)

  eps <- as.numeric(y - x %*% beta)

  iv_coeff <- GetIVCoefficients(x, z, y, w0)
  z_w <- iv_coeff$z_w
  zwz <- iv_coeff$zwz
  zwx_qr <- iv_coeff$zwx_qr
  betahat <- iv_coeff$betahat

  if (compute_derivs) {
    # This derivative is the same for both SE methods.
    dbetahat_dw <- solve(zwx_qr, t(z * eps))
  } else {
    dbetahat_dw <- NA
  }

  if (!is.null(se_group)) {
    ##############################################
    # Derivatives for grouped standard errors:

    # Enforces that the grouping variable is sequential and zero-indexed.
    se_group <- as.numeric(factor(se_group)) - 1

    s_mat <- GroupedSum(z * eps * w0, se_group)

    # colMeans(s_mat) is zero at the weights used for regression, but include
    # it so we can test partial derivatives.
    s_mat <- s_mat - rep(colMeans(s_mat), each=nrow(s_mat))

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


# Safely extract the IV variables from a fit object.
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



######################################################
# Functions for re-running the IV and OLS regressions.


###############
# OLS

#' Run a regression using zaminfluence code.  This should be identical
#' to ordinary regression.
#'
#' @param lm_result `r docs$lm_result`
#' @param weights `r docs$weights`
#' @param se_group `r docs$se_group`
#'
#' @return `r docs$rerun_return`
#'
#' @export
ComputeRegressionResults <- function(lm_result, weights=NULL, se_group=NULL) {
  reg_vars <- GetRegressionVariables(lm_result)
  if (is.null(weights)) {
    weights <- reg_vars$w0
  }

  # TODO: re-use the QR decomposition in GetRegressionSEDerivs
  reg_coeff <- GetRegressionCoefficients(
    x=reg_vars$x, y=reg_vars$y, w0=weights)

  reg_grad_list <- GetRegressionSEDerivs(
    x=reg_vars$x,
    y=reg_vars$y,
    beta=reg_coeff$betahat,
    w0=weights,
    se_group=se_group,
    testing=FALSE,
    compute_derivs=FALSE)
  return(list(
    betahat=reg_grad_list$betahat,
    se=reg_grad_list$se,
    se_mat=reg_grad_list$se_mat
  ))
}


#' Rerun the regression with a new subset of rows.
#'@param w_bool A boolean vector of rows to keep in the original dataframe.
#'@param lm_result `r docs$lm_result`
#'@param se_group `r docs$se_group`
#'@param save_w Optional. If `TRUE`, save the new weight vector in the output.
#'
#'@return `r docs$rerun_return`
#'@export
RerunRegression <- function(w_bool, lm_result, se_group=NULL, save_w=FALSE) {
  new_w <- rep(0.0, length(w_bool))
  if ("weights" %in% names(lm_result)) {
    new_w[w_bool] <- lm_result$weights[w_bool]
  } else {
    new_w[w_bool] <- 1.0
  }

  # Rerun using my own code; I don't want to deal with how R handles the
  # scoping of the weight variables in the regression.
  ret_list <-
    ComputeRegressionResults(lm_result, weights=new_w, se_group=se_group)
  if (save_w) {
    ret_list$w <- new_w
  }
  return(ret_list)
}


###############
# IV


#' Compute the standard error matrix for an IV regression.  Deprecated --- use
#' [ComputeIVRegressionResults()] instead.
#'
#' @param iv_res `r docs$iv_res`
#' @param se_group `r docs$se_group`
#'
#' @return `r docs$grad_return`
#'
#' @export
ComputeIVRegressionErrorCovariance <- function(iv_res, se_group=NULL) {
  warning(paste0(
    "ComputeIVRegressionErrorCovariance is deprecated; use the ",
    "se_mat output of ComputeIVRegressionResults instead."))
  iv_vars <- GetIVVariables(iv_res)
  iv_grad_list <- GetIVSEDerivs(
    x=iv_vars$x, z=iv_vars$z, y=iv_vars$y,
    beta=iv_vars$betahat, w0=iv_vars$w0, se_group=se_group)
  return(iv_grad_list$se_mat)
}


#' Run an IV regression using zaminfluence code.  This should be identical
#' to ordinary regression.
#'
#' @param lm_result `r docs$lm_result`
#' @param weights `r docs$weights`
#' @param se_group `r docs$se_group`
#'
#' @return A list containing the regression coefficients and standard errors.
#' @export
ComputeIVRegressionResults <- function(iv_res, weights=NULL, se_group=NULL) {
  iv_vars <- GetIVVariables(iv_res)
  if (is.null(weights)) {
    weights <- iv_vars$w0
  }

  # TODO: re-use the QR decomposition in GetRegressionSEDerivs
  iv_coeff <- GetIVCoefficients(
    x=iv_vars$x,
    z=iv_vars$z,
    y=iv_vars$y,
    w0=weights)

  iv_grad_list <- GetIVSEDerivs(
    x=iv_vars$x,
    z=iv_vars$z,
    y=iv_vars$y,
    beta=iv_coeff$betahat,
    w0=weights,
    se_group=se_group,
    testing=FALSE,
    compute_derivs=FALSE)

  return(list(
    betahat=iv_grad_list$betahat,
    se=iv_grad_list$se,
    se_mat=iv_grad_list$se_mat
  ))
}



#' Rerun the regression with a new subset of rows.
#'@param w_bool A boolean vector of rows to keep in the original dataframe.
#'@param iv_res `r docs$iv_res`
#'@param se_group `r docs$se_group`
#'@param save_w Optional. If TRUE, save the new weight vector in the output.
#'
#'@return `r docs$rerun_return`
#'
#'@export
RerunIVRegression <- function(w_bool, iv_res, se_group=NULL, save_w=FALSE) {
  # Rerun using my own code; I don't want to deal with how R handles the
  # scoping of the weight variables in the regression.
  if (length(iv_res$y) != length(w_bool)) {
    stop(paste0("``w_bool`` is not the same length as the regression data. ",
                "Note that re-running regression with aggregated ",
                "influence functions is not implemented."))
  }

  new_w <- rep(0.0, length(w_bool))
  if (is.null(iv_res[["weights"]])) {
    new_w[w_bool] <- 1.0
  } else {
    new_w[w_bool] <- iv_res$weights[w_bool]
  }

  ret_list <-
    ComputeIVRegressionResults(iv_res, weights=new_w, se_group=se_group)
  if (save_w) {
    ret_list$w <- new_w
  }
  return(ret_list)
}


#########################################################################
# The remaining functions wrap the actual derivative computation into
# a format that can be used by the rest of zaminfluence.  In principle,
# everything above this section could be replaced by wrappers that use
# automatic differentiation, perhaps in a different programming language.


#' Compute all influence scores for a regression.
#' @param lm_result `r docs$lm_result`
#' @param se_group `r docs$se_group`
#'
#' @return `r docs$grad_return`
#'
#' @export
ComputeRegressionInfluence <- function(lm_result, se_group=NULL) {
  reg_vars <- GetRegressionVariables(lm_result)
  reg_grad_list <- GetRegressionSEDerivs(
    x=reg_vars$x, y=reg_vars$y, beta=reg_vars$betahat,
    w0=reg_vars$w0, se_group=se_group)

  RerunFun <- function(model_fit, w_bool) {
    RerunRegression(w_bool=w_bool, lm_result=model_fit, se_group=se_group)
  }

  return(ModelGrads(model_fit=lm_result,
              n_obs=reg_vars$num_obs,
              parameter_names=colnames(lm_result$x),

              betahat=reg_vars$betahat,
              se=reg_grad_list$se,
              weights=reg_vars$w0,
              se_group=se_group,

              beta_grad=reg_grad_list$dbetahat_dw,
              se_grad=reg_grad_list$dse_dw,

              RerunFun=RerunFun))
}


#' Compute all influence scores for an IV regression.
#' @param iv_res `r docs$iv_res`
#' @param se_group `r docs$se_group`
#'
#' @return `r docs$grad_return`
#'
#' @export
ComputeIVRegressionInfluence <- function(iv_res, se_group=NULL) {
    iv_vars <- GetIVVariables(iv_res)
    iv_grad_list <- GetIVSEDerivs(
      x=iv_vars$x, z=iv_vars$z, y=iv_vars$y,
      beta=iv_vars$betahat, w0=iv_vars$w0, se_group=se_group)

    RerunFun <- function(model_fit, w_bool) {
      RerunIVRegression(w_bool=w_bool, iv_res=model_fit, se_group=se_group)
    }

    # Note that the standard errors may not match iv_res when using se_group.
    return(ModelGrads(model_fit=iv_res,
                n_obs=iv_vars$num_obs,
                parameter_names=colnames(iv_res$x$regressors),

                betahat=iv_vars$betahat,
                se=iv_grad_list$se,
                weights=iv_vars$w0,
                se_group=se_group,

                beta_grad=iv_grad_list$dbetahat_dw,
                se_grad=iv_grad_list$dse_dw,

                RerunFun=RerunFun
              ))
}


#' Compute the influence functions for all regressors given a model fit.
#' @param model_fit `r docs$model_fit`
#' @param se_group `r docs$se_group`
#'
#' @return `r docs$model_grads`
#'
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
