
#############################################################
# Functions to help re-running to check the approximations.

#######################################
# Functions for ordinary regression


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


#######################################
# Functions for IV regression




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
