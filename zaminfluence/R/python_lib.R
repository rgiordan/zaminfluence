library(reticulate)
library(broom)


# The role of this library is to calculate regression and sensitivity estimates
# for a given dataset.  The two main functions are RegressWithInfluence,
# which runs the regression and calculates influence vectors, and
# GetTargetRegressorGrads, which processes the output of RegressWithInfluence
# for a particular target regressor and confidence interval width.
#
# I will call the output of GetTargetRegressorGrads, which is a dataframe
# with certain columns and attributes, a "gradient dataframe".  The gradient
# dataframe is further processed in sorting_lib.R to calculate adversarial
# removal sets.


#' Start a Python session using the given Python 3 binary.
#'
#' @param python_path The path to a Python 3 binary.
#' @return The `reticulate::py_main` object.
#' @export
InitializePython <- function(python_path="/usr/bin/python3") {
    if (!file.exists(python_path)) {
        stop(sprintf("Python path %s missing.", python_path))
    }

    reticulate::use_python(python_path)
    reticulate::py_run_string("
import autograd
import autograd.numpy as np
import scipy as sp
from copy import deepcopy
import paragami
import vittles
import regsens_rgiordandev
from regsens_rgiordandev import iv_sensitivity_lib as iv_lib
")
    return(reticulate::import_main())
}


#######################
# Ordinary regression

GetLMWeights <- function(lm_res) {
  n_obs <- nrow(lm_res$x)
  if (is.null(lm_res$weights)) {
    weights <- rep(1.0, n_obs)
  } else {
    weights <- lm_res$weights
  }
  return(weights)
}


#' Set python variables corresponding to the output of lm.
SetPythonRegressionVariables <- function(lm_res, se_group=NULL) {
  if (!(("x" %in% names(lm_res)) &
        ("y" %in% names(lm_res)))) {
    stop("You must run lm with the arguments x=TRUE and y=TRUE.")
  }
  n_obs <- nrow(lm_res$x)
  weights <- GetLMWeights(lm_res)
  py_main <- reticulate::import_main()
  py_main$w0 <- as.array(as.numeric(weights))
  py_main$x <- as.array(as.matrix(lm_res$x))
  py_main$y <- as.array(as.numeric(lm_res$y))
  py_main$beta <- as.array(as.numeric(lm_res$coefficients))
  if (is.null(se_group)) {
      reticulate::py_run_string("se_group = None")
  } else {
      py_main$se_group <- as.array(as.integer(se_group))
  }
  return(py_main)
}


#' Compute all influence scores for a regression.
#' @param lm_result The output of a call to lm.
#' @param se_group Optional, a vector of integers defining a standard error grouping.
#' @return A list containing the regression and influence result.
#' @export
ComputeRegressionInfluence <- function(lm_result, se_group=NULL) {
  py_main <- SetPythonRegressionVariables(lm_result, se_group=se_group)
  reg <- broom::tidy(lm_result)
  reticulate::py_run_string("
betahat = regsens_rgiordandev.reg(y, x, w=w0)
se, betahat_grad, se_grad = regsens_rgiordandev.get_regression_w_grads(
    betahat, y, x, w0, se_group=se_group)
")
  if (max(abs(py_main$betahat - reg$estimate)) > 1e-8) {
      warning("Regression coefficients do not match.")
  }

  # Note that the standard errors may not match lm_result when using se_group.
  return(list(model_fit=lm_result,
              n_obs=nrow(lm_result$x),
              regressor_names=colnames(lm_result$x),
              grad_fun="get_regression_w_grads",

              betahat=py_main$betahat,
              se=py_main$se,
              weights=py_main$w0,

              beta_grad=py_main$betahat_grad,
              se_grad=py_main$se_grad)
          )
}

#################################
# Regression moment sensitivity

#' Compute the sensitivity of a regression to the moment condition on the
#' residuals.  That is, compute the gradients with respect to `offset`,
#' where the regression is fit with the moment condition Cov(x, resid) = offset.
#' @param lm_result The output of a call to lm.
#' @param se_group Optional, a vector of integers defining a standard error grouping.
#' @return A list containing the regression and influence result.
#' @export
ComputeRegressionMomentSensitivity <- function(lm_result, se_group=NULL) {
  weights <- GetLMWeights(lm_result)
  if (max(abs(weights - 1)) > 1e-8) {
    stop("ComputeRegressionMomentSensitivity only works for identity weights.")
  }
  py_main <- SetPythonRegressionVariables(lm_result, se_group=se_group)

  # Set the base offset, which must be zero, since lm does not allow for
  # a more flexible offset specification.
  py_main <- reticulate::import_main()
  x_dim <- ncol(lm_result$x)
  offset0 <- rep(0.0, x_dim)
  py_main$offset0 <- as.array(as.numeric(offset0))

  reg <- broom::tidy(lm_result)
  reticulate::py_run_string("
betahat = regsens_rgiordandev.reg(y, x, w=w0)
se, betahat_grad, se_grad = regsens_rgiordandev.get_regression_offset_grads(
    betahat, y, x, offset0, se_group=se_group)
")
  if (max(abs(py_main$betahat - reg$estimate)) > 1e-8) {
      warning("Regression coefficients do not match.")
  }

  # Note that the standard errors may not match lm_result when using se_group.
  return(list(model_fit=lm_result,
              n_obs=nrow(lm_result$x),
              regressor_names=colnames(lm_result$x),
              grad_fun="get_regression_offset_grads",

              betahat=py_main$betahat,
              se=py_main$se,
              weights=py_main$w0,
              offset=py_main$offset0,

              beta_grad=py_main$betahat_grad,
              se_grad=py_main$se_grad)
          )
}



#' Compute a regression using the moment condition Cov(x, resid) = offset.
#' @param lm_result The output of a call to lm.
#' @param offset The target moment condition.
#' @param se_group Optional, a vector of integers defining a standard error grouping.
#' @return A list containing the regression and influence result.
#' @export
RegressWithOffset <- function(lm_result, offset, se_group=NULL) {
  weights <- GetLMWeights(lm_result)
  if (max(abs(weights - 1)) > 1e-8) {
    stop("ComputeRegressionMomentSensitivity only works for identity weights.")
  }
  py_main <- SetPythonRegressionVariables(lm_result, se_group=se_group)

  # Set the base offset, which must be zero, since lm does not allow for
  # a more flexible offset specification.
  py_main <- reticulate::import_main()
  x_dim <- ncol(lm_result$x)
  if (length(offset) != x_dim) {
    stop("`offset` must be a vector of the same length as the dimension of x.")
  }
  py_main$offset0 <- as.array(as.numeric(offset))

  reg <- broom::tidy(lm_result)
  reticulate::py_run_string("
betahat = regsens_rgiordandev.reg(y, x, offset=offset0)
# Offsets do not affect the standard errors.
se_cov = regsens_rgiordandev.get_standard_error_matrix(
    betahat, y, x, w=np.ones(x.shape[0]), se_group=se_group)
")

  return(list(betahat=py_main$betahat,
              se_cov=py_main$se_cov,
              offset=py_main$offset0)
          )
}


#######################
# IV regression

#' Set python variables corresponding to the output of AER::ivreg.
#' @export
SetPythonIVRegressionVariables <- function(iv_res, se_group=NULL) {
    if (!(("x" %in% names(iv_res)) &
          ("y" %in% names(iv_res)))) {
        stop("You must run ivreg with the arguments x=TRUE and y=TRUE.")
    }
    x <- iv_res$x$regressors
    z <- iv_res$x$instruments
    n_obs <- nrow(x)
    if (is.null(iv_res$weights)) {
        weights <- rep(1.0, n_obs)
    } else {
        weights <- iv_res$weights
    }
    py_main <- reticulate::import_main()
    py_main$w0 <- as.array(as.numeric(weights))
    py_main$x <- as.array(as.matrix(x))
    py_main$z <- as.array(as.matrix(z))
    py_main$y <- as.array(as.numeric(iv_res$y))
    py_main$beta <- as.array(as.numeric(iv_res$coefficients))
    if (is.null(se_group)) {
        reticulate::py_run_string("se_group = None")
    } else {
        py_main$se_group <- as.array(as.integer(se_group))
    }
    return(py_main)
}


#' Compute all influence scores for an IV regression.
#' @param df A dataframe with regression data.
#' @param iv_res The output of an IV regression computed with AER::ivreg.
#' @param se_group Optional, a vector of integers defining a standard error grouping.
#' @return A list containing the regression and influence result..
#' @export
ComputeIVRegressionInfluence <- function(iv_res, se_group=NULL) {
    py_main <- SetPythonIVRegressionVariables(iv_res, se_group=se_group)
    reg <- broom::tidy(iv_res)
    reticulate::py_run_string("
betahat = iv_lib.iv_reg(y, x, z, w=w0)
se, betahat_grad, se_grad = iv_lib.get_iv_regression_w_grads(
    betahat, y, x, z, w0, se_group=se_group)
")
    if (max(abs(py_main$betahat - reg$estimate)) > 1e-8) {
        warning("Regression coefficients do not match.")
    }

    # Note that the standard errors may not match iv_res when using se_group.
    return(list(model_fit=iv_res,
                n_obs=length(iv_res$y),
                regressor_names=colnames(iv_res$x$regressors),
                grad_fun="get_iv_regression_w_grads",

                betahat=py_main$betahat,
                se=py_main$se,
                weights=py_main$w0,

                beta_grad=py_main$betahat_grad,
                se_grad=py_main$se_grad)
    )
}


#' @export
ComputeIVRegressionErrorCovariance <- function(iv_res, se_group=NULL) {
  py_main <- SetPythonIVRegressionVariables(iv_res, se_group=se_group)
  reticulate::py_run_string("
betahat = iv_lib.iv_reg(y, x, z, w=w0)
se2 = iv_lib.get_iv_standard_error_matrix(betahat, y, x, z, w0, se_group=se_group)
")
  return(py_main$se2)
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
