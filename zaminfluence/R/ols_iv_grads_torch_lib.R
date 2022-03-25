######################################################3
# Ordinary least squares
library(torch)


ValidateKeepInds <- function(keep_inds, x) {
  # If keep_inds is unspecified, compute derivatives for all coefficients.
  if (is.null(keep_inds)) {
    keep_inds <- 1:ncol(x)
  }
  keep_inds <- as.integer(keep_inds)
  if (min(keep_inds) < 1) {
    stop("keep_inds must be >= 1")
  }
  if (max(keep_inds) > ncol(x)) {
    stop("keep_inds must be no greater than the number of x columns")
  }
  if (length(unique(keep_inds)) != length(keep_inds)) {
    stop("keep_inds must not contain repeats")
  }
  return(keep_inds)
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
  parameter_names <- colnames(lm_res$x)

  return(list(x=x, y=y, num_obs=num_obs, w0=w0,
              z_equals_x=TRUE,
              betahat=betahat, parameter_names=parameter_names))
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
  parameter_names <- colnames(x)
  return(list(x=x, y=y, z=z, num_obs=num_obs, w0=w0,
              z_equals_x=FALSE,
              betahat=betahat, parameter_names=parameter_names))
}



######################################################
# Torch stuff.


TorchGroupedAggregate <- function(src_mat, inds) {
    # https://discuss.pytorch.org/t/groupby-aggregate-mean-in-pytorch/45335
    max_ind <- inds$max() %>% as.integer()
    zero_mat <- torch_tensor(
        matrix(0, nrow=max_ind, ncol=src_mat$shape[2]),
        dtype=src_mat$dtype)
    inds_rep <- inds[["repeat"]](list(1, src_mat$shape[2]))
    agg_mat <- zero_mat$scatter_add(dim=1, index=inds_rep, src=src_mat)
    return(agg_mat)
}

# inds <- torch_tensor(c(1, 1, 3, 2, 2) %>% matrix(5, 1), dtype=torch_int64())
# src <- torch_tensor(c(rep(0.1, 5), rep(0.2, 5)) %>% matrix(5, 2))
# TorchGroupedAggregate(src, inds)


DefineTorchVars <- function(iv_vars) {
    if (iv_vars$z_equals_x) {
      stopifnot(all(c("x", "y", "w0") %in% names(iv_vars)))
    } else {
      stopifnot(all(c("x", "y", "z", "w0") %in% names(iv_vars)))
    }

    num_obs <- nrow(iv_vars$x)
    num_cols <- ncol(iv_vars$x)

    stopifnot(length(iv_vars$y) == num_obs)
    stopifnot(nrow(iv_vars$x) == num_obs)
    stopifnot(ncol(iv_vars$x) == num_cols)

    if (!iv_vars$z_equals_x) {
      stopifnot(nrow(iv_vars$z) == num_obs)
      stopifnot(ncol(iv_vars$z) == num_cols)
    }

    x <- torch_tensor(iv_vars$x, requires_grad=FALSE, dtype=torch_double())
    if (iv_vars$z_equals_x) {
      z <- x
    } else {
      z <- torch_tensor(iv_vars$z, requires_grad=FALSE, dtype=torch_double())
    }
    y <- torch_tensor(matrix(iv_vars$y, ncol=1),
                      requires_grad=FALSE, dtype=torch_double())
    w <- torch_tensor(matrix(iv_vars$w0, ncol=1),
                      requires_grad=TRUE, dtype=torch_double())

    z_w <- z * w
    z_w_t <- z_w$transpose(2, 1)
    zwx <- torch_matmul(z_w_t, x)
    #t(iv_vars$z) %*% iv_vars$z - zwz # Why 1e-6?
    betahat <- torch::linalg_solve(zwx, torch_matmul(z_w_t, y))
    eps <- y - torch_matmul(x, betahat)

    return(list(
        num_obs=num_obs,
        num_cols=num_cols,
        x=x,
        z=z,
        y=y,
        w=w,
        zwx=zwx,
        zw_t=z_w_t,
        eps=eps,
        betahat=betahat
    ))
}


#'@export
GetIVRegressionSEDerivsTorch <- function(
      iv_vars, se_group=NULL, keep_inds=NULL, compute_derivs=TRUE) {

    keep_inds <- ValidateKeepInds(keep_inds, iv_vars$x)
    tv <- DefineTorchVars(iv_vars)

    if (!is.null(se_group)) {
        tv$score_mat <- tv$z * tv$eps * tv$w
        tv$se_group <-
            torch_tensor(
                as.integer(factor(se_group)) %>%
                matrix(ncol=1),
                dtype=torch_int64())
        tv$score_sum <- TorchGroupedAggregate(
          src_mat=tv$score_mat, inds=tv$se_group)
        tv$s_mat <- tv$score_sum - tv$score_sum$mean(dim=1, keepdim=TRUE)

        num_groups <- tv$s_mat$shape[1]
        tv$v_mat <- torch_matmul(tv$s_mat$transpose(2, 1), tv$s_mat) / num_groups

        tv$zwx_inv_vmat <- linalg_solve(tv$zwx, tv$v_mat)
        tv$se_cov_mat <- linalg_solve(
          tv$zwx, torch_transpose(tv$zwx_inv_vmat, 2, 1)) * num_groups
    } else {
        tv$sig2_hat <- torch_sum(
          tv$w * (tv$eps ** 2)) / (tv$num_obs - tv$num_cols)

        if (iv_vars$z_equals_x) {
            # Regression is when Z == X and we can save some computation
            tv$se_cov_mat <- tv$sig2_hat * torch_inverse(tv$zwx)
        } else {
            # IV is when Z != X
            tv$zwz <- torch_matmul(tv$zw_t, tv$z)
            tv$zwx_inv_zwz <- linalg_solve(tv$zwx, tv$zwz)
            tv$se_cov_mat <- tv$sig2_hat * linalg_solve(
              tv$zwx, tv$zwx_inv_zwz$transpose(2, 1))
        }
    }
    tv$betahat_se <- torch_sqrt(torch_diag(tv$se_cov_mat))

    return_list <- list(
        tv=tv,
        betahat=as.numeric(tv$betahat),
        betahat_se=as.numeric(tv$betahat_se))

    if (compute_derivs) {
      # betahat
      betahat_infl_mat <- matrix(NA, nrow=length(keep_inds), ncol=tv$num_obs)
      for (di in 1:length(keep_inds)) {
          d <- keep_inds[di]
          betahat_infl_mat[di, ] <-
              torch::autograd_grad(
                  tv$betahat[d], tv$w, retain_graph=TRUE)[[1]] %>% as.numeric()
      }

      # standard errors
      se_infl_mat <- matrix(NA, nrow=length(keep_inds), ncol=tv$num_obs)
      for (di in 1:length(keep_inds)) {
          d <- keep_inds[di]
          se_infl_mat[di, ] <-
              torch::autograd_grad(
                  tv$betahat_se[d],
                  tv$w, retain_graph=TRUE)[[1]] %>% as.numeric()
      }

      return_list$betahat_infl_mat <- betahat_infl_mat
      return_list$betahat_se_infl_mat <- se_infl_mat
    }

    return(return_list)
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
  if (!is.null(weights)) {
    reg_vars$w0 <- weights
  }

  reg_grad_list <- GetIVRegressionSEDerivsTorch(
    iv_vars=reg_vars,
    se_group=se_group,
    keep_inds=NULL,
    compute_derivs=FALSE)

  return(list(
    betahat=reg_grad_list$betahat,
    se=reg_grad_list$betahat_se,
    se_mat=reg_grad_list$tv$se_cov_mat
  ))
}




###############
# IV


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
  if (!is.null(weights)) {
    iv_vars$w0 <- weights
  }

  reg_grad_list <- GetIVRegressionSEDerivsTorch(
    iv_vars=iv_vars,
    se_group=se_group,
    keep_inds=NULL,
    compute_derivs=FALSE)

  return(list(
    betahat=reg_grad_list$betahat,
    se=reg_grad_list$betahat_se,
    se_mat=reg_grad_list$tv$se_cov_mat
  ))
}





#########################################################################
# The remaining functions wrap the actual derivative computation into
# a format that can be used by the rest of zaminfluence.  In principle,
# everything above this section could be replaced by wrappers that use
# automatic differentiation, perhaps in a different programming language.


GetKeepInds <- function(coeff_names, keep_pars=NULL) {
    if (is.null(keep_pars)) {
      return(1:length(coeff_names))
    }
    missing_pars <- setdiff(keep_pars, coeff_names)
    if (length(missing_pars) > 0) {
        stop(sprintf("Parameters %s are not present in model", paste(missing_pars, collapse=", ")))
    }
    inds <- setNames(1:length(coeff_names), coeff_names)
    return(inds[keep_pars] %>% unname())
}


#' Compute all influence scores for a regression.
#' @param lm_result `r docs$lm_result`
#' @param se_group `r docs$se_group`
#'
#' @return `r docs$grad_return`
#'
#' @export
ComputeRegressionInfluence <- function(
    lm_result, se_group=NULL, keep_pars=NULL) {

  all_par_names <- names(coefficients(lm_result))
  if (is.null(keep_pars)) {
    keep_pars <- all_par_names
  }
  keep_inds <- GetKeepInds(all_par_names, keep_pars)

  reg_vars <- GetRegressionVariables(lm_result)
  # reg_grad_list <- GetRegressionSEDerivs(
  #   x=reg_vars$x, y=reg_vars$y, beta=reg_vars$betahat,
  #   w0=reg_vars$w0, se_group=se_group)
  reg_grad_list <- GetIVRegressionSEDerivsTorch(
    iv_vars=reg_vars,
    se_group=se_group,
    keep_inds=keep_inds,
    compute_derivs=TRUE)

  RerunFun <- function(weights) {
    ret_list <-
      ComputeRegressionResults(lm_result, weights=weights, se_group=se_group)
    return(ModelFit(
      fit_object=ret_list,
      num_obs=reg_vars$num_obs,
      param=ret_list$betahat,
      se=ret_list$se,
      parameter_names=reg_vars$parameter_names,
      weights=weights,
      se_group=se_group))
  }

  model_fit <- ModelFit(
    fit_object=lm_result,
    num_obs=reg_vars$num_obs,
    parameter_names=reg_vars$parameter_names,
    param=reg_vars$betahat,
    se=reg_grad_list$betahat_se,
    weights=reg_vars$w0,
    se_group=se_group)

  rownames(reg_grad_list$betahat_infl_mat) <- keep_pars
  rownames(reg_grad_list$betahat_se_infl_mat) <- keep_pars

  return(ModelGrads(model_fit=model_fit,
                    param_grad=reg_grad_list$betahat_infl_mat,
                    se_grad=reg_grad_list$betahat_se_infl_mat,
                    RerunFun=RerunFun))

}


#' Compute all influence scores for an IV regression.
#' @param iv_res `r docs$iv_res`
#' @param se_group `r docs$se_group`
#'
#' @return `r docs$grad_return`
#'
#' @export
ComputeIVRegressionInfluence <- function(
      iv_res, se_group=NULL, keep_pars=NULL) {

    all_par_names <- names(coefficients(iv_res))
    if (is.null(keep_pars)) {
      keep_pars <- all_par_names
    }
    keep_inds <- GetKeepInds(all_par_names, keep_pars)

    iv_vars <- GetIVVariables(iv_res)
    # iv_grad_list <- GetIVSEDerivs(
    #   x=iv_vars$x, z=iv_vars$z, y=iv_vars$y,
    #   beta=iv_vars$betahat, w0=iv_vars$w0, se_group=se_group)
    iv_grad_list <- GetIVRegressionSEDerivsTorch(
      iv_vars=iv_vars,
      se_group=se_group,
      keep_inds=keep_inds,
      compute_derivs=TRUE)

    RerunFun <- function(weights) {
        ret_list <-
          ComputeIVRegressionResults(iv_res, weights=weights, se_group=se_group)
        return(ModelFit(
          fit_object=ret_list,
          num_obs=iv_vars$num_obs,
          param=ret_list$betahat,
          se=ret_list$se,
          parameter_names=iv_vars$parameter_names,
          weights=weights,
          se_group=se_group))
      }

    model_fit <- ModelFit(
      fit_object=iv_res,
      num_obs=iv_vars$num_obs,
      parameter_names=iv_vars$parameter_names,
      param=iv_vars$betahat,
      se=iv_grad_list$betahat_se,
      weights=iv_vars$w0,
      se_group=se_group)

    rownames(iv_grad_list$betahat_infl_mat) <- keep_pars
    rownames(iv_grad_list$betahat_se_infl_mat) <- keep_pars

    # Note that the standard errors may not match iv_res when using se_group.
    return(ModelGrads(model_fit=model_fit,
                      param_grad=iv_grad_list$betahat_infl_mat,
                      se_grad=iv_grad_list$betahat_se_infl_mat,
                      RerunFun=RerunFun))
}


#' Compute the influence functions for all regressors given a model fit.
#' @param model_fit `r docs$model_fit`
#' @param se_group `r docs$se_group`
#'
#' @return `r docs$model_grads`
#'
#' @export
ComputeModelInfluence <- function(fit_object, se_group=NULL, keep_pars=NULL) {
  valid_classes <- c("lm", "ivreg")
  model_class <- class(fit_object)
  if (!(model_class %in% valid_classes)) {
    stop(sprintf("The class of `model_fit` must be one of %s",
                 paste(valid_classes, collapse=", ")))
  }
  if (model_class == "lm") {
    return(ComputeRegressionInfluence(
      fit_object, se_group=se_group, keep_pars=keep_pars))
  } else if (model_class == "ivreg") {
    return(ComputeIVRegressionInfluence(
      fit_object, se_group=se_group, keep_pars=keep_pars))
  } else {
    # Redundant, so sue me.
    stop(sprint("Unknown model class %s", model_class))
  }
}
