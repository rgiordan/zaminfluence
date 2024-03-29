
new_ModelFit <- function(
  fit_object, num_obs, parameter_names, param, se, weights, se_group) {
  return(structure(
    list(fit_object=fit_object,
         num_obs=num_obs,
         parameter_names=as.character(parameter_names),
         param=param,
         se=se,
         weights=weights,
         parameter_dim=length(param),
         se_group=se_group),
    class="ModelFit"
    ))
}


validate_ModelFit <- function(model_fit) {
  stopifnot(class(model_fit) == "ModelFit")
  StopIfNotNumericScalar(model_fit$num_obs)
  StopIfNotNumericScalar(model_fit$parameter_dim)

  num_obs <- model_fit$num_obs
  stopifnot(length(model_fit$weights) == num_obs)
  if (!is.null(model_fit$se_group)) {
    stopifnot(length(model_fit$se_group) == num_obs)
  }

  parameter_dim <- model_fit$parameter_dim
  stopifnot(length(model_fit$parameter_names) == parameter_dim)
  stopifnot(length(model_fit$param) == parameter_dim)
  stopifnot(length(model_fit$se) == parameter_dim)
  return(invisible(model_fit))
}



#'@export
ModelFit <- function(fit_object, num_obs, param, se,
                     parameter_names=NULL, weights=NULL, se_group=NULL) {
    if (is.null(weights)) {
        weights <- rep(1.0, num_obs)
    }
    if (is.null(parameter_names)) {
        parameter_names <- sprintf("theta%d", 1:length(param))
    }
    return(validate_ModelFit(new_ModelFit(
      fit_object=fit_object,
      num_obs=num_obs,
      parameter_names=parameter_names,
      param=param,
      se=se,
      weights=weights,
      se_group=se_group
    )))
}


# Define S3 class for ModelGrads

new_ModelGrads <- function(
    model_fit,
    parameter_names,
    param_grad,
    se_grad,
    param_infls,
    RerunFun,
    PredictFun) {
  return(structure(
    list(model_fit=model_fit,

         parameter_names=parameter_names,
         param_grad=param_grad,
         se_grad=se_grad,

         param_infls=param_infls,

         RerunFun=RerunFun),
    class="ModelGrads"
  ))
}


validate_ModelGrads <- function(model_grads) {
  stopifnot(class(model_grads) == "ModelGrads")
  validate_ModelFit(model_grads$model_fit)
  model_fit <- model_grads$model_fit

  grad_pars <- model_grads$parameter_names
  stopifnot(all(grad_pars %in% model_fit$parameter_names))

  CheckGradDim <- function(grad_mat) {
    stopifnot(length(dim(grad_mat)) == 2)
    stopifnot(ncol(grad_mat) == model_fit$num_obs)
    stopifnot(nrow(grad_mat) == length(grad_pars))
  }

  CheckGradDim(model_grads$param_grad)
  CheckGradDim(model_grads$se_grad)

  stopifnot(class(model_grads$param_infls) == "list")
  stopifnot(all(names(model_grads$param_infls) %in% grad_pars))
  for (param_infl in model_grads$param_infls) {
    stopifnot(class(param_infl) == "ParameterInferenceInfluence")
  }

  return(invisible(model_grads))
}


#'@export
PredictModelFit <- function(model_grads, weights) {
    stopifnot(class(model_grads) == "ModelGrads")
    stopifnot(is.numeric(weights))

    model_fit <- model_grads$model_fit
    stopifnot(length(weights) == model_fit$num_obs)

    weight_diff <- weights - model_fit$weights
    kept_indices <- GetParameterIndex(model_fit, model_grads$parameter_names)
    param_pred <- model_fit$param[kept_indices] + model_grads$param_grad %*% weight_diff
    se_pred <- model_fit$se[kept_indices] + model_grads$se_grad %*% weight_diff
    pred_fit <-
        ModelFit(
            fit_object="prediction",
            num_obs=model_grads$model_fit$num_obs,
            param=param_pred,
            se=se_pred,
            parameter_names=model_grads$parameter_names,
            weights=weights,
            se_group=model_fit$se_group)
    return(pred_fit)
}



#'@export
ModelGrads <- function(
    model_fit,
    param_grad,
    se_grad,
    RerunFun) {

  parameter_names <- rownames(param_grad)
  if (any(rownames(se_grad) != parameter_names)) {
    stop(paste0(
      "The rownames of se_grad and param_grad must match.",
      "rownames(param_grad) = (",
      paste(rownames(param_grad), collapse=","), ")\n",
      "rownames(se_grad) = (",
      paste(rownames(se_grad), collapse=","), ")\n",
    ))
  }
  return(validate_ModelGrads(new_ModelGrads(
      model_fit=model_fit,
      parameter_names=parameter_names,
      param_grad=param_grad,
      se_grad=se_grad,
      RerunFun=RerunFun,
      param_infls=list())))
}



#'@export
GetParameterIndex <- function(m, par_name) {
  UseMethod("GetParameterIndex")
}


#'@export
GetParameterIndex.ModelGrads <- function(model_grads, par_names) {
  return(GetParameterIndexLocal(
    model_grads$parameter_names, par_names, object_class="ModelGrads"))
}


#'@export
GetParameterIndex.ModelFit <- function(model_fit, par_names) {
  return(GetParameterIndexLocal(
    model_fit$parameter_names, par_names, object_class="ModelFit"))
}


GetParameterIndexLocal <- function(all_par_names, par_names, object_class) {
  missing_names <- setdiff(par_names, all_par_names)
  if (length(missing_names) > 0) {
    stop(paste0(
      "Parameter names ",
      paste(missing_names, collapse=", "), " not found in ",
      object_class, "."))
  }
  target_indices <- setNames(1:length(all_par_names), all_par_names)[par_names]
  return(target_indices)
}
