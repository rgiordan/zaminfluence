
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
    param_grad,
    se_grad,
    param_infls,
    RerunFun,
    PredictFun) {
  return(structure(
    list(model_fit=model_fit,

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

  CheckGradDim <- function(grad_mat) {
    stopifnot(length(dim(grad_mat)) == 2)
    stopifnot(ncol(grad_mat) == model_fit$num_obs)
    stopifnot(nrow(grad_mat) == model_fit$parameter_dim)
  }

  CheckGradDim(model_grads$param_grad)
  CheckGradDim(model_grads$se_grad)

  stopifnot(class(model_grads$param_infls) == "list")
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
    param_pred <- model_fit$param + model_grads$param_grad %*% weight_diff
    se_pred <- model_fit$se + model_grads$se_grad %*% weight_diff
    pred_fit <-
        ModelFit(
            fit_object="prediction",
            num_obs=model_grads$model_fit$num_obs,
            param=param_pred,
            se=se_pred,
            parameter_names=model_fit$parameter_names,
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

  return(validate_ModelGrads(new_ModelGrads(
      model_fit=model_fit,
      param_grad=param_grad,
      se_grad=se_grad,
      RerunFun=RerunFun,
      param_infls=list())))
}
