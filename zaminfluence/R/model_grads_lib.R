
new_ModelFit <- function(
  fit_object, parameter_names, betahat, se, weights, se_group) {
  return(structure(
    list(fit_object=fit_object,
         n_obs=n_obs,
         parameter_names=as.character(parameter_names),
         betahat=betahat,
         se=se,
         weights=weights,
         parameter_dim=length(betahat),
         se_group=se_group),
    class="ModelFit"
    ))
}


validate_ModelFit <- function(model_fit) {
  stopifnot(class(model_fit) == "ModelFit")
  StopIfNotNumericScalar(model_fit$n_obs)
  StopIfNotNumericScalar(model_fit$parameter_dim)

  n_obs <- model_fit$n_obs
  stopifnot(length(model_fit$weights) == n_obs)
  if (!is.null(model_fit$se_group)) {
    stopifnot(length(model_fit$se_group) == n_obs)
  }

  parameter_dim <- model_fit$parameter_dim
  stopifnot(length(parameter_names) == parameter_dim)
  stopifnot(length(betahat) == parameter_dim)
  stopifnot(length(se) == parameter_dim)
}


#'@export
ModelFit <- function(fit_object, betahat, se,
                     parameter_names=NULL, weights=NULL, se_group=NULL) {
    if (is.null(weights)) {
        weights <- rep(1.0, n_obs)
    }
    if (is.null(parameter_names)) {
        parameter_names <- sprintf("theta%d", 1:length(betahat))
    }
    return(validate_ModelFit(new_ModelFit(
      fit_object=fit_object,
      parameter_names=parameter_names,
      betahat=betahat,
      se=se,
      weights=weights,
      se_group=se_group
    )))
)


# Define S3 class for ModelGrads

new_ModelGrads <- function(
    model_fit,
    beta_grad,
    se_grad,
    RerunFun) {
  return(structure(
    list(model_fit=model_fit,

         beta_grad=beta_grad,
         se_grad=se_grad,

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
    stopifnot(ncol(grad_mat) == model_fit$n_obs)
    stopifnot(nrow(grad_mat) == model_fit$parameter_dim)
  }

  CheckGradDim(model_grads$beta_grad)
  CheckGradDim(model_grads$se_grad)

  return(invisible(model_grads))
}


#'@export
ModelGrads <- function(
    model_fit,
    beta_grad,
    se_grad,
    RerunFun) {

  return(validate_ModelGrads(new_ModelGrads(
      model_fit=model_fit,
      beta_grad=beta_grad,
      se_grad=se_grad,
      RerunFun=RerunFun)))
}
