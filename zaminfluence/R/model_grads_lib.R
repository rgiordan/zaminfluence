# Define S3 class for ModelGrads

new_ModelGrads <- function(
    model_fit,
    n_obs,
    parameter_names,
    betahat,
    se,
    beta_grad,
    se_grad,
    RerunFun,
    weights,
    se_group) {
  return(structure(
    list(model_fit=model_fit,
                n_obs=n_obs,
                parameter_names=as.character(parameter_names),

                betahat=betahat,
                se=se,
                weights=weights,
                se_group=se_group,

                beta_grad=beta_grad,
                se_grad=se_grad,

                RerunFun=RerunFun),
    class="ModelGrads"
  ))
}


validate_ModelGrads <- function(model_grads) {
  stopifnot(class(model_grads) == "ModelGrads")
  stopifnot(length(model_grads$weights) == model_grads$n_obs)
  if (!is.null(model_grads$se_group)) {
    stopifnot(length(model_grads$se_group) == model_grads$n_obs)
  }

  CheckGradDim <- function(grad_mat) {
    stopifnot(length(dim(grad_mat)) == 2)
    stopifnot(ncol(grad_mat) == model_grads$n_obs)
    stopifnot(nrow(grad_mat) == length(model_grads$parameter_names))
  }

  CheckGradDim(model_grads$beta_grad)
  CheckGradDim(model_grads$se_grad)

  return(invisible(model_grads))
}


#'@export
ModelGrads <- function(
    model_fit,
    n_obs,
    parameter_names,
    betahat,
    se,
    beta_grad,
    se_grad,
    RerunFun,
    weights=NULL,
    se_group=NULL) {

  if (is.null(weights)) {
    weights <- rep(1.0, n_obs)
  }

  return(validate_ModelGrads(new_ModelGrads(
      model_fit=model_fit,
      n_obs=n_obs,
      parameter_names=parameter_names,
      betahat=betahat,
      se=se,
      beta_grad=beta_grad,
      se_grad=se_grad,
      RerunFun=RerunFun,
      weights=weights,
      se_group=se_group)))
}
