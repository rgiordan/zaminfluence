# A list of strings for arguments and return values that repeatedly recur
# in the docs.

docs <- list(
  model_grads=paste0(
      "A model gradient object (e.g. as produced by [ComputeModelInfluence])"
  ),
  param_infl=paste0(
    "A parameter influence object (e.g. as produced by ",
    "[AppendTargetRegressorInfluence])."
  ),
  qoi=paste0(
    "A quantity of interest object (e.g. as produced by ",
    "[QOIInfluence])."
  ),
  signal=paste0("A signal object (e.g. one of the elements of the list ",
    "produced by [GetInferenceSignalsForParameter])."),
  apip=paste0("An approximation perturbation inducing proportion object ",
"(e.g. as produced by [GetAPIPForQOI])"),
drop_inds="The indices to drop (in the order of the original data)",

  model_fit="The fit from [lm()] or [AER::ivreg()].",
  lm_result="The regression result, i.e, the output of [lm()].",
  iv_res="The iv regression result, i.e, the output of [ivreg()].",
  se_group="Optional. The standard error grouping variable.",
  rerun_return=paste0(
    "A list containing the new regression estimate, standard error ",
    "covariance (`se_mat`), and standard errors (`se`).  The standard errors ",
    "are the square root of the diagonal of `se_mat`, and can be interpreted ",
    "as an estimate of the standard deviation of the coefficient estimates."
  ),
  grad_return=paste0(
    "A list containing the sensitivity of the coefficients and standard ",
    "errors to data re-weighting.  The list entries are\n",
    "\\describe{ \n",
    "  \\item{model_fit}{The original fit as returned by the fitting function} \n",
    "  \\item{num_obs}{The number of observations given weights} \n",
    "  \\item{parameter_names}{The names of the regressors as strings} \n",
    "  \\item{grad_fun}{The function used to compute the derivatives} \n",
    "  \\item{param}{The estimated regression coefficient ",
    "at the original weights} \n",
    "  \\item{se}{The standard errors at the original weights} \n",
    "  \\item{weights}{The original weights} \n",
    "  \\item{param_grad}{The gradient of param with respect to the weights} \n",
    "  \\item{se_grad}{The gradient of se with respect to the weights} \n",
    "}\n"
  ),
  weights="Optional.  A vector of weights.  If unset, use the original weights.",
  sig_num_ses="How wide the confidence interval is, as a number of standard errors."
)

# Usage example:
##' @param lm_result `r docs$lm_result`
