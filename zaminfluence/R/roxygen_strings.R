# A list of strings for arguments and return values that repeatedly recur
# in the docs.

docs <- list(
  infl_df="A single influence dataframe, dude!",
  alpha_colname=paste0(
      "A string with the name of the target alpha column (typically ",
      "`prop_removed` or `num_removed`)"
  ),
  alpha_val="The alpha value to target corresponding to `alpha_colname`.",
  influence_dfs="Influence dataframes, dude!",
  lm_result="The regression result, i.e, the output of [lm()].",
  iv_res="The iv regression result, i.e, the output of [ivreg()].",
  se_group="Optional. The standard error grouping variable.",
  rerun_return=paste0(
    "A list containing the new regression estimate, standard error ",
    "covariance (`se_mat`), and standard errors (`se`).  The standard errors ",
    "are the square root of the diagonal of `se_mat`, and can be interpreted ",
    "as an estimate of the standard deviation of the coefficient estimates."
  ),
  model_fit="The fit from [lm()] or [AER::ivreg()].",
  grad_return=paste0(
    "A list containing the sensitivity of the coefficients and standard ",
    "errors to data re-weighting.  The list entries are\n",
    "\\describe{ \n",
    "  \\item{model_fit}{The original fit as returned by the fitting function} \n",
    "  \\item{n_obs}{The number of observations given weights} \n",
    "  \\item{regressor_names}{The names of the regressors as strings} \n",
    "  \\item{grad_fun}{The function used to compute the derivatives} \n",
    "  \\item{betahat}{The estimated regression coefficient at the original weights} \n",
    "  \\item{se}{The standard errors at the original weights} \n",
    "  \\item{weights}{The original weights} \n",
    "  \\item{beta_grad}{The gradient of betahat with respect to the weights} \n",
    "  \\item{se_grad}{The gradient of se with respect to the weights} \n",
    "}\n"
  ),
  weights="Optional.  A vector of weights.  If unset, use the original weights.",
  sig_num_ses="How wide the confidence interval is, as a number of standard errors.",
  grad_df=paste0(
    "The gradient dataframe with attributes for a particular coefficient, ",
    "e.g., as returned by [GetTargetRegressorGrads()].")
)

# Usage example:
##' @param lm_result `r docs$lm_result`
