# A list of strings for arguments and return values that repeatedly recur
# in the docs.

docs <- list(
  infl_df=paste0(
      "A single influence dataframe, e.g., a single member of the list ",
      "returned by [SortAndAccumulate()]."
  ),
  alpha_colname=paste0(
      "A string with the name of the target alpha column (typically ",
      "`prop_removed` or `num_removed`)"
  ),
  alpha_val="The alpha value to target corresponding to `alpha_colname`.",
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
    "e.g., as returned by [GetTargetRegressorGrads()]."),
  influence_cols_default=paste0(
      "Optional.  Additional influence columns to be ",
      "cumulatively summed after sorting.  ",
      "Effectively, a linear approximation to the change in the ",
      "values in `influence_cols` is computed.  ",
      "By default all the `base_vals` ",
      "gradients are summed."),
  influence_cols=paste0(
      "Which influence columns are to be cumulatively summed after sorting."),
  influence_dfs=paste0(
    "A list of sorted influence dataframes as returned, e.g., by ",
    "[SortAndAccumulate()].  The list is expected to have results sorted ",
    "for both negative and positive changes in both sign and signifiance."),
  target_change=paste0(
    "A dataframe of target changes, e.g., as returned by ",
    "[GetRegressionTargetChange()] ")
)

# Usage example:
##' @param lm_result `r docs$lm_result`
