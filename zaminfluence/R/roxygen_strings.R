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
  model_fit="The fit from [lm()] or [AER::ivreg()]."
)
