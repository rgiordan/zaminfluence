
######################
# Functions to help re-running to check the approximations.


#' Get the row of an influence dataframe that cooresponds to a particular alpha.
#' @param infl_df A single influence dataframe.
#' @param alpha_colname A string with the name of the target alpha column
#' @param alpha_val The alpha value to target
#'@export
GetAlphaRow <- function(infl_df, alpha_colname, alpha_val) {
  # This will round up the number of rows to be removed.
  return(min(which(infl_df[[alpha_colname]] >= alpha_val)))
}


#' Get a set of rows targeting a particular value of alpha.
#' @param infl_df A single influence dataframe.
#' @param alpha_colname A string with the name of the target alpha column
#' @param alpha_val The alpha value to target
#' @param boolean Optional.  If TRUE, return a boolean vector.  If FALSE,
#' return a vector of integer indices.
#' @param rows_to_keep Optional.  If TRUE, return a vector of rows to keep
#' for the target level of alpha.  If FALSE, return a vector of rows to drop.
#' @return A row vector for the original dataframe corresponding to the target
#' level of alpha.  Note that the number of left-out rows will be rounded up.
#'@export
GetWeightForAlpha <- function(infl_df, alpha_colname, alpha_val,
                              boolean=TRUE, rows_to_keep=TRUE) {
  # Some checks.
  if (min(diff(infl_df[[alpha_colname]])) < 0) {
    stop("infl_df must be sorted by the column <alpha_colname>.")
  }
  n_grad_rows <- attr(infl_df, "n_grad_rows")
  data_row_cols <- attr(infl_df, "data_row_cols")

  # Get how many rows need to be removed for the target alpha.
  num_rows_removed <- GetAlphaRow(infl_df, alpha_colname, alpha_val)

  # Get the indices of the removed rows from the original data.
  row_inds <-
    lapply(data_row_cols,
           function(row_col) { infl_df[[row_col]][1:num_rows_removed] }) %>%
    unlist()

  # Recall that there is an NA for no rows removed in the influence dataframe.
  # Let's remove it.
  row_inds <- row_inds[!is.na(row_inds)]

  # Process the output.
  if (boolean) {
      w_keep_bool <- rep(TRUE, n_grad_rows)
      w_keep_bool[row_inds] <- FALSE
      if (rows_to_keep) {
        # Return a boolean vector selecting rows to keep from the original dataset.
        return(w_keep_bool)
      } else {
        # Return a boolean vector selecting rows to drop from the original dataset.
        return(!w_keep_bool)
      }
  } else {
    if (rows_to_keep) {
      # Return a vector of indices to keep from the original dataset.
      w_int <- setdiff(1:n_grad_rows, row_inds)
      return(w_int)
    } else {
      # Return a vector of indices to drop from the original dataset.
      return(row_inds)
    }
  }
}


#######################################
# Functions for ordinary regression

#' Rerun the regression with a new subset of rows.
#'@param w_bool A boolean vector of rows to keep in the original dataframe.
#'@param lm_result The original regression result.
#'@param se_group Optional. The standard error grouping variable.
#'@param save_w Optional. If TRUE, save the new weight vector in the output.
#'@return A list containing the new regression estimate, standard error
#' covariance, and standard errors.
#'@export
RerunRegression <- function(w_bool, lm_result, se_group=NULL, save_w=FALSE) {
  if (length(lm_result$y) != length(w_bool)) {
    stop(paste0("``w_bool`` is not the same length as the regression data. ",
                "Note that re-running regression with aggregated ",
                "influence functions is not implemented."))
  }
  py_main <- SetPythonRegressionVariables(lm_result, se_group=se_group)

  new_w <- rep(0.0, nrow(py_main$x))
  if ("weights" %in% names(lm_result)) {
    new_w[w_bool] <- lm_result$weights[w_bool]
  } else {
    new_w[w_bool] <- 1.0
  }
  py_main$new_w <- as.array(as.numeric(new_w))
  reticulate::py_run_string("betahat_w = regsens_rgiordandev.reg(y, x, w=new_w)")
  reticulate::py_run_string("
se_cov = regsens_rgiordandev.get_standard_error_matrix(
  betahat_w, y, x, w=new_w, se_group=se_group)")
  ret_list <- list(
    betahat=py_main$betahat_w,
    se_cov=py_main$se_cov,
    se=sqrt(diag(py_main$se_cov)))
  if (save_w) {
    ret_list$w <- py_main$new_w
  }
  return(ret_list)
}


#' Rerun the target regression for a particular target alpha.
#'@param infl_df A particular influence dataframe.
#'@param lm_result The original regression result.
#'@param alpha_colname A string with the name of the target alpha column
#'@param alpha_val The alpha value to target
#'@param se_group Optional. The standard error grouping variable.
#'@return A dataframe with variables from the re-run regression.
#'@export
RerunTargetRegressorForAlpha <- function(
    infl_df, lm_result, alpha_colname, alpha_val, se_group=NULL) {

  rerun_result <-
    GetWeightForAlpha(infl_df, alpha_colname, alpha_val) %>%
    RerunRegression(lm_result, se_group=se_group)
  target_index <- attr(infl_df, "target_index")
  sig_num_ses <- attr(infl_df, "sig_num_ses")
  result_df <-
    data.frame(target_index=target_index,
               beta=rerun_result$betahat[target_index],
               se=rerun_result$se[target_index]) %>%
    mutate(beta_pzse=beta + !!sig_num_ses * se,
           beta_mzse=beta - !!sig_num_ses * se)
  result_df[[alpha_colname]] <- alpha_val
  return(result_df)
}


#######################################
# Functions for IV regression


#' Rerun the regression with a new subset of rows.
#'@param w_bool A boolean vector of rows to keep in the original dataframe.
#'@param iv_res The original IV regression result.
#'@param se_group Optional. The standard error grouping variable.
#'@param save_w Optional. If TRUE, save the new weight vector in the output.
#'@return A list containing the new regression estimate, standard error
#' covariance, and standard errors.
#'@export
RerunIVRegression <- function(w_bool, iv_res, se_group=NULL, save_w=FALSE) {
  if (length(iv_res$y) != length(w_bool)) {
    stop(paste0("``w_bool`` is not the same length as the regression data. ",
                "Note that re-running regression with aggregated ",
                "influence functions is not implemented."))
  }
  py_main <- SetPythonIVRegressionVariables(iv_res, se_group=se_group)

  new_w <- rep(0.0, nrow(py_main$x))
  if (is.null(iv_res[["weights"]])) {
    new_w[w_bool] <- 1.0
  } else {
    new_w[w_bool] <- iv_res$weights[w_bool]
  }
  py_main$new_w <- as.array(as.numeric(new_w))
  reticulate::py_run_string("
betahat_w = iv_lib.iv_reg(y, x, z, w=new_w)")
  reticulate::py_run_string("
se_cov = iv_lib.get_iv_standard_error_matrix(
  betahat_w, y, x, z, w=new_w, se_group=se_group)")
  ret_list <- list(
    betahat=py_main$betahat_w,
    se_cov=py_main$se_cov,
    se=sqrt(diag(py_main$se_cov)))
  if (save_w) {
    ret_list$w <- py_main$new_w
  }
  return(ret_list)
}



#' Rerun the target regression for a particular target alpha.
#'@param infl_df A particular influence dataframe.
#'@param iv_res The original IV regression result.
#'@param alpha_colname A string with the name of the target alpha column
#'@param alpha_val The alpha value to target
#'@param se_group Optional. The standard error grouping variable.
#'@return A dataframe with variables from the re-run regression.
#'@export
RerunTargetIVRegressionForAlpha <- function(
    infl_df, iv_res, alpha_colname, alpha_val, se_group=NULL) {

  rerun_result <-
    GetWeightForAlpha(infl_df, alpha_colname, alpha_val) %>%
    RerunIVRegression(iv_res, se_group=se_group)
  target_index <- attr(infl_df, "target_index")
  sig_num_ses <- attr(infl_df, "sig_num_ses")
  result_df <-
    data.frame(target_index=target_index,
               beta=rerun_result$betahat[target_index],
               se=rerun_result$se[target_index]) %>%
    mutate(beta_pzse=beta + !!sig_num_ses * se,
           beta_mzse=beta - !!sig_num_ses * se)
  result_df[[alpha_colname]] <- alpha_val
  return(result_df)
}



####################################
# A wrapper that works for both

#' Rerun the target regression for a particular target alpha.
#'@param infl_df A particular influence dataframe.
#'@param model_fit The fit from `lm` or `AER:ivreg`.
#'@param alpha_colname A string with the name of the target alpha column
#'@param alpha_val The alpha value to target
#'@param se_group Optional. The standard error grouping variable.
#'@return A dataframe with variables from the re-run model.
#'@export
RerunTargetModelForAlpha <- function(
    infl_df, model_fit, alpha_colname, alpha_val, se_group=NULL) {

    model_class <- class(model_fit)
    if (model_class == "lm") {
      return(RerunTargetRegressorForAlpha(
        infl_df, model_fit, alpha_colname, alpha_val, se_group))
    } else if (model_class == "ivreg") {
      return(RerunTargetIVRegressionForAlpha(
        infl_df, model_fit, alpha_colname, alpha_val, se_group))
    } else {
      stop(sprint("Unknown model class %s.  Valid classes are %s",
                  model_class,
                  paste(c("lm", "ivreg", collapse=", "))))
    }
}


#' Rerun the target regression for every change in target_change.
#'@param influence_dfs A list of influence dataframes
#'@param target_change The output of GetRegressionTargetChange
#'@param model_fit The fit from `lm` or `AER:ivreg`.
#'@return A dataframe with variables from all the re-run models.
#'@export
RerunForTargetChanges <- function(influence_dfs, target_change, model_fit,
                                  se_group=NULL) {
    # Get whether proportion or number removed was used.
    alpha_col <- intersect(c("num_removed", "prop_removed"), names(target_change))
    if (length(alpha_col) >= 2) {
        stop("target_change should have only one of `num_removed` or `prop_removed`.")
    }
    if (length(alpha_col) == 0) {
        stop("target_change must have a column `num_removed` or `prop_removed`.")
    }

    # Which influence dataframe to get for each change.
    target_df <- list()
    target_df[["sign"]] <- "sign"
    target_df[["significance"]] <- "sig"
    target_df[["sign and significance"]] <- "sig"

    changes <- target_change[ !is.na(target_change[[alpha_col]]), "change"]
    if (length(changes) == 0) {
        return(NULL)
    }

    # base_vals are the original values, and should be the same for every
    # influence function.
    base_vals <- attr(influence_dfs$sign$pos, "base_vals")
    target_index <- attr(influence_dfs$sign$pos, "target_index")
    results_list <- list()
    results_list[[1]] <-
        as.data.frame(t(base_vals)) %>%
        mutate(target_index=!!target_index, change="original")
    for (change in changes) {
        this_target_change <- filter(target_change, change == !!change)
        alpha <- this_target_change[[alpha_col]]
        direction <- this_target_change$direction
        rerun_result <- RerunTargetModelForAlpha(
              infl_df=influence_dfs[[ target_df[[change]] ]][[direction]],
              model_fit=model_fit,
              alpha_colname=alpha_col,
              alpha_val=alpha,
              se_group=se_group) %>%
            mutate(change=!!change)
        results_list[[length(results_list) + 1]] <- rerun_result
    }
    results <- do.call(bind_rows, results_list)
    return(results)
}
