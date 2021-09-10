
######################
# Functions to help re-running to check the approximations.



#######################################
# Functions for ordinary regression

#' Rerun the regression with a new subset of rows.
#'@param w_bool A boolean vector of rows to keep in the original dataframe.
#'@param lm_result `r docs$lm_result`
#'@param se_group `r docs$se_group`
#'@param save_w Optional. If `TRUE`, save the new weight vector in the output.
#'
#'@return `r docs$rerun_return`
#'@export
RerunRegression <- function(w_bool, lm_result, se_group=NULL, save_w=FALSE) {
  new_w <- rep(0.0, length(w_bool))
  if ("weights" %in% names(lm_result)) {
    new_w[w_bool] <- lm_result$weights[w_bool]
  } else {
    new_w[w_bool] <- 1.0
  }

  # Rerun using my own code; I don't want to deal with how R handles the
  # scoping of the weight variables in the regression.
  ret_list <-
    ComputeRegressionResults(lm_result, weights=new_w, se_group=se_group)
  if (save_w) {
    ret_list$w <- new_w
  }
  return(ret_list)
}


#######################################
# Functions for IV regression


#' Rerun the regression with a new subset of rows.
#'@param w_bool A boolean vector of rows to keep in the original dataframe.
#'@param iv_res `r docs$iv_res`
#'@param se_group `r docs$se_group`
#'@param save_w Optional. If TRUE, save the new weight vector in the output.
#'
#'@return `r docs$rerun_return`
#'
#'@export
RerunIVRegression <- function(w_bool, iv_res, se_group=NULL, save_w=FALSE) {
  # Rerun using my own code; I don't want to deal with how R handles the
  # scoping of the weight variables in the regression.
  if (length(iv_res$y) != length(w_bool)) {
    stop(paste0("``w_bool`` is not the same length as the regression data. ",
                "Note that re-running regression with aggregated ",
                "influence functions is not implemented."))
  }

  new_w <- rep(0.0, length(w_bool))
  if (is.null(iv_res[["weights"]])) {
    new_w[w_bool] <- 1.0
  } else {
    new_w[w_bool] <- iv_res$weights[w_bool]
  }

  ret_list <-
    ComputeIVRegressionResults(iv_res, weights=new_w, se_group=se_group)
  if (save_w) {
    ret_list$w <- new_w
  }
  return(ret_list)
}
