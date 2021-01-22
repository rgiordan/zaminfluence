library(dplyr)
library(stringr)

#' Calcualate the scale (L2 norm) of an influence vector.
#'
#' @param influence_vec: A vector of influence scores.
#' @return $sigma_xi$, the "scale" of the influence function.
#' @export
GetInfluenceScale <- function(influence_vec) {
  return(sqrt(length(influence_vec) * sum(influence_vec ^ 2)))
}


#' Return the $Gamma_alpha$ quantity.
#'
#' @param influence_df `r docs$infl_df`
#' @param alpha `r docs$alpha_val`
#' @return $Gamma_alpha$, the shape coefficient.
#' @export
GetGammaAlpha <- function(influence_df, alpha, alpha_col, gamma_alpha_col) {
  return(approx(x=influence_df[[alpha_col]],
                y=influence_df[[gamma_alpha_col]],
                xout=alpha)$y)
}


FindFirstValueIndex <- function(vec, target) {
  indices <- which(diff(sign(vec - target)) != 0)
  if (length(indices) == 0) {
    return(NA)
  }
  else {
    return(min(indices) + 1)
  }
}


GetAlphaForTarget <- function(influence_df, estimate_col, alpha_col, target) {
  if (influence_df[[estimate_col]][1] == target) {
    # In the unlikely event that the estimate is already equal to the target,
    # return 0.
    alpha <- 0
  } else {
    index <- FindFirstValueIndex(influence_df[[estimate_col]], target)
    if (is.na(index)) {
      # No alpha can achieve the target (up to linearity assumptions, of course)
      alpha <- NA
    } else {
      alpha <- influence_df[[alpha_col]][index]
    }
  }
  result <-
      data.frame(target=target,
                 estimate_col=estimate_col,
                 stringsAsFactors=FALSE)
  result[[alpha_col]] <- alpha
  return(result)
}


GetAlphaForEstimateTarget <- function(influence_dfs, influence_col,
                                      alpha_col, target) {
  base_vals <- SafeGetBaseVals(influence_dfs)
  base_val <- base_vals[influence_col]
  delta <- target - base_val
  delta_sign <- sign(delta)
  direction <- if (delta_sign > 0) "pos" else "neg"
  # influence_df <- if (delta_sign > 0) influence_dfs$pos else influence_dfs$neg
  alpha_df <- GetAlphaForTarget(
    influence_df=influence_dfs[[direction]],
    estimate_col=paste0(influence_col, "_est"),
    alpha_col=alpha_col,
    target=target) %>% mutate(direction=!!direction)
  return(alpha_df)
}


#' @export
GetAlphaForSignChange <- function(influence_dfs, alpha_col) {
  return(GetAlphaForEstimateTarget(
    influence_dfs=influence_dfs,
    influence_col="beta",
    alpha_col=alpha_col,
    target=0.0) %>% mutate(change="sign"))
}


#' @export
GetAlphaForSignificanceChange <- function(influence_dfs, alpha_col, target) {
  LocalGetAlphaForTarget <- function(influence_df, influence_col) {
    GetAlphaForTarget(
      influence_df=influence_df,
      estimate_col=paste0(influence_col, "_est"),
      alpha_col=alpha_col,
      target=0.0)
  }
  base_vals <- SafeGetBaseVals(influence_dfs)
  beta_pzse <- base_vals["beta_pzse"]
  beta_mzse <- base_vals["beta_mzse"]
  beta_est <- base_vals["beta"]

  # Get the change directions
  pzse_dir <- if(beta_pzse > 0) "neg" else "pos"
  mzse_dir <- if(beta_mzse > 0) "neg" else "pos"
  # pzse_df <- if(beta_pzse > 0) influence_dfs$neg else influence_dfs$pos
  # mzse_df <- if(beta_mzse > 0) influence_dfs$neg else influence_dfs$pos

  # Make descriptive names corresponding to each kind of change.
  sig_desc <- "significance"
  #insig_desc <- "insignificance"
  insig_desc <- sig_desc # At Rachael's request
  sign_and_sig_desc <- "sign and significance"

  if (sign(beta_pzse) == sign(beta_mzse)) {
    # The result is significant.
    if (beta_est > 0) {
      # If beta_est > 0 then
      # setting pzse to zero changes the sign of the point estimate
      pzse_desc <- sign_and_sig_desc
      mzse_desc <- insig_desc
    } else {
      # If beta_est < 0 then
      # setting mzse to zero changes the sign of the point estimate
      pzse_desc <- insig_desc
      mzse_desc <- sign_and_sig_desc
    }
  } else {
    # The result is insignificant.
    if (beta_est > 0) {
      # If beta_est > 0 then setting
      # pzse to zero changes the sign of the point estimate
      pzse_desc <- sign_and_sig_desc
      mzse_desc <- sig_desc
    } else {
      # If beta_est < 0 then setting
      # mzse to zero changes the sign of the point estimate
      pzse_desc <- sig_desc
      mzse_desc <- sign_and_sig_desc
    }
  }

  result <- bind_rows(
    LocalGetAlphaForTarget(influence_dfs[[pzse_dir]], "beta_pzse") %>%
      mutate(change=!!pzse_desc, direction=!!pzse_dir),
    LocalGetAlphaForTarget(influence_dfs[[mzse_dir]], "beta_mzse") %>%
      mutate(change=!!mzse_desc, direction=!!mzse_dir)
  )

  return(result)
}

#' Estimate the number and proportion of datapoints needed to effect various
#' changes.
#'
#' @param influence_dfs `r docs$influence_dfs`
#' @param alpha_col `r docs$alpha_col`
#'
#' @return A dataframe summarizing the estimated number of points needed
#' to effect changes in sign, significance, or both sign and significance,
#' or `NA` if the linear approximation cannot produce such a change.
#'
#' @export
GetRegressionTargetChange <- function(influence_dfs, alpha_col) {
    return(bind_rows(
      GetAlphaForSignChange(influence_dfs$sign, alpha_col),
      GetAlphaForSignificanceChange(influence_dfs$sig, alpha_col)))
}
