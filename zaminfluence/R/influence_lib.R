
#' Process a vector of influence scores to produce sorted influence scores.
#' @param infl A vector of influence scores for a quantity of interest,
#' in the same order as the original data.
#' @param base_value The value of the quantity of interest at the original fit.
#'
#' @return See "Quantity of Interest" in README.md
#' @export
ProcessInfluenceVector <- function(infl, base_value, num_obs=NULL, obs_per_row=1) {
    if (is.null(num_obs)) {
        num_obs <- length(infl)
    }
    infl_pos <- infl > 0
    infl_neg <- infl < 0

    inds_pos <- (1:length(infl))[infl_pos]
    inds_neg <- (1:length(infl))[infl_neg]

    ordered_inds_pos <- inds_pos[order(-1 * infl[infl_pos])]
    ordered_inds_neg <- inds_neg[order(infl[infl_neg])]

    infl_pos <- infl[ordered_inds_pos]
    infl_cumsum_pos <- cumsum(infl_pos)

    infl_neg <- infl[ordered_inds_neg]
    infl_cumsum_neg <- cumsum(infl_neg)
    return(list(
        neg=list(infl_inds=ordered_inds_neg,
                 infl_cumsum=infl_cumsum_neg,
                 num_obs=num_obs,
                 obs_per_row=obs_per_row),
        pos=list(infl_inds=ordered_inds_pos,
                 infl_cumsum=infl_cumsum_pos,
                 num_obs=num_obs,
                 obs_per_row=obs_per_row),
        base_value=base_value,
        infl=infl  # In the original order
        ))
}


#' Compute the approximate perturbation-inducing proportion (APIP).
#' @param qoi ``r docs$qoi``
#' @param signal `r docs$signal`
#'
#' @return A list containing
#' `n`: The number of points to drop
#' `prop`: The proportion of points to drop
#' `inds`: `r docs$drop_inds`
#' @export
GetAPIP <- function(qoi, signal) {
    # To produce a negative change, drop observations with positive influence
    # scores, and vice-versa.
    qoi_sign <- if (signal < 0) qoi$pos else qoi$neg
    n_vec <- 1:length(qoi_sign$infl_cumsum)
    # TODO: do this more efficiently using your own routine, since
    # we know that infl_cumsum is increasing?
    n_drop <- approx(x=-1 * c(0, qoi_sign$infl_cumsum),
                     y=c(0, n_vec),
                     xout=signal)$y %>% ceiling()
    if (is.na(n_drop)) {
        drop_inds <- NA
    } else if (n_drop == 0) {
        drop_inds <- c()
    } else {
        drop_inds <- qoi_sign$infl_inds[1:n_drop]
    }
    return(list(
        n=n_drop,
        prop=n_drop / qoi_sign$num_obs,
        inds=drop_inds
    ))
}


#' Compute a weight vector for a set of dropped indices.
#'
#' @param drop_inds ``r docs$drop_inds``
#' @param num_obs The number of observations in the original data.
#' @param bool (Optional)  If true, return a boolean vector.  Otherwise,
#' return a numeric vector (with ones and zeros).
#' @param invert (Optional) If TRUE, return a vector that retains the
#' observations in `drop_inds`.  Default is FALSE.
#'
#' @return A vector of weights in the order of the original data.
#' @export
GetWeightVector <- function(drop_inds, num_obs, bool=FALSE, invert=FALSE) {
  if (bool) {
    w <- rep(TRUE, num_obs)
    w[drop_inds] <- FALSE
    if (invert) {
      return(!w)
    } else {
      return(w)
    }
  } else { # Integers, not boolean weights
    w <- rep(1, num_obs)
    w[drop_inds] <- 0
    if (invert) {
      return(1 - w)
    } else {
      return(w)
    }
  }
}


#' Compute the approximate maximally-influential set (AMIS).
#' @param qoi ``r docs$qoi``
#' @param n_drop The number of points to drop (we will round up).
#'
#' @return `r docs$drop_inds`
#' @export
GetAMIS <- function(qoi, sign, n_drop) {
  if (!(sign %in% c("pos", "neg"))) {
    stop("Sign must be either `pos` or `neg`.")
  }
  if (n_drop < 0) {
    stop("`n_drop` must be non-negative.")
  }
  n_drop <- ceiling(n_drop)
  n_scores <- length(qoi[[sign]]$infl_inds)
  if (n_drop > n_scores) {
    warning(sprintf(paste0(
      "More dropped observations were requested than influence scores",
      "present (%d > %d).  Returning all influence scores ",
      "of the specified sign."), n_drop, n_scores))
      return(qoi[[sign]]$infl_inds)
  }
  if (n_drop > 0) {
    return(qoi[[sign]]$infl_inds[1:n_drop])
  } else {
    return(c())
  }
}


#' Compute the approximate maximum influence perturbation (AMIP).
#' @param qoi ``r docs$qoi``
#' @param n_drop The number of points to drop (we will round up).
#'
#' @return The approximate largest change that can be produced by dropping
#' the specified number of points.
#' @export
GetAMIP <- function(qoi, sign, n_drop) {
  if (n_drop == 0) {
    return(0)
  }
  amis <- GetAMIS(qoi, sign, n_drop)
  return(PredictChange(qoi, amis))
}

#' Predict the effect of dropping points.
#' @param qoi ``r docs$qoi``
#' @param drop_inds ``r docs$drop_inds``
#'
#' @return A linear approximation to the effect of dropping the specified
#' observations.
#'@export
PredictChange <- function(qoi, drop_inds) {
    return(-1 * sum(qoi$infl[drop_inds]))
}