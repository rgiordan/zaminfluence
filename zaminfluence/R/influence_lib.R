
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
    qoi <- if (signal < 0) qois$pos else qois$neg
    n_vec <- 1:length(qoi$infl_cumsum)
    # TODO: do this more efficiently using your own routine, since
    # we know that infl_cumsum is increasing?
    n_drop <- approx(x=-1 * c(0, qoi$infl_cumsum),
                     y=c(0, n_vec),
                     xout=signal)$y %>% ceiling()
    if (is.na(n_drop)) {
        drop_inds <- NA
    } else if (n_drop == 0) {
        drop_inds <- c()
    } else {
        drop_inds <- qoi$infl_inds[1:n_drop]
    }
    return(list(
        n=n_drop,
        prop=n_drop / qoi$num_obs,
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
#' Exactly one of `n_drop` and `prop_drop` must be specified.
#' @param prop_drop The proportion of points to drop (we will round up).
#' Exactly one of `n_drop` and `prop_drop` must be specified.
#'
#' @return `r docs$drop_inds`
#' @export
GetAMIS <- function(qoi, n_drop=NULL, prop_drop=NULL) {
  if (sum(c(!is.null(n_drop), !is.null(prop_drop)) == 1)) {
    stop("Exactly one of `n_drop` and `prop_drop` can be specified.")
  }
  if (!is.null(prop_drop)) {
    n_drop <- ceiling(prop_drop * qoi$num_obs)
  }
  return(qoi$inds[1:n_drop])
}


#' Compute the approximate maximum influence perturbation (AMIP).
#' @param qoi ``r docs$qoi``
#' @param n_drop The number of points to drop (we will round up).
#' Exactly one of `n_drop` and `prop_drop` must be specified.
#' @param prop_drop The proportion of points to drop (we will round up).
#' Exactly one of `n_drop` and `prop_drop` must be specified.
#'
#' @return The approximate largest change that can be produced by dropping
#' the specified number of points.
#' @export
GetAMIP <- function(qoi, n_drop=NULL, prop_drop=NULL) {
  if (sum(c(!is.null(n_drop), !is.null(prop_drop)) == 1)) {
    stop("Only one of `n_drop` and `prop_drop` can be specified.")
  }
  if (!is.null(prop_drop)) {
    n_drop <- ceiling(prop_drop * qoi$num_obs)
  }
  if (n_drop == 0) {
    return(0)
  } else if (n_drop > qoi$num_obs){
    stop(sprintf(
      paste0("The number of observations to drop (%d) must be no greater ",
             "than the total number of observations (%d)."),
      n_drop, qoi$num_obs))
  }
  return(sum(qoi$infl_cumsum[n_drop]))
}

#' Predict the effect of dropping points.
#' @param qoi ``r docs$qoi``
#' @param drop_inds ``r docs$drop_inds``
#'
#' @return A linear approximation to the effect of dropping the specified
#' observations.
#'@export
PredictChange <- function(qoi, drop_inds) {
    return(sum(qoi$infl[drop_inds]))
}
