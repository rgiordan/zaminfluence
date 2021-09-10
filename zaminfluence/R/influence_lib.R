

#' @export
ProcessInfluenceVector <- function(infl, base_value, num_obs, obs_per_row=1) {
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


#' @export
GetAPIP <- function(infl_lists, signal) {
    # To produce a negative change, drop observations with positive influence
    # scores, and vice-versa.
    infl_list <- if (signal < 0) infl_lists$pos else infl_lists$neg
    n_vec <- 1:length(infl_list$infl_cumsum)
    # TODO: do this more efficiently using your own routine, since
    # we know that infl_cumsum is increasing?
    n_drop <- approx(x=-1 * c(0, infl_list$infl_cumsum),
                     y=c(0, n_vec),
                     xout=signal)$y %>% ceiling()
    if (is.na(n_drop)) {
        drop_inds <- NA
    } else if (n_drop == 0) {
        drop_inds <- c()
    } else {
        drop_inds <- infl_list$infl_inds[1:n_drop]
    }
    return(list(
        n=n_drop,
        prop=n_drop / infl_list$num_obs,
        inds=drop_inds
    ))
}


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


#' @export
GetAMIS <- function(infl_list, n_drop=NULL, prop_drop=NULL) {
  if (sum(c(!is.null(n_drop), !is.null(prop_drop)) == 1)) {
    stop("Exactly one of `n_drop` and `prop_drop` can be specified.")
  }
  if (!is.null(prop_drop)) {
    n_drop <- ceiling(prop_drop * infl_list$num_obs)
  }
  return(infl_list$inds[1:n_drop])
}


#' @export
GetAMIP <- function(infl_list, n_drop=NULL, prop_drop=NULL) {
  if (sum(c(!is.null(n_drop), !is.null(prop_drop)) == 1)) {
    stop("Only one of `n_drop` and `prop_drop` can be specified.")
  }
  if (!is.null(prop_drop)) {
    n_drop <- ceiling(prop_drop * infl_list$num_obs)
  }
  if (n_drop == 0) {
    return(0)
  } else if (n_drop > infl_list$num_obs){
    stop(sprintf(
      paste0("The number of observations to drop (%d) must be no greater ",
             "than the total number of observations (%d)."),
      n_drop, infl_list$num_obs))
  }
  return(sum(infl_list$infl_cumsum[n_drop]))
}


#'@export
PredictChange <- function(qoi_infl, drop_inds) {
    return(sum(qoi_infl$infl[drop_inds]))
}
