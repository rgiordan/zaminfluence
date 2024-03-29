#library(broom)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(reticulate)
library(latex2exp)
library(testthat)


StopIfNotNumericScalar <- function(x) {
  stopifnot(is.numeric(x))
  stopifnot(length(x) == 1)
}


AssertNearlyEqual <- function(x, y, tol=1e-9, desc=NULL) {
  diff_norm <- max(abs(x - y))
  if (is.null(desc)) {
    info_str <- sprintf("%e > %e", diff_norm, tol)
  } else {
    info_str <- sprintf("%s: %e > %e", desc, diff_norm, tol)
  }
  expect_true(diff_norm < tol, info=info_str)
}


AssertNearlyZero <- function(x, tol=1e-15, desc=NULL) {
  x_norm <- max(abs(x))
  if (is.null(desc)) {
    info_str <- sprintf("%e > %e", x_norm, tol)
  } else {
    info_str <- sprintf("%s: %e > %e", desc, x_norm, tol)
  }
  expect_true(x_norm < tol, info=info_str)
}



# This is the version of the sandwich covariance that we compute
vcovWrap <- function(obj, cluster=NULL) {
  vcovCL(obj, cluster=cluster, type="HC0", cadjust=FALSE)
}


# Wrap vcovWrap so that se_group can be null.
GetFitCovariance <- function(fit, se_group=NULL) {
  if (is.null(se_group)) {
    return(vcov(fit))
  } else {
    return(vcovCL(fit, cluster=se_group, type="HC0", cadjust=FALSE))
  }
}



GenerateRandomEffects <- function(num_obs, num_groups=NULL) {
  if (!is.null(num_groups)) {
    # Add random effects.  se_group must be zero-indexed.
    se_group <- floor(num_groups * runif(num_obs))
    se_group[se_group == num_groups] <- 0  # Better safe than sorry

    re <- rnorm(num_groups)[se_group + 1]
    return(data.frame(se_group=se_group, re=re))
  } else {
    return(data.frame(re=rep(0, num_obs)))
  }
}


#' @export
GenerateRegressionData <- function(num_obs, param_true, x=NULL, num_groups=NULL) {
  x_dim <- length(param_true)
  if (is.null(x)) {
    x <- matrix(runif(num_obs * x_dim), num_obs, x_dim)
    x <- x - rep(colMeans(x), each=num_obs)
  } else {
    if (!(nrow(x) == num_obs)) {
      stop("Wrong number of x rows.")
    }
  }
  eps <- rnorm(num_obs)
  re_df <- GenerateRandomEffects(num_obs, num_groups)
  y <- x %*% param_true + eps + re_df$re
  df <- data.frame(x)
  names(df)  <- paste0("x", 1:x_dim)
  df$y <- y
  df$eps <- eps
  df$re <- re_df$re
  if (!is.null(num_groups)) {
    df$se_group <- re_df$se_group
  }
  return(df)
}


#' @export
GenerateIVRegressionData <- function(num_obs, param_true, num_groups=NULL) {
  # Simulate some IV data

  x_dim <- length(param_true)
  x <- rnorm(num_obs * x_dim) %>% matrix(nrow=num_obs)
  x_rot <- diag(x_dim) + rep(0.2, x_dim ^ 2) %>% matrix(nrow=x_dim)
  x <- x %*% x_rot
  x <- x - rep(colMeans(x), each=num_obs)
  colMeans(x)

  z <- rnorm(num_obs * x_dim) %>% matrix(nrow=num_obs)
  z_rot <- diag(x_dim) + rep(0.2, x_dim ^ 2) %>% matrix(nrow=x_dim)
  z <- z %*% z_rot + x
  z <- z - rep(colMeans(z), each=num_obs)
  colMeans(z)

  Project <- function(z, vec) {
      num_obs <- nrow(z)
      ztz <- t(z) %*% z / num_obs
      ztv <- t(z) %*% vec / num_obs
      return(z %*% solve(ztz, ztv))
  }

  ProjectPerp <- function(z, vec) {
      return(vec - Project(z, vec))
  }

  re_df <- GenerateRandomEffects(num_obs, num_groups)

  sigma_true <- 2.0
  eps_base <- rnorm(num_obs) + rowSums(x) + re_df$re
  eps_true <- ProjectPerp(z, eps_base)
  y <- x %*% param_true + eps_true


  x_df <- data.frame(x)
  x_names <- sprintf("x%d", 1:x_dim)
  names(x_df) <- x_names

  z_df <- data.frame(z)
  z_names <- sprintf("z%d", 1:x_dim)
  names(z_df) <- z_names

  df <- bind_cols(x_df, z_df) %>% mutate(y=!!y)

  if (!is.null(num_groups)) {
    df$se_group <- re_df$se_group
  }
  return(df)
}
