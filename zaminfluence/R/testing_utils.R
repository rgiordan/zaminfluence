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



GenerateRandomEffects <- function(n_obs, num_groups=NULL) {
  if (!is.null(num_groups)) {
    # Add random effects.  se_group must be zero-indexed.
    se_group <- floor(num_groups * runif(n_obs))
    se_group[se_group == num_groups] <- 0  # Better safe than sorry

    re <- rnorm(num_groups)[se_group + 1]
    return(data.frame(se_group=se_group, re=re))
  } else {
    return(data.frame(re=rep(0, n_obs)))
  }
}


#' @export
GenerateRegressionData <- function(n_obs, beta_true, x=NULL, num_groups=NULL) {
  x_dim <- length(beta_true)
  if (is.null(x)) {
    x <- matrix(runif(n_obs * x_dim), n_obs, x_dim)
    x <- x - rep(colMeans(x), each=n_obs)
  } else {
    if (!(nrow(x) == n_obs)) {
      stop("Wrong number of x rows.")
    }
  }
  eps <- rnorm(n_obs)
  re_df <- GenerateRandomEffects(n_obs, num_groups)
  y <- x %*% beta_true + eps + re_df$re
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
GenerateIVRegressionData <- function(n_obs, beta_true, num_groups=NULL) {
  # Simulate some IV data

  x_dim <- length(beta_true)
  x <- rnorm(n_obs * x_dim) %>% matrix(nrow=n_obs)
  x_rot <- diag(x_dim) + rep(0.2, x_dim ^ 2) %>% matrix(nrow=x_dim)
  x <- x %*% x_rot
  x <- x - rep(colMeans(x), each=n_obs)
  colMeans(x)

  z <- rnorm(n_obs * x_dim) %>% matrix(nrow=n_obs)
  z_rot <- diag(x_dim) + rep(0.2, x_dim ^ 2) %>% matrix(nrow=x_dim)
  z <- z %*% z_rot + x
  z <- z - rep(colMeans(z), each=n_obs)
  colMeans(z)

  Project <- function(z, vec) {
      n_obs <- nrow(z)
      ztz <- t(z) %*% z / n_obs
      ztv <- t(z) %*% vec / n_obs
      return(z %*% solve(ztz, ztv))
  }

  ProjectPerp <- function(z, vec) {
      return(vec - Project(z, vec))
  }

  re_df <- GenerateRandomEffects(n_obs, num_groups)

  sigma_true <- 2.0
  eps_base <- rnorm(n_obs) + rowSums(x) + re_df$re
  eps_true <- ProjectPerp(z, eps_base)
  y <- x %*% beta_true + eps_true


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
