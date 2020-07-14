
#' @export
GetBoundQuantities <- function(reg_infl, sorted_rows) {
  x_mat <- reg_infl$lm_result$x[sorted_rows, ]
  x_dim <- ncol(x_mat)
  n_row <- nrow(x_mat)
  eps_hat <- reg_infl$lm_result$residuals[sorted_rows]

  xtx <- t(x_mat) %*% x_mat
  hhat <- xtx / nrow(x_mat)
  xtx_cumsum <- matrix(0, x_dim, x_dim)
  xeps_cumsum <- rep(0, x_dim)
  s_vec <- c()
  r_vec <- c()
  for (n in 1:nrow(x_mat)) {
    this_x <- x_mat[n, ]
    xtx_cumsum <- xtx_cumsum + outer(this_x, this_x)
    xeps_cumsum <- xeps_cumsum + this_x * eps_hat[n]
    s_vec <- c(s_vec, sqrt(sum(xtx_cumsum^2)) / n)
    r_vec <- c(r_vec, sqrt(sum(xeps_cumsum^2)) / n)
  }

  c_op <- 1 / min(eigen(hhat)$values)
  m1 <- sqrt(sum(hhat^2))

  # Get the IJ estimates
  beta_ij_diff <- -1 * apply(reg_infl$beta_grad[, sorted_rows],
                             cumsum, MARGIN=1)
  betahat <- reg_infl$lm_result$coefficients
  beta_ij_est <- beta_ij_diff + rep(betahat, each=nrow(beta_ij_diff))
  colnames(beta_ij_est) <- paste0("betahat_", 1:x_dim)
  diff_l2_norm <- sapply(1:nrow(beta_ij_diff),
                         function(n) { sqrt(sum(beta_ij_diff[n, ]^2)) })

  # Get a function to return the bound for a particular b, rho, and maximum
  # number of left-out points.
  GetBound <- function(b, rho, max_n_lo) {
    n_lo <- 1:max_n_lo
    s <- s_vec[n_lo]
    r <- r_vec[n_lo]
    alpha_vec <- n_lo / n_row
    delta0_vec <- alpha_vec * (r + sqrt(x_dim) * b * s)
    delta0_lb_vec <- alpha_vec * r
    delta1_vec <- alpha_vec * s

    c_set <- c_op * alpha_vec * s
    ctwiddle_op <- c_op / (1 - rho)
    ij_err <- 2 * (ctwiddle_op^2) * delta1_vec * delta0_vec
    ij_err_lb <- 2 * (ctwiddle_op^2) * delta1_vec * delta0_lb_vec

    ball_width <- diff_l2_norm[1:max_n_lo] + ij_err

    max_valid_n_lo_b <- max(which(ball_width < b))
    max_valid_n_lo_set <- max(which(c_set < rho))
    max_valid_n_lo <- min(c(max_valid_n_lo_set, max_valid_n_lo_b))

    return(list(
      alpha_vec=alpha_vec,
      c_set=c_set,
      ij_err=ij_err,
      ij_err_lb=ij_err_lb,
      ball_width=ball_width,
      diff_l2_norm=diff_l2_norm[1:max_n_lo],

      max_valid_n_lo_b=max_valid_n_lo_b,
      max_valid_n_lo_set=max_valid_n_lo_set,
      max_valid_n_lo=max_valid_n_lo,

      r=r,
      s=s,
      valid=c_set < rho,
      delta0=delta0_vec,
      delta1=delta1_vec,
      ctwiddle_op=ctwiddle_op,
      rho=rho,
      b=b))
  }

  return(list(
    GetBound=GetBound,
    beta_ij_est=beta_ij_est,
    beta_ij_diff=beta_ij_diff,
    diff_l2_norm=diff_l2_norm,
    c_op=c_op,
    m1=m1,
    s_vec=s_vec,
    r_vec=r_vec,
    hhat=hhat,
    x_dim=x_dim
  ))
}
