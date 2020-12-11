
GetIVSEDerivs <- function(x, z, y, beta, w0, se_group=NULL, testing=FALSE) {
    if (!is.null(se_group)) {
        stop("Not implemented.")
    }

    z_w <- z * w0
    zwz <- t(z_w) %*% z
    zwx <- t(z_w) %*% x

    eps <- as.numeric(y - x %*% beta)
    sig2_hat <- sum(w0 * eps^2) / (num_obs - length(beta))

    zwx_inv_zwz <- solve(zwx, zwz)
    sand_mat <- solve(zwx, t(zwx_inv_zwz))
    se_mat <- sig2_hat * sand_mat

    ##############################
    # Derivatives

    dbetahat_dw <- solve(zwx, t(z * eps))

    # standard error partial derivatives
    dsig2_hat_dw_partial <- eps^2 / (num_obs - length(beta))
    dsig2_hat_dbeta <- -2 * colSums(w0 * eps * x) / (num_obs - length(beta))

    # Derivative of the diagonal of the sandwich matrix, which does not
    # depend on beta.
    # See notes for the definition of these terms and the derivation.
    AR_x <- zwx_inv_zwz %*% solve(t(zwx), t(x))
    R_z <- solve(zwx, t(z))
    dsand_mat_diag_dw_partial <- R_z * R_z - 2 * AR_x * R_z

    dse_mat_diag_dw_partial <-
        dsand_mat_diag_dw_partial * sig2_hat +
        outer(diag(sand_mat), dsig2_hat_dw_partial)

    # combine with the chain rule.
    # Note that dsig2_hat_dbeta gets broadcasted like a column vector.
    dsig2_hat_dw <- colSums(dsig2_hat_dbeta * dbetahat_dw) + dsig2_hat_dw_partial
    dse_mat_diag_dw <-
        dse_mat_diag_dw_partial +
        outer(diag(sand_mat), colSums(dsig2_hat_dbeta * dbetahat_dw))

    dse_dw <- 0.5 * dse_mat_diag_dw / sqrt(diag(se_mat))

    # Specify return values
    ret_list <- list(
        se_mat=se_mat,
        se=sqrt(diag(se_mat)),
        sig2_hat=sig2_hat,
        dbetahat_dw=dbetahat_dw,
        dsig2_hat_dw=dsig2_hat_dw,
        dse_mat_diag_dw=dse_mat_diag_dw,
        dse_dw=dse_dw)

    if (testing) {
        # For testing and debugging, it's useful to get the intermediate results.
        ret_list$betahat <- solve(zwx, t(z_w) %*% y) %>% as.numeric()
        ret_list$sand_mat <- sand_mat
        ret_list$dsig2_hat_dbeta <- dsig2_hat_dbeta # tested
        ret_list$dsig2_hat_dw_partial <- dsig2_hat_dw_partial # tested
        ret_list$dsand_mat_diag_dw_partial <- dsand_mat_diag_dw_partial # tested
        ret_list$dse_mat_diag_dw_partial <- dse_mat_diag_dw_partial # tested
    }

    return(ret_list)
}




ComputeIVRegressionInfluenceR <- function(iv_res, se_group=NULL) {
    x <- iv_res$x$regressors
    num_obs <- nrow(x)
    z <- iv_res$x$instruments
    y <- as.numeric(iv_res$y)
    betahat <- as.numeric(iv_res$coefficients)
    if (is.null(iv_res$weights)) {
        w0 <- rep(1.0, num_obs)
    } else {
        w0 <- iv_res$weights
    }


    iv_se_list <- GetIVSEDerivs(
      x=x, z=z, y=y, beta=betahat, w0=w0, se_group=se_group)


    # Note that the standard errors may not match iv_res when using se_group.
    return(list(model_fit=iv_res,
                n_obs=length(iv_res$y),
                regressor_names=colnames(iv_res$x$regressors),
                grad_fun="hand_coded_derivatives",

                betahat=betahat,
                se=iv_se_list$se_mat,
                weights=py_main$w0,

                beta_grad=py_main$betahat_grad,
                se_grad=py_main$se_grad)
    )
}
