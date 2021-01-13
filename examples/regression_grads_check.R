#base_dir <- Sys.getenv("REPO")
base_dir  <- "/home/rgiordan/Documents/git_repos/zaminfluence"
setwd(base_dir)
source(file.path(base_dir, "zaminfluence/tests/testthat/utils.R"))

library(ggplot2)
library(dplyr)
library(reshape2)
library(gridExtra)
library(sandwich)
library(zaminfluence)
library(AER)
library(numDeriv)

library(devtools)
devtools::load_all("zaminfluence")

py_main <- InitializePython(file.path(base_dir, "venv/bin/python3"))
compare <- function(x, y) { return(max(abs(x - y))) }
check_equivalent  <- function(x, y) { stopifnot(compare(x, y) < 1e-8) }

n_obs <- 10000

# The test utilities can simulate data.
source(file.path(base_dir, "zaminfluence/tests/testthat/utils.R"))

set.seed(42)


AssertNearlyZero <- function(x, tol=1e-15) {
    stopifnot(max(abs(x)) < tol)
}

#######################
# Let's do some checks

do_iv <- TRUE

x_dim <- 3
beta_true <- runif(x_dim)

if (do_iv) {
    df <- GenerateIVRegressionData(20, beta_true, num_groups=5)
} else {
    df <- GenerateRegressionData(20, beta_true, num_groups=5)
}
df$weights <- runif(nrow(df)) + 1

# Fit a model.

if (do_iv) {
    # IV:
    x_names <- sprintf("x%d", 1:x_dim)
    z_names <- sprintf("z%d", 1:x_dim)
    iv_form <- formula(sprintf("y ~ %s - 1 | %s - 1",
                               paste(x_names, collapse=" + "),
                               paste(z_names, collapse=" + ")))
    reg_fit <- ivreg(data=df, formula = iv_form, x=TRUE, y=TRUE, weights=weights)
    reg_cov <- sandwich::vcovCL(reg_fit, cluster=df$se_group, type="HC0", cadjust=FALSE)
} else {
    # Regression:
    x_names <- sprintf("x%d", 1:x_dim)
    reg_form <- formula(sprintf("y ~ %s - 1",
                                paste(x_names, collapse=" + ")))
    reg_fit <- lm(data=df, formula=reg_form, x=TRUE, y=TRUE, weights=weights)
    reg_cov <- sandwich::vcovCL(reg_fit, cluster=df$se_group, type="HC0", cadjust=FALSE)
}

#source(file.path(base_dir, "zaminfluence/R/ols_iv_grads_lib.R"))

# Get influence.

if (do_iv) {
    # iv_infl <- ComputeModelInfluence(reg_fit)
    # grad_df <- GetTargetRegressorGrads(iv_infl, "x1")
    # influence_dfs <- SortAndAccumulate(grad_df)
    reg_se_list <- GetIVSEDerivs(
        x=df[, x_names] %>% as.matrix(),
        z=df[, z_names] %>% as.matrix(),
        y=df$y,
        beta=reg_fit$coefficients,
        w0=df$weights,
        se_group=df$se_group,
        testing=TRUE)
} else {
    reg_se_list <- GetRegressionSEDerivs(
        x=df[, x_names] %>% as.matrix(),
        y=df$y,
        beta=reg_fit$coefficients,
        w0=df$weights,
        se_group=df$se_group,
        testing=TRUE)
}


###################################

LocalGetRegressionSEDerivs <- function(w=df$weights, beta=reg_fit$coefficients) {
    if (do_iv) {
        return(GetIVSEDerivs(
            x=df[, x_names] %>% as.matrix(),
            z=df[, z_names] %>% as.matrix(),
            y=df$y,
            beta=beta,
            w0=w,
            se_group=df$se_group,
            testing=TRUE))
    } else {
        return(GetRegressionSEDerivs(
            x=df[, x_names] %>% as.matrix(),
            y=df$y,
            beta=beta,
            w0=w,
            se_group=df$se_group,
            testing=TRUE))
    }
}



w0 <- df$weights
beta <- reg_fit$coefficients

AssertNearlyZero(reg_se_list$betahat - reg_fit$coefficients, tol=1e-11)
AssertNearlyZero(colMeans(reg_se_list$s_mat), tol=1e-11)
num_groups <- max(df$se_group) + 1
AssertNearlyZero(cov(reg_se_list$s_mat) * (num_groups - 1) / num_groups -
                 reg_se_list$v_mat, tol=1e-12)

vcov_se_cov <- vcovCL(reg_fit, cluster=df$se_group, type="HC0", cadjust=FALSE)
AssertNearlyZero(vcov_se_cov / num_groups - reg_se_list$se_mat, tol=1e-12)
AssertNearlyZero(sqrt(diag(vcov_se_cov)/ num_groups) - reg_se_list$se, tol=1e-12)

#######
# Check that aggregation is doing the right thing

# Passes
# Test that the s_mat_expanded worked correctly.
for (n in 1:length(df$se_group)) {
    gn <- df$se_group[n] + 1
    AssertNearlyZero(reg_se_list$s_mat_expanded[n, ] - reg_se_list$s_mat[gn, ], tol=1e-11)
}

#######

# Passes
ddiag_semat_dw_partial_num <-
    numDeriv::jacobian(function(w) {
        LocalGetRegressionSEDerivs(w=w)$se_mat %>% diag()
    }, w0)
AssertNearlyZero(ddiag_semat_dw_partial_num - reg_se_list$ddiag_semat_dw_partial, tol=1e-10)

# Passes
ddiag_semat_dbeta_partial_num <-
    numDeriv::jacobian(function(beta) {
        LocalGetRegressionSEDerivs(beta=beta)$se_mat %>% diag()
    }, beta)
AssertNearlyZero(ddiag_semat_dbeta_partial_num - reg_se_list$ddiag_semat_dbeta_partial, tol=1e-10)

#######################
# Full derivatives

GetRegTestResults <- function(w) {
    df_test <- df
    df_test$weights <- w
    beta_test <- LocalGetRegressionSEDerivs(w=w)$betahat
    reg_se_test_list <- LocalGetRegressionSEDerivs(beta=beta_test, w=w)
    return(
        list(beta=reg_se_test_list$betahat,
             se=reg_se_test_list$se,
             se_mat_diag=diag(reg_se_test_list$se_mat))
    )
}

# Passes
dbetahat_dw_num <-
    numDeriv::jacobian(function(w) { GetRegTestResults(w)$beta }, w0)
AssertNearlyZero(dbetahat_dw_num - reg_se_list$dbetahat_dw, tol=1e-10)

# Passes
dse_mat_diag_dw_num <-
    numDeriv::jacobian(function(w) { GetRegTestResults(w)$se_mat_diag }, w0)
AssertNearlyZero(dse_mat_diag_dw_num - reg_se_list$dse_mat_diag_dw, tol=1e-10)

# Passes
dse_dw_num <-
    numDeriv::jacobian(function(w) { GetRegTestResults(w)$se }, w0)
AssertNearlyZero(dse_dw_num - reg_se_list$dse_dw, tol=1e-10)


##############################
# OK!  What now?

