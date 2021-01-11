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


####################
# Sanity checks

y <- reg_fit$y
x <- reg_fit$x$regressors
z <- reg_fit$x$instruments


num_obs <- length(y)

eps <- as.numeric(y - x %*% beta)

z_w <- z * w0
z_eps <- z * eps
z_w_eps <- z * eps * w0

zwz <- t(z_w) %*% z
zwx <- t(z_w) %*% x


# grouped aggregation
z <- reg_fit$x$instruments
group_rows <- split(1:nrow(z), df$se_group)
s_mat <- lapply(group_rows,
       function(rows) { colSums(z_w[rows, , drop=FALSE] * eps[rows])}) %>%
    do.call(rbind, .)
s_mat <- s_mat - colMeans(s_mat)
# Check for well-formedness of the groups.
all(as.numeric(row.names(s_mat))  == unique(df$se_group) %>% sort())
all(as.numeric(row.names(s_mat))  == (min(df$se_group):max(df$se_group)))
v_mat <- t(s_mat) %*% s_mat / nrow(s_mat)

s_mat_expanded <- matrix(NA, dim(z)[1], dim(z)[2])
for (row in names(group_rows)) {
    s_mat_expanded[group_rows[[row]], ] <- s_mat[as.numeric(row) + 1]
}

##################
# New code

source(file.path(base_dir, "zaminfluence/R/ols_iv_grads_lib.R"))

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
        w0=df$weights, testing=TRUE)
}

reg_se_list
df$se_group

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



######

ddiag_semat_dw_partial_num <-
    numDeriv::jacobian(function(w) {
        LocalGetRegressionSEDerivs(w=w)$se_mat %>% diag()
    }, w0)
    
reg_se_list$ddiag_semat_dw_partial
ddiag_semat_dw_partial_num

#AssertNearlyZero(ddiag_vmat_dw_partial_num - reg_se_list$ddiag_semat_dw_partial, tol=1e-8)
plot(ddiag_semat_dw_partial_num, reg_se_list$ddiag_semat_dw_partial); abline(0, 1)







######### old:

###########################
# Load and test

AssertNearlyZero(reg_se_list$sig2_hat - sigma(reg_fit)^2)
AssertNearlyZero(reg_se_list$se_mat - vcov(reg_fit), tol=1e-11)
AssertNearlyZero(reg_se_list$se - vcov(reg_fit) %>% diag() %>% sqrt(), tol=1e-11)


# Test the partial derivatives

w0 <- df$weights
beta <- reg_fit$coefficients

# Passes
dsig2_hat_dw_num <-
    numDeriv::jacobian(function(w) {
        LocalGetRegressionSEDerivs(w=w)$sig2_hat
    }, w0) %>%
    as.numeric()
AssertNearlyZero(dsig2_hat_dw_num - reg_se_list$dsig2_hat_dw_partial, tol=1e-8)

# Passes
dsig2_hat_dbeta_num <-
    numDeriv::jacobian(function(beta) {
        LocalGetRegressionSEDerivs(beta=beta)$sig2_hat
    }, beta)
AssertNearlyZero(dsig2_hat_dbeta_num - reg_se_list$dsig2_hat_dbeta, tol=1e-8)

# Passes
dsand_mat_diag_dw_num <-
    numDeriv::jacobian(function(w) {
        LocalGetRegressionSEDerivs(w=w)$sand_mat %>% diag()
    }, w0)
AssertNearlyZero(dsand_mat_diag_dw_num - reg_se_list$dsand_mat_diag_dw_partial, tol=1e-6)

# Passes
dse_mat_diag_dw_num <-
    numDeriv::jacobian(function(w) {
        LocalGetRegressionSEDerivs(w=w)$se_mat %>% diag()
    }, w0)
AssertNearlyZero(dse_mat_diag_dw_num - reg_se_list$dse_mat_diag_dw_partial, tol=5e-6)





#############################
# Test the full derivatives

GetRegTestResults <- function(w) {
    df_test <- df
    df_test$weights <- w
    beta_test <- lm(data=df_test, formula=reg_form,
                       x=TRUE, y=TRUE, weights=weights)$coefficients
    reg_se_test_list <- LocalGetRegressionSEDerivs(beta=beta_test, w=w)
    return(
        list(beta=reg_se_test_list$betahat,
             se=reg_se_test_list$se,
             sig2_hat=reg_se_test_list$sig2_hat,
             se_mat_diag=diag(reg_se_test_list$se_mat))
    )
}


dbetahat_dw_num <-
    numDeriv::jacobian(function(w) { GetRegTestResults(w)$beta }, w0)
AssertNearlyZero(dbetahat_dw_num - reg_se_list$dbetahat_dw, tol=1e-8)
#plot(dbetahat_dw_num, iv_se_list$dbetahat_dw); abline(0, 1)

dsig2_hat_num <-
    numDeriv::jacobian(function(w) { GetRegTestResults(w)$sig2_hat }, w0)
#plot(dsig2_hat_num, iv_se_list$dsig2_hat_dw); abline(0, 1)
AssertNearlyZero(dsig2_hat_num - reg_se_list$dsig2_hat_dw, tol=1e-6)


dse_mat_diag_dw_num <-
    numDeriv::jacobian(function(w) { GetRegTestResults(w)$se_mat_diag }, w0)
#plot(dse_mat_diag_dw_num, iv_se_list$dse_mat_diag_dw); abline(0, 1)

# Check the relative error for this one
AssertNearlyZero((dse_mat_diag_dw_num - reg_se_list$dse_mat_diag_dw) /
                     (reg_se_list$dse_mat_diag_dw + 1e-3), tol=1e-6)

dse_dw_num <- numDeriv::jacobian(function(w) { GetRegTestResults(w)$se }, w0)
#plot(dse_dw_num, iv_se_list$dse_dw); abline(0, 1)
# Check the relative error for this one
AssertNearlyZero((dse_dw_num - reg_se_list$dse_dw) / (reg_se_list$dse_dw + 1e-3), tol=1e-6)

