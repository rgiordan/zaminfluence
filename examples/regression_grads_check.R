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


###################

do_iv <- FALSE

x_dim <- 3
beta_true <- runif(x_dim)

if (do_iv) {
    df <- GenerateIVRegressionData(20, beta_true, num_groups=5)
} else {
    df <- GenerateRegressionData(20, beta_true, num_groups=5)
}
df$weights <- runif(nrow(df)) + 1
df$se_group <- 1:nrow(df) - 1
df$se_group <- NULL

# Fit a model.
if (do_iv) {
    # IV:
    x_names <- sprintf("x%d", 1:x_dim)
    z_names <- sprintf("z%d", 1:x_dim)
    iv_form <- formula(sprintf("y ~ %s - 1 | %s - 1",
                               paste(x_names, collapse=" + "),
                               paste(z_names, collapse=" + ")))
    reg_fit <- ivreg(data=df, formula = iv_form,
                     x=TRUE, y=TRUE, weights=weights)
    # reg_cov <- sandwich::vcovCL(
    #   reg_fit, cluster=df$se_group, type="HC0", cadjust=FALSE)
} else {
    # Regression:
    x_names <- sprintf("x%d", 1:x_dim)
    reg_form <- formula(sprintf("y ~ %s - 1",
                                paste(x_names, collapse=" + ")))
    reg_fit <- lm(data=df, formula=reg_form,
                  x=TRUE, y=TRUE, weights=weights)
}


LocalGetRegressionSEDerivs <- function(w=df$weights,
                                       beta=reg_fit$coefficients) {
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


reg_se_list <- LocalGetRegressionSEDerivs()
reg <- tidy(reg_fit)
reg_cov <- sandwich::vcovCL(reg_fit, type="HC0", cadjust=FALSE)

sqrt(diag(reg_cov))
reg$std.error
sqrt(diag(reg_se_list$se_mat))



#######################
# Let's do some checks

iv_infl <- ComputeModelInfluence(reg_fit)
grad_df <- GetTargetRegressorGrads(iv_infl, "x1")
influence_dfs <- SortAndAccumulate(grad_df)
