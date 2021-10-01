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

library(Matrix)

library(devtools)
devtools::load_all("zaminfluence")

py_main <- InitializePython(file.path(base_dir, "venv/bin/python3"))
compare <- function(x, y) { return(max(abs(x - y))) }
check_equivalent  <- function(x, y) { stopifnot(compare(x, y) < 1e-8) }

num_obs <- 10000

# The test utilities can simulate data.
setwd(file.path(base_dir, "zaminfluence/tests/testthat/"))
source(file.path(base_dir, "zaminfluence/tests/testthat/utils.R"))
source(file.path(base_dir, "zaminfluence/tests/testthat/test_derivs.R"))


debug(TestGroupedRegressionDerivatives)
TestGroupedRegressionDerivatives(do_iv = FALSE)



ExpandGroupedSum <- function(s_mat, group) {
    group_rows <- split(1:length(group), group)
    s_mat_expanded <- matrix(NA, length(group), ncol(s_mat))
    for (g in names(group_rows)) {
        for (row in group_rows[[g]]) {
            s_mat_expanded[row, ] <- s_mat[as.numeric(g) + 1, ]
        }
    }
    return(s_mat_expanded)
}

s_mat <- 1:15 %>% matrix(5, 3)
group <- rep(1:5, each=3) - 1

ExpandGroupedSum(s_mat, group) -  s_mat[group + 1, ]

#####################



set.seed(42)


AssertNearlyZero <- function(x, tol=1e-9) {
    x_norm <- max(abs(x))
    info_str <- sprintf("%e > %e", x_norm, tol)
    expect_true(x_norm < tol, info=info_str)
}


###############


###################

do_iv <- FALSE

x_dim <- 3
param_true <- runif(x_dim)
num_groups <- NULL

if (do_iv) {
    df <- GenerateIVRegressionData(20, param_true, num_groups=num_groups)
} else {
    df <- GenerateRegressionData(20, param_true, num_groups=num_groups)
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


reg_fit



#######################
# Let's do some checks

reg_infl <- ComputeModelInfluence(reg_fit)
grad_df <- GetTargetRegressorGrads(iv_infl, "x1")
influence_dfs <- SortAndAccumulate(grad_df)

if (do_iv) {
    reg_infl_python <- ComputeIVRegressionInfluencePython(reg_fit)
} else {
    reg_infl_python <- ComputeRegressionInfluencePython(reg_fit)
}

AssertNearlyZero(reg_infl_python$param - reg_infl$param)
AssertNearlyZero(reg_infl_python$param_grad - reg_infl$param_grad)
AssertNearlyZero(reg_infl_python$se_grad - reg_infl$se_grad)


