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


#######################
# Let's do some checks


x_dim <- 3
beta_true <- runif(x_dim)
#df <- GenerateIVRegressionData(10, beta_true, num_groups=NULL)
df <- GenerateRegressionData(10, beta_true, num_groups=NULL)
df$weights <- runif(nrow(df)) + 1

# Fit a model.
# x_names <- sprintf("x%d", 1:x_dim)
# z_names <- sprintf("z%d", 1:x_dim)
# iv_form <- formula(sprintf("y ~ %s - 1 | %s - 1",
#                            paste(x_names, collapse=" + "),
#                            paste(z_names, collapse=" + ")))
# iv_fit <- ivreg(data = df, formula = iv_form, x=TRUE, y=TRUE, weights=weights)
x_names <- sprintf("x%d", 1:x_dim)
reg_form <- formula(sprintf("y ~ %s - 1",
                           paste(x_names, collapse=" + ")))
reg_fit <- lm(data=df, formula=reg_form, x=TRUE, y=TRUE, weights=weights)

# Get influence.
# iv_infl <- ComputeModelInfluence(iv_fit)
# grad_df <- GetTargetRegressorGrads(iv_infl, "x1")
# influence_dfs <- SortAndAccumulate(grad_df)

reg_se_list <- GetRegressionSEDerivs(
    x=df[, x_names] %>% as.matrix(),
    y=df$y,
    beta=reg_fit$coefficients,
    w0=df$weights, testing=TRUE)



##################
# New code

AssertNearlyZero <- function(x, tol=1e-15) {
    stopifnot(max(abs(x)) < tol)
}


###########################
# Load and test

AssertNearlyZero(reg_se_list$betahat - reg_fit$coefficients, tol=1e-11)
AssertNearlyZero(reg_se_list$sig2_hat - sigma(reg_fit)^2)
AssertNearlyZero(reg_se_list$se_mat - vcov(reg_fit), tol=1e-11)
AssertNearlyZero(reg_se_list$se - vcov(reg_fit) %>% diag() %>% sqrt(), tol=1e-11)


