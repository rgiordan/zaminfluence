#base_dir <- Sys.getenv("REPO")
base_dir  <- "/home/rgiordan/Documents/git_repos/zaminfluence"
setwd(base_dir)

library(ggplot2)
library(dplyr)
library(reshape2)
library(gridExtra)
library(sandwich)
library(zaminfluence)
library(AER)

py_main <- InitializePython(file.path(base_dir, "venv/bin/python3"))

reticulate::use_python(file.path(base_dir, "venv/bin/python3"))
py_main <- reticulate::import_main()

reticulate::py_run_string("import regsens_rgiordandev")

reticulate::py_run_string("x = paragami.__file__")
py_main$x

reticulate::py_run_string("x = paragami.__version__")
py_main$x

compare <- function(x, y) { return(max(abs(x - y))) }
check_equivalent  <- function(x, y) { stopifnot(compare(x, y) < 1e-8) }

n_obs <- 10000

# The test utilities can simulate data.
source(file.path(base_dir, "zaminfluence/tests/testthat/utils.R"))

set.seed(42)
