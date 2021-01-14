# ZAM influence

This repo contains code to perform the analyses in our paper,
"Automatic Finite-Sample Robustness Metrics".
The repo name comes from "Z-estimator approximate maximal influence".

# Installation

1. You can install the R library directly from github.
```
library(devtools)
devtools::install_github("https://github.com/rgiordan/zaminfluence/",
                         ref="reg_derivs",
                         force=TRUE)
```

2. In R, during each session where you want to use python, you must first run
```
library(zaminfluence)
repo_loc <- Sys.getenv("REPO") # Or just set to the correct directory path
InitializePython(file.path(repo_loc, "venv/bin/python"))
```
This will tell R to use your virtual environment's python, where you've
installed the necessary python libraries.  You can then
use the `zaminfluence` functions.

3. You can now run the R tests to make sure everything is working.
```
cd $REPO/zaminfluence/tests
./testthat.R
```

## Done, hopefully!

You should now be able to run the script in `examples/simple_examples.R`.

Please submit an issue or email us if you have any questions or comments!
