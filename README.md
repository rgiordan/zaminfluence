# ZAM influence

This repo contains code to perform the analyses in our paper,
"Automatic Finite-Sample Robustness Metrics".
The repo name comes from "Z-estimator approximate maximal influence".

# Installation

This package has two components which need to work together: a Python component
and an R component.  So let's not kid ourselves: installation probably won't be
a cakewalk.

As of now, the easiest and safest way to install the package is by creating
a clone of the repository.  Change to your favorite directory for git
repos and run

```
git clone https://github.com/rgiordan/zaminfluence
cd zaminfluence
```

Optionally, set the `REPO` environment variable.  For example, for `bash` users,
run

```
export REPO=$(pwd)
```

For the rest of these instructions I'll use `$REPO` to refer to the path
to the directory with the git repository.

## Python part (`regression_sensitivity`)

Everything assumes you're using an updated version of Python 3.

1. Let's begin by creating a virtual environment.  In `$REPO`, run the following
commands:
```
python -m venv venv         # Create a virtual environment
source venv/bin/activate    # Activate the virtual environment.
```
Note that you need to activate the virtual environment in every new shell.

1. Install the package in your virtual environment.  Ideally, all dependencies
will be handled automatically, though we'll manually install the latest master
branches of `paragami` and `vittles`.
```
pip3 install --upgrade pip
pip3 install git+https://github.com/rgiordan/vittles.git@master
pip3 install git+https://github.com/rgiordan/paragami.git@master
pip3 install -e $REPO/regression_sensitivity
```

You can now run the python tests to make sure everything is working.

```
cd $REPO/inst/regression_sensitivity
python3 -m pytest
```

## R part (`zaminfluence`)

1. Install the R library from the ``zaminfluence`` directory.
```
library(devtools)
repo_loc <- Sys.getenv("REPO") # Or just set to the correct directory path
devtools::install_local(file.path(repo_loc, "zaminfluence"), force=TRUE)
```

1. In R, during each session where you want to use python, you must first run
```
library(zaminfluence)
repo_loc <- Sys.getenv("REPO") # Or just set to the correct directory path
InitializePython(file.path(repo_loc, "venv/bin/python"))
```
This will tell R to use your virtual environment's python.

## Done, possibly!

You should now be able to use the functions on your own regressions. You can
look at the `examples/regression_intuition/intuition.ipynb` notebook and
associated library, `intuition_2_lib.R`, for examples. Some paltry documentation
can be found in the function comments. You'll want to first run your regression,
then call `GetInfluence` first, then `AnalyzeInfluence`, and then the other
functions as necessary.
