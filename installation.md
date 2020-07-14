# AdversarialInfluenceWorkbench

This private repo contains our work on adversarial robustness.

This regression sensitivity package has two components: a Python component
and an R component.  So don't expect installation to be a cakewalk.

## Python part

Everything assumes you're using an updated version of Python 3.  Let's
assume that
`$REPO` is the path to the `AdversarialInfluenceWorkbench/`
directory in your cloned git repo.

1. Let's begin by creating a virtual environment.  In `$REPO`, run the following
commands:
```
python -m venv venv         # Create a virtual environment
source venv/bin/activate    # Activate the virtual environment.
```
Note that you need to activate the virtual environment in every new shell.

1. Install the package in your virtual environment.  Ideally, all dependencies
will be handled automatically, though we'll manually install the latest master
branches of `paragami` and `vittles`.  Ask Ryan if there's a problem.
```
pip3 install --upgrade pip
pip3 install git+https://github.com/rgiordan/vittles.git@master
pip3 install git+https://github.com/rgiordan/paragami.git@master
pip3 install -e $REPO/regression_sensitivity
```

1. Optional: create an ipython kernel that uses the virtual environment.
```
python -m ipykernel install --user --name=regsens_rgiordandev
```

## R part


1. Install the R library from the ``rminfluence`` directory.
```
library(devtools)
devtools::install_local("$REPO/rminfluence")
```
Of course, you need to manually replace `$REPO` with the actual path.

1. In R, during each session where you want to use python, you must first run
```
library(rminfluence)
InitializePython($REPO/venv/bin/python)
```
Of course, you need to manually replace `$REPO` with the actual path.  This will
tell R to use your virtual environment's python.

## Done, possibly!

You should now be able to use the functions on your own regressions. You can
look at the `examples/regression_intuition/intuition.ipynb` notebook and
associated library, `intuition_2_lib.R`, for examples. Some paltry documentation
can be found in the function comments. You'll want to first run your regression,
then call `GetInfluence` first, then `AnalyzeInfluence`, and then the other
functions as necessary.
