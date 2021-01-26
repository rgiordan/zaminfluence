# ZAM influence

This repo contains code to perform analyses described in our paper, "Automatic
Finite-Sample Robustness Metrics," for the case of ordinary least squares and
instrumental variable regressions. The repo name comes from "Z-estimator
approximate maximal influence".

# Installation

1. You can install the R library directly from github.
```
> library(devtools)
> devtools::install_github("https://github.com/rgiordan/zaminfluence/",
                           ref="master",
                           force=TRUE)
```

You can install different branches using the `ref` argument.

2. Done, hopefully!  You should now be able to run the script in
   `examples/simple_examples.R`.

Please submit an issue or email us if you have any questions or comments!

# Python backend

The original version of `zaminfluence` depended on a Python backend.  The
Python part is no longer required.

However, if you know what you're doing and are interested in defining
sensitivity to custom objectives, you may want to use Python automatic
differentiation tools. To do that, follow the (more complicated) [Python
installation instructions](python_installation.md).  Currently the interface for
custom objectives is not documented, but the authors could be prompted to write
such documentation if there is a need.  In other words, send us an email!
