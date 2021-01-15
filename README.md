# ZAM influence

This repo contains code to perform the analyses in our paper,
"Automatic Finite-Sample Robustness Metrics".
The repo name comes from "Z-estimator approximate maximal influence".

# Installation

1. You can install the R library directly from github.
```
> library(devtools)
> devtools::install_github("https://github.com/rgiordan/zaminfluence/",
                           ref="reg_derivs",
                           upgrade="never",
                           force=TRUE)
```

You can install different branches using the `ref` argument.

2. Run the R tests to make sure everything is working.
```
> devtools::test("zaminfluence")
```

3. Done, hopefully!  You should now be able to run the script in
   `examples/simple_examples.R`.

Please submit an issue or email us if you have any questions or comments!

For defining sensitivity to custom objectives, you may want to use
Python automatic differentiation tools.  To do that, follow the
(more complicated) [Python installation instructions](python_installation.md).
