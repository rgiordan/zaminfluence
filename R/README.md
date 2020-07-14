The organization of the code is roughly as follows.

- `python_lib.R`: Functions for interfacing with Python.  Most of the actual
influence functions are estimated using Python's automatic differentiation
capacities.  This code mostly interfaces with the python package
currently called `regsens_rgiordandev` in the `regression_sensitivity`
folder of this repository.

- `sorting_lib.R`: Functions for taking the raw influence functions and sorting
them for particular kinds of adversarial perturbations.

- `paired_lib.R`: Similar to `sorting_lib.R`, but only removing pairs of
observations according to some category, typically treatment and control.

- `regression_sensitivity_lib.R`: Functions to compute diagnostics
from a given sorting (i.e. a given set of adversarial changes), for example
the number of observations to remove in order to change sign or significance.
Despite the name, these functions work generally and not only for regression.

- `plotting_lib.R`: Plotting functions.

- `rerun_lib.R`: Helper functions for re-running regression or IV regression
with sets of points left out.

- `error_bounds_lib.R`: Experimental library for computing the finite-sample
error bounds for regression.
