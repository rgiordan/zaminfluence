
# Overall organization

The code is organized into a hierarchy of objects, which I try to refer
to using common variable names.

TODO: Define rerun objects, too.


- **Model Gradients** objects (`model_grads`).  They contain information for an
entire model, including all the information needed to re-run it, and all the
gradients that will be needed to compute influence functions.
A Model Gradients object will typically
contain a number of Parameter Influence objects.
- **Parameter Influence**  objects (`param_infl`).  A parameter influence object
stores a number of quantities of interest and their influence functions for
a particular parameter in the model.
A Parameter Influence object will contain a number of
Quantity of Interest objects.  Currently library, the three
quantities of interest stored for each parameter are its point estimate
and the two ends of a confidence interval.
- **Quantity of Interest (QOI)** objects (`qoi`).  This contains
information about a
particular quantity of interest, including sorted influence functions and
the base value.  The key quantites of our paper --- the approximate
perturbation-inducing proportion (APIP), approximate maximally influential
set (AMIS), and approximate maximumally influential perturbation (AMIP) ---
are computed for and from Quantity of Interest objects.
- **Signal** objects (`signal`).  A signal object stores information about a
particular change in a quantity of interest, including a description of
what the change means, the approximate
perturbation-inducing proportion (APIP), the corresponding approximate maximally
influential set (AMIS), and, optionally, the result of a re-run.  A signal
only uses a single Quantity of Interest, but a Quantity of Interest may
appear in more than one signal, if changes of different magnitudes have
different meanings.

The code has an organization that corresponds to the above structures.

- Model Gradient objects for regression are computed in `ols_iv_grads_lib.R`.
Most of `ols_iv_grads_lib.R` is for computing the acutal influence scores,
the actual creation of the Model Gradient objects is in a couple thin wrappers
at the end.
- QOI objects are computed and wrangled in `influence_lib.R`.
- Parameter influence and signal objects dealing with changes of sign,
significance, and both sign and significance are constructed and wrangled in
`inference_targets_lib.R`.

# Data Structure Details

### Model Gradients

The most capacious object is a Model Gradients object.  They are
created with `ComputeModelInfluence`.  A model is a fit
of a `D`-dimensional parameter using `n_obs` datapoints  It is a list and must
have the following fields:
- `weights`:    The original data weights (length `n_obs`)
- `betahat`:   The original estimator (length `D`)
- `se`:      The original standard errors (length `D`)
- `n_obs`:      The original number of observations
- `beta_grad`:  A matrix (`n_obs` x `D`) of gradients of the point estimates
- `se_grad`:    A matrix (`n_obs` x `D`) of gradients of the standard errors
- `regressor_names`:    The names of the regressors
- `model_fit`:          Everything you need to re-fit the model
- `param_infl_list`:  An (optional) list of Parameter Influence objects.
- `RerunFun`:  A function taking the arguments `model_fit` and a vector of
boolean weights and returning a refit object.

### Parameter Influence

A Parameter Influence object contains influence scores for inference concerning
a single parameter form a Model Gradients object.  They are
created with `AppendTargetRegressorInfluence`.
It is a list and must have the following fields:
- `target_index`:   The index into betahat
- `sig_num_ses`:    The number of ses that form a confidence interval
- `beta`:  A QOI object for the point estimate
- `beta_mzse`:  A QOI object for the lower end of a confidence
interval ("minus z standard errors").
- `beta_pzse`:  A QOI object for the upper end of a confidence
interval ("plus z standard errors").

### Quantity of Interest

QOI objects contain a processed influence vector for a particular scalar-valued
quantity of interest.
They are created with `ProcessInfluenceVector`.  Then contain
- `base_value`:     The original value of the quantity of interest
- `infl`: The unsorted influence scores (in the same order as the original data)
- `num_obs`:        The total number of observations (is this necessary?)
- `obs_per_row`:    The number of observations per row (is this necessary?)
- `neg`, `pos`:       Sorted influence scores for the negative and positive
influence scores, where sorted influence scores have:
  - `infl_inds`:      Indices into the original data that sort the influence
scores of the corresponding sign.  For example, `infl_inds[1]`
for the `neg` entry is the index in the order of the original
data of the most negative influence score.  Equivalently,
`infl[neg$infl_inds[1]] == min(infl)`, and
`infl[pos$infl_inds[1]] == max(infl)`.
  - `infl_cumsum`:    The cumulative sum of the sorted influences scores with
the specified sign.

QOI objects can be maniupated using functions in `influence_lib.R`.

### Signal

A signal records a target change in a QOI.  A signal is a list and must have
- `qoi_name`:         The name of the quantity of interest
- `signal`:         The amount to change the quantity of interest
- `description`:    A plain language description of what the change means
- `apip`:           The APIP for this particular change, which contains
    - `n`: The number of points to remove
    - `prop`: The proportion of points to remove
    - `inds`: The indices (in the original data order) of the data dropped in the corresponding AMIS.

Currently, signals are created with `GetInferenceSignals`, which defines signals
to change the sign, significance, and both sign and significance.

Optionally, a signal may also contain:
- `rerun`: The complete result of a re-run at the AMIS.
- `rerun_df`: A tidy summary dataframe of the predicted and actual changes
in the parameter and its confidence interval.


# TODOs

- Rename beta to theta everywhere
- Consistently use "hat" or not in the names (I think no hat)
- Document rerun objects
- Consistently use "num" or "n" everywhere (I think num)
