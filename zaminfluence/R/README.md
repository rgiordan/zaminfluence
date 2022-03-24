
# Overall organization

The code is organized into a hierarchy of objects, which I try to refer
to using common variable names, and which each have their own S3 classes.

- **Model Fit** objects
(variable `model_fit`, S3 class `ModelFit`,
see `model_grads_lib.R`).
The information for a single fit of a model, including parameter and
standard error estimates, as well as the information
needed to re-run it at different weights.

- **Model Gradients** objects
(variable `model_grads`, S3 class `ModelGrads`, see `model_grads_lib.R`).
All the gradients for model evaluated at a particular fit.
A Model Gradients object will typically
contain a number of Parameter Influence objects.

- **Parameter Inference Influence**  objects
(variable `param_infl`, S3 class `ParameterInferenceInfluence`,
see `inference_targets_lib.R`).
A number of quantities of interest and their
influence functions for performing frequentist inference on a particular
parameter in the model. A Parameter Influence object will contain a number of
Quantity of Interest objects --- namely, for the point estimate, the
standard error, and the two ends of a
confidence interval.

- **Quantity of Interest (QOI) influence** objects
(variable `qoi`, S3 class `QOIInfluence`, see `influence_lib`).
Information about a
particular quantity of interest, including sorted influence functions and
the base value.  The key quantites of our paper --- the approximate
perturbation-inducing proportion (APIP), approximate maximally influential
set (AMIS), and approximate maximumally influential perturbation (AMIP) ---
are computed for and from `QOIInfluence` objects.

- **QOI Signal** objects
(variable `signal`, S3 class `QOISignal`, see `inference_targets_lib.R`).
Information about a particular change in a quantity of interest,
including a description of what the change means, the approximate
perturbation-inducing proportion
(APIP), the corresponding approximate maximally influential set (AMIS).

- **Approximate Maximum Influence Perturbation** objects
(variable `apip`, S3 class `APIP`, see `influence_lib.R`).  An estimate,
based on the influence function, of the data points to drop to effect a
particular change in a particular quantity of interest.

- **Parameter and Target Change lists**.  Though not an S3 class, the
function `GetInferenceSignals` computes a nested list of QOI Signals
for all parameters and target changes in a particular Model Grads object.
The outer index is into parameter names, and the inner index is target
changes.  The functions `RerunForSignals` and `PredictForSignals` take
in a `signals` list so structured and return lists of `ModelFit` objects
with the same structure, which can then be passed to
`GetSignalsAndRerunsDataframe` to get a tidy summary of refits and predictions.


# Data Structure Details

### ModelFit

A model is a fit of a `D`-dimensional parameter using `num_obs` datapoints.
- `fit_object`:          Everything you need to re-fit the model
- `num_obs`:      The original number of observations
- `parameter_names`:    The names of the parameters
- `param`:   The estimator of the parameter (length `D`)
- `se`:      The estimate of the standard errors (length `D`)
- `weights`:    The original data weights (length `num_obs`)
- `parameter_dim`:    The dimension `D`
- `se_group`: If non-null, a vector of length `num_obs` used for grouped
standard errors.

### ModelGrads

A `ModelFit` together with gradients of the parameter and standard errors
with respect to the data weights.
- `model_fit`:  A `ModelFit` at the original weights.
- `parameter_names`:  The names of parameters for which gradients are computed.
- TODO: `D` should be replaced for all grads by select parameters,
and the columns should be named.
- TODO: is it actually num_obs x D or D x num_obs?
- `param_grad`:  A matrix (`num_obs` x `D`) of gradients of the point estimates
- `se_grad`:    A matrix (`num_obs` x `D`) of gradients of the standard errors
- `param_infls`:  An (optional) list of `ParameterInferenceInfluence` objects.
- `RerunFun`:  A function taking new weights and returning a new `ModelFit`
object evaluated at the new weights.

### ParameterInferenceInfluence

A Parameter Influence object contains influence scores for inference concerning
a single parameter form a Model Gradients object.  They are
created with `AppendTargetRegressorInfluence` and are typically stored in
the `param_infls` field of a `ModelGrads` object.
- TODO: Can we get rid of target_index and get parameters by name?
- `target_index`:   The index into param for this parameter
- `target_parameter`:   The  parameter name (as in `ModelFit$parameter_names`)
- `sig_num_ses`:    The number of ses that form a confidence interval
- `se`: A `QOIInfluence` object for the standard error
- `param`:  A  `QOIInfluence` object for the point estimate
- `param_mzse`:  A  `QOIInfluence` object for the lower end of a confidence
interval ("minus z standard errors").
- `param_pzse`:  A  `QOIInfluence` object for the upper end of a confidence
interval ("plus z standard errors").

### QOIInfluence

QOI objects contain a processed influence vector for a particular scalar-valued
quantity of interest.
- `base_value`:     The original value of the quantity of interest
- `infl`: The unsorted influence scores (in the same order as the original data)
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
  - `num_obs`:   The total number of observations in the original dataset

### QOISignal

A signal records a target change in a QOI.  A signal is a list and must have
- `qoi`:         The target `QOIInfluence` object
- `signal`:         The numerical amount to change the quantity of interest
- `description`:    A plain language description of what the change means
- `apip`:    An `APIP` object for this particular change.

### APIP
The APIP for a particular change, which contains
- `n`: The number of points to remove
- `prop`: The proportion of points to remove
- `inds`: The indices (in the original data order) of the data dropped in
the corresponding Approximate Maximally Influential Set (AMIS).


# Extending `zaminfluence`

Currently, we have implemented OLS and IV regression.
In order to use `zaminfluence` for your own models, you need only to
provide your own instances of the `ModelFit` and `ModelGrads` classes,
as instantiated by the corresponding S3 functions `ModelGrads` and `ModelFit`.
We recommend using automatic differentiation in Python together with a
wrapper in `reticulate` that translates your results into valid
`ModelGrads` and `ModelFit` objects.
