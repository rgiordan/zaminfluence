
# The main object is a Model Gradients object.  It is a list and must have
# - weights:    The original data weights
# - betahat:   The original estimator (length D)
# - sehat:      The original standard errors (length D)
# - n_obs:      The original number of observations
#
# - beta_grad:  A matrix (n_obs x D) of beta gradients
# - se_grad:    A matrix (n_obs x D) of se gradients
#
# For a regression, it must also have
# - regressor_names:    The names of the regressors
# - model_fit:          Everything you need to re-run stuff

# Then maybe a parameter inference object?  Need a better name.
# It has a:
# - target_index:   The index into betahat
# - sig_num_ses:    The number of ses that form a confidence interval
# - beta, beta_mzse, beta_pzse:  Influence summary objects (quantity of interest objects?)

# Influence summary objects contain an processed influence vector.
# (Quantity of interest objects?)
# They are the influence scores for a single quantity of interest.
# They are created with ProcessInfluenceVector().  Then contain
# - base_value:     The original value of the quantity of interest
# - neg, pos:       Sorted influence scores for the negative and positive
#                    influence scores, where.
# Sorted influence scores (neg and pos) have
# - infl_inds:      Indices into the original data that sort the influence
#                   scores of the corresponding sign.  E.g., infl_inds[1]
#                   for the `neg` entry is the index of the most negative
#                   influence score.
# - infl_cumsum:    The cumulative sum of the sorted influences scores with
#                   the specified sign.
# - num_obs:        The total number of observations (is this necessary?)
# - obs_per_row:    The number of observations per row (is this necessary?)



# Then there are signals, which are changes in quantities of interest
# that produce certain changes.  A signal can have
# - metric:         The name of the quantity of interest
# - signal:         The amount to change the quantity of interest
# - description:    A summary of what the change means
# - apip:           The APIP for this particular change
