import autograd
import autograd.numpy as np
import scipy as sp

from copy import deepcopy

import paragami
import vittles

import time


def _validate(y=y, x=x):
    n_obs = x.shape[0]
    x_dim = x.shape[1]
    if len(y) != n_obs:
        raise ValueError(
            'The length of ``y`` must match the number of rows in ``x``.')
    return n_obs, x_dim

def reg(y=y, x=x,
        w=None,
        offset=None):
    """The regression parameter at the given perturbations.  This should
    be the optimum of ``reg_obj`` with the corresponding parameters.
    """

    n_obs, x_dim = _validate(y=y, x=x)
    if w is not None:
        x_w = x * np.expand_dims(w, axis=1)
    else:
        x_w = x
    if offset is None:
        offset = np.zeros(x_dim)
    x_wt_x = x_w.T @ x
    x_wt_y = x_w.T @ y
    return np.linalg.solve(x_wt_x, x_wt_y - offset).flatten()


def reg_obj(beta, y=y, x=x, w=None, offset=None):
    """The objective function for linear regression.
    """
    n_obs, x_dim = _validate(y=y, x=x)
    resid = y - (x @ beta)
    x_resid = (x * np.expand_dims(w, axis=1)).T @ resid - offset
    result = x_resid.T @ x_resid
    assert result.shape == ()
    # For numerical stability, it is important to include n_obs in the objective.
    return np.sum(result) / n_obs


# Estimate how many datapoints we would have to remove to effect a change of delta.
def inds_to_effect_change(leverage, desired_delta):
    # Argsort sorts low to high.
    # We are removing points, so multiply by -1.
    sort_inds = np.argsort(leverage * np.sign(desired_delta))
    deltas = -1 * np.cumsum(leverage[sort_inds])
    change_sign_inds = np.argwhere(
        np.sign(desired_delta) * (desired_delta - deltas) <= 0.)
    if len(change_sign_inds) > 0:
        first_ind_change_sign = np.min(change_sign_inds)
        remove_inds = sort_inds[:(first_ind_change_sign + 1)]
        return remove_inds
    else:
        return None


def print_change_results(inds, effect_str, lev_len):
    print('Assuming linearity, which may not hold for large numbers of points,')
    if inds is not None:
        print('removing {} observations ({:0.2f}%) would {}.'.format(
               len(inds), 100 * len(inds) / lev_len, effect_str))
    else:
        print('no number of observations would {}.'.format(effect_str))
