import autograd
import autograd.numpy as np
import scipy as sp

from copy import deepcopy

import paragami
from paragami.autograd_supplement_lib import grouped_sum
import vittles

import time


def _validate(y, x):
    n_obs = x.shape[0]
    x_dim = x.shape[1]
    if len(y) != n_obs:
        raise ValueError(
            'The length of ``y`` must match the number of rows in ``x``.')
    return n_obs, x_dim


def reg(y, x,
        w=None,
        offset=None):
    """The regression parameter at the given perturbations.  This should
    be the optimum of ``reg_obj`` with the corresponding parameters.
    """

    n_obs, x_dim = _validate(y=y, x=x)

    if offset is None:
        offset = np.zeros(x_dim)

    if w is not None:
        x_w = x * np.expand_dims(w, axis=1)
    else:
        x_w = x

    x_wt_x = x_w.T @ x / n_obs
    x_wt_y = x_w.T @ y / n_obs
    return np.linalg.solve(x_wt_x, x_wt_y - offset).flatten()


def reg_obj(beta, y, x, w=None, offset=None):
    """The objective function for linear regression.

    Here, I use the weighted method-of-moments objective function defined as
    follows.  Let epsilon = Y - X beta, and let X^T epsilon = m.
    The objective function is to set m equal to offset using a weighted
    quadratic loss, which is

    (m - offset)^T (X^T X)^{-1} (m - offset)

    The weighting matrix (X^T X)^{-1} improves numerical stability and gives
    a loss function that is identical to OLS when offset = 0.
    """
    n_obs, x_dim = _validate(y=y, x=x)

    if offset is None:
        offset = np.zeros_like(beta)

    if w is not None:
        xw = x * w[:, None]
        xtx = xw.T @ x / n_obs
        xty = xw.T @ y / n_obs
    else:
        xtx = x.T @ x / n_obs
        xty = x.T @ y / n_obs

    # This is the method of moments objective after expanding, dropping
    # terms that do not depend on beta, and collecting.
    result = \
        np.dot(beta, xtx @ beta) + \
        2 * np.dot(beta, offset - xty)

    assert result.shape == ()
    return result


def get_standard_error_matrix(betahat, y, x, w, se_group=None):
    """Return the standard error matrix for the regression estimate betahat.

    If se_group is None, compute the ordinary regression standard error.
    Otherwise, compute the robust standard errors using the grouping given by
    se_group, which is assumed to be integers 0:(num_groups - 1).

    Note that se_group must be zero-indexed, and the number of groups is taken
    to be the largest index plus one.  (This behavior is implicitly assumed in
    group_sum.)

    With the se_group option, no finite-sample bias adjustment is applied.
    For example, the resulting ses should be equivalent to calling the
    R function

    sandwich::vcovCL(..., cluster=se_group, type="HC0", cadjust=FALSE)
    """

    resid = y - x @ betahat

    # For now, I am taking the weights to parameterize a change to the
    # objective function rather than a change to the empirical distribution.
    # See email from me to Rachael and Tamara on Jan 31, 2020, 2:50 PM
    # for more discussion of this subtle point.
    if se_group is None:
        # I am using num_obs instead of np.sum(w) because w does not
        # parameterize the empirical distribution.
        num_obs = len(y)
        xtx_bar = np.einsum('ni,nj,n->ij', x, x, w) / num_obs
        sigma2hat = np.sum(w * (resid ** 2)) / (num_obs - len(betahat))
        xtx_inv = np.linalg.inv(xtx_bar)
        se2 = sigma2hat * xtx_inv / num_obs
        return se2
    else:

        if len(se_group) != len(y):
            raise ValueError("se_group must be the same length as the data.")
        #resid =  y - x @ betahat

        if np.min(se_group) != 0:
            raise ValueError('se_group must be zero-indexed ' +
                             '(its minimum must be zero)')

        # Calculate the sample variance of the gradient where each group
        # is treated as a single observation.
        grad = w[:, None] * resid[:, None] * x
        grad_grouped = grouped_sum(grad, se_group)
        num_groups = grad_grouped.shape[0]
        grad2_mean = np.einsum('gi,gj->ij',
                               grad_grouped, grad_grouped) / num_groups
        grad_mean = np.einsum('gi->i', grad_grouped) / num_groups
        grad_cov = grad2_mean - np.outer(grad_mean, grad_mean)

        # Weight by the Hessian.
        xtx_bar = np.einsum('ni,nj,n->ij', x, x, w) / num_groups
        hinv_grad_cov = np.linalg.solve(xtx_bar, grad_cov)
        se2 = np.linalg.solve(xtx_bar, hinv_grad_cov.T) / num_groups
        return se2


def get_regression_w_grads(beta, y, x, w0, se_group=None):
    sens_reg_obj = lambda beta, w: reg_obj(beta, y=y, x=x, w=w)
    obs_w_sens = vittles.HyperparameterSensitivityLinearApproximation(
        objective_fun=sens_reg_obj,
        opt_par_value=beta,
        hyper_par_value=w0,
        validate_optimum=True,
        grad_tol=1e-08)

    get_betahat = obs_w_sens.get_opt_par_function()

    def get_se(w):
        betahat = get_betahat(w)
        se_cov = get_standard_error_matrix(betahat, y, x, w=w, se_group=se_group)
        return np.sqrt(np.diag(se_cov))

    se = get_se(w0)
    betahat_grad = obs_w_sens.get_dopt_dhyper()
    se_grad = autograd.jacobian(get_se)(w0)

    return se, betahat_grad, se_grad


#############################################################
# Sensitivity to the `offset`, i.e to the moment condition E[X eps] = offset.

# This is actually a little silly, since the regression solution is
# linear in the offset.  But this shows how you would to it in general and
# it isn't expensive.
def get_regression_offset_grads(beta, y, x, offset0, se_group=None):
    sens_reg_obj = lambda beta, offset: reg_obj(beta, y=y, x=x, offset=offset)
    offset_sens = vittles.HyperparameterSensitivityLinearApproximation(
        objective_fun=sens_reg_obj,
        opt_par_value=beta,
        hyper_par_value=offset0,
        validate_optimum=True,
        grad_tol=1e-08)

    get_betahat = offset_sens.get_opt_par_function()

    # I believe that using an offset should not affect the values of the
    # standard errors.
    def get_se(offset):
        betahat = get_betahat(offset)
        se_cov = get_standard_error_matrix(
            betahat, y, x, w=np.ones(x.shape[0]), se_group=se_group)
        return np.sqrt(np.diag(se_cov))

    se = get_se(offset0)
    betahat_grad = offset_sens.get_dopt_dhyper()
    se_grad = autograd.jacobian(get_se)(offset0)

    return se, betahat_grad, se_grad



##########################################################
# The below functions are now being done in the R library.

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
