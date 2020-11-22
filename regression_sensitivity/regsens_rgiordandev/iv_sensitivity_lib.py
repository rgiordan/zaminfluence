# Note that these are tested in the R package zaminfluence to make sure they
# match the output of `ivreg`.

import autograd
import autograd.numpy as np
from copy import deepcopy

import paragami
from paragami.autograd_supplement_lib import grouped_sum
import vittles
from vittles import solver_lib
from vittles.sensitivity_lib import EstimatingEquationLinearApproximation

def iv_reg(y, x, z, w=None):
    n_obs = x.shape[0]
    assert len(y) == n_obs
    assert z.shape[0] == n_obs
    assert x.shape[1] == z.shape[1]
    if w is None:
        w = np.ones(n_obs)
    zx_cov = np.einsum('ni,nj,n->ij', z, x, w) / n_obs
    zy_cov = np.einsum('ni,n,n->i', z, y, w) / n_obs
    return np.linalg.solve(zx_cov, zy_cov)


def iv_reg_estimating_equation(beta, y, x, z, w=None):
    n_obs = x.shape[0]
    assert len(y) == n_obs
    assert z.shape[0] == n_obs
    assert x.shape[1] == z.shape[1]
    if w is None:
        w = np.ones(n_obs)
    zx_cov = np.einsum('ni,nj,n->ij', z, x, w) / n_obs
    zy_cov = np.einsum('ni,n,n->i', z, y, w) / n_obs
    return zx_cov @ beta - zy_cov


def get_iv_standard_error_matrix(betahat, y, x, z, w, se_group=None):
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
        ztz_bar = np.einsum('ni,nj,n->ij', z, z, w) / num_obs
        ztx_bar = np.einsum('ni,nj,n->ij', z, x, w) / num_obs
        sigma2hat = np.sum((w * resid) ** 2) / (num_obs - len(betahat))
        ztx_inv = np.linalg.inv(ztx_bar)
        se2 = sigma2hat * np.linalg.solve(ztx_bar, ztz_bar) @ ztx_inv.T / num_obs
        return se2
    else:
        if len(se_group) != len(y):
            raise ValueError('se_group must be the same length as the data.')

        if np.min(se_group) != 0:
            raise ValueError('se_group must be zero-indexed ' +
                             '(its minimum must be zero)')

        # Calculate the sample variance of the gradient where each group
        # is treated as a single observation.
        grad = w[:, None] * resid[:, None] * z
        grad_grouped = grouped_sum(grad, se_group)
        num_groups = grad_grouped.shape[0]
        grad2_mean = np.einsum('gi,gj->ij',
                               grad_grouped, grad_grouped) / num_groups
        grad_mean = np.einsum('gi->i', grad_grouped) / num_groups
        grad_cov = grad2_mean - np.outer(grad_mean, grad_mean)

        # Weight by the Hessian.
        ztx_bar = np.einsum('ni,nj,n->ij', z, x, w) / num_groups
        hinv_grad_cov = np.linalg.solve(ztx_bar, grad_cov)
        se2 = np.linalg.solve(ztx_bar, hinv_grad_cov.T) / num_groups
        return se2


def get_iv_regression_w_grads(beta, y, x, z, w0, se_group=None):
    sens_iv_reg_obj = \
        lambda beta, w: iv_reg_estimating_equation(beta, y=y, x=x, z=z, w=w)

    hess0 = autograd.jacobian(sens_iv_reg_obj, argnum=0)(beta, w0)

    # The `Hessian` is not symmetric; you cannot use a Cholesky solver!
    def hess_solver(v):
        return np.linalg.solve(hess0, v)

    obs_w_sens = EstimatingEquationLinearApproximation(
        estimating_equation=sens_iv_reg_obj,
        input_val0=beta,
        hyper_val0=w0,
        hess_solver=hess_solver,
        validate_solution=True,
        solution_tol=1e-08)
    betahat_grad = obs_w_sens.get_dinput_dhyper()

    # Standard errors
    get_betahat = obs_w_sens.get_input_par_function()
    def get_se(w):
        betahat = get_betahat(w)
        se_cov = get_iv_standard_error_matrix(
            betahat=betahat, y=y, x=x, z=z, w=w, se_group=se_group)
        return np.sqrt(np.diag(se_cov))

    se = get_se(w0)
    se_grad = autograd.jacobian(get_se)(w0)

    return se, betahat_grad, se_grad
