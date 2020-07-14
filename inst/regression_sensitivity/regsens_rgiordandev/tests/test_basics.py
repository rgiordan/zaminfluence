#!/usr/bin/env python3

import autograd
import autograd.numpy as np
from autograd.test_util import check_grads
from numpy.testing import assert_array_almost_equal
import regsens_rgiordandev
import unittest

class TestRegression(unittest.TestCase):
    def generate_data(self, x, beta):
        n_obs = x.shape[0]
        eps = 1.0 * np.random.normal(size=n_obs)
        y = x @ beta + eps
        return y, eps

    def test_regression(self):
        np.random.seed(42)

        dim = 3
        n_obs = 3000
        x_cov = np.eye(dim) + np.full((dim, dim), 0.8)
        x = np.random.random((n_obs, dim)) @ x_cov
        x = x - np.mean(x, axis=0)
        beta_true = np.random.random(dim)
        eps = np.random.normal(n_obs)

        y, eps = self.generate_data(x, beta_true)
        w = np.ones(n_obs)

        assert_array_almost_equal(
            regsens_rgiordandev.reg(y, x),
            np.linalg.solve(x.T @ x, x.T @ y))

        reg_obj_grad = autograd.grad(regsens_rgiordandev.reg_obj)

        # Check that the gradient of the obejctive is zero at the
        # optimal value for a few different weights
        for w in [None, np.ones(n_obs), np.random.random(n_obs) + 0.5]:
            assert_array_almost_equal(
                reg_obj_grad(regsens_rgiordandev.reg(y, x, w=w), y, x, w=w),
                np.zeros(dim))

        # Test the covariance.  Up to numerical error these should all be
        # the same because the data is IID.
        betahat = regsens_rgiordandev.reg(y, x)

        cov_hat0 = regsens_rgiordandev.get_standard_error_matrix(
            betahat, y, x, w, se_group=None)

        cov_hat1 = regsens_rgiordandev.get_standard_error_matrix(
            betahat, y, x, w, se_group=np.arange(n_obs))

        groups = np.floor((n_obs / 2) * np.random.random(n_obs)).astype(int)
        cov_hat2 = regsens_rgiordandev.get_standard_error_matrix(
            betahat, y, x, w, se_group=groups)

        n_sim = 1000
        betas = []
        for n in range(n_sim):
            ysim, _ = self.generate_data(x, beta_true)
            betas.append(regsens_rgiordandev.reg(ysim, x))
        cov_sim = np.cov(np.vstack(betas).T)

        def _assert_cov_approx_eq(cov1, cov2, tol):
            cov_diff = (cov1 - cov2) / np.abs(cov1 + 1e-5)
            print(cov_diff)
            self.assertTrue(np.all(np.abs(cov_diff) < tol))

        _assert_cov_approx_eq(cov_hat0, cov_hat1, 0.4)
        _assert_cov_approx_eq(cov_hat0, cov_hat2, 0.4)
        _assert_cov_approx_eq(cov_hat0, cov_sim, 0.5)

        # For now let's just test that these run.
        betahat = regsens_rgiordandev.reg(y, x)
        w = np.ones(n_obs)

        se, betahat_grad, se_grad = \
            regsens_rgiordandev.get_regression_w_grads(
                betahat, y, x, w)
        se, betahat_grad, se_grad1 = \
            regsens_rgiordandev.get_regression_w_grads(
                betahat, y, x, w, se_group=np.arange(n_obs))
        se, betahat_grad, se_grad2 = \
            regsens_rgiordandev.get_regression_w_grads(
                betahat, y, x, w, se_group=groups)


    def test_offset(self):
        np.random.seed(42)

        dim = 3
        n_obs = 3000
        x_cov = np.eye(dim) + np.full((dim, dim), 0.8)
        x = np.random.random((n_obs, dim)) @ x_cov
        x = x - np.mean(x, axis=0)
        beta_true = np.random.random(dim)

        y, eps = self.generate_data(x, beta_true)
        betahat = regsens_rgiordandev.reg(y, x)

        def offset_reg_obj(beta, offset):
            return regsens_rgiordandev.reg_obj(beta, y, x, offset=offset)
        reg_obj_grad = autograd.grad(offset_reg_obj, argnum=0)
        reg_obj_offset_jac = autograd.jacobian(reg_obj_grad, argnum=1)

        num_tests = 20
        for _ in range(num_tests):
            offset = np.random.random(dim)
            betahat_offset = regsens_rgiordandev.reg(y, x, offset=offset)

            # Check that the gradient is zero at the offset.
            assert_array_almost_equal(
                reg_obj_grad(betahat_offset, offset),
                np.zeros(dim))

            # Check that the offset is actually determining the optimal
            # covariance between the residuals and regressors
            resid = y - x @ betahat_offset
            assert_array_almost_equal(
                offset,
                np.mean(resid[:, None] * x, axis=0))


if __name__ == '__main__':
    unittest.main()
