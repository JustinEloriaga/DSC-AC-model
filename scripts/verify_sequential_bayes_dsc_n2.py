#!/usr/bin/env python3
"""Check the algebra in sequential_bayes_dsc_n2.tex without fitting a model.

Synthetic evaluation points exercise covariance/density identities only. No
numerical hyperparameters or posterior estimates are supplied to the document.
"""

import json
from pathlib import Path

import numpy as np
from scipy.linalg import block_diag, logm
from scipy.stats import invgamma, invwishart, multivariate_normal, norm


ROOT = Path(__file__).resolve().parents[1]


def rw_covariance(length, initial_scale, variance):
    times = np.arange(1, length + 1)
    return variance * (initial_scale + np.minimum.outer(times, times))


def main():
    rng = np.random.default_rng(928)
    errors = {}

    def check(name, left, right):
        np.testing.assert_allclose(left, right, atol=2e-10, rtol=2e-10)
        errors[name] = max(errors.get(name, 0.0), float(np.max(np.abs(np.asarray(left) - right))))

    for _ in range(64):
        m_b = rng.normal(size=2)
        initial_factor = rng.normal(size=(2, 2))
        c_b = initial_factor @ initial_factor.T + np.eye(2)
        scale_factor = rng.normal(size=(2, 2))
        psi = scale_factor @ scale_factor.T + np.eye(2)
        nu, a = rng.uniform(4, 9), rng.uniform(2, 7)
        b_h, b_r, c_h, c_r = np.exp(rng.uniform(-1, 1, size=4))
        q_h = invgamma.rvs(a, scale=b_h, size=2, random_state=rng)
        q_r = invgamma.rvs(a, scale=b_r, random_state=rng)
        v = invwishart.rvs(nu, psi, random_state=rng)
        mu_h, mu_r = rng.normal(size=2), rng.normal()
        b_path = np.empty((4, 2))
        b_path[0] = rng.multivariate_normal(m_b, c_b)
        for t in range(1, 4):
            b_path[t] = rng.multivariate_normal(b_path[t - 1], v)
        h_path = np.column_stack([
            rng.multivariate_normal(np.full(4, mu_h[j]), rw_covariance(4, c_h, q_h[j]))
            for j in range(2)
        ])
        r_path = rng.multivariate_normal(np.full(4, mu_r), rw_covariance(4, c_r, q_r))
        # Keep likelihood checks away from floating-point correlation saturation.
        r_path = np.clip(r_path, -2.5, 2.5)
        y = b_path + rng.normal(size=(4, 2)) * np.exp(h_path / 2)

        likelihood = []
        for t in range(4):
            rho = np.tanh(r_path[t])
            correlation = np.array([[1, rho], [rho, 1]])
            check("matrix_log_coordinate", logm(correlation)[1, 0], r_path[t])
            sd = np.exp(h_path[t] / 2)
            sigma = np.outer(sd, sd) * correlation
            z = (y[t] - b_path[t]) / sd
            explicit = (-np.log(2 * np.pi) - h_path[t].sum() / 2
                        - np.log1p(-rho**2) / 2
                        - (z[0]**2 - 2 * rho * z[0] * z[1] + z[1]**2) / (2 * (1 - rho**2)))
            dense = multivariate_normal.logpdf(y[t], b_path[t], sigma)
            check("bivariate_likelihood", explicit, dense)
            likelihood.append(dense)

        mean_covariance = np.block([
            [c_b + min(t, s) * v for s in range(4)] for t in range(4)
        ])
        b_prior = multivariate_normal.logpdf(b_path.ravel(), np.tile(m_b, 4), mean_covariance)
        b_initial = multivariate_normal.logpdf(b_path[0], m_b, c_b)
        b_transitions = [multivariate_normal.logpdf(b_path[t], b_path[t - 1], v) for t in range(1, 4)]
        check("mean_path_covariance_factorization", b_prior, b_initial + sum(b_transitions))

        coordinates = [
            (h_path[:, j], q_h[j], mu_h[j], c_h, b_h) for j in range(2)
        ] + [(r_path, q_r, mu_r, c_r, b_r)]
        initial_logpdf = b_initial
        transition_logpdfs = np.array(b_transitions)
        dense_prior = b_prior
        old_dense_prior = multivariate_normal.logpdf(
            b_path[:3].ravel(), np.tile(m_b, 3), mean_covariance[:6, :6]
        )
        for path, q, mu, c, scale in coordinates:
            cov = rw_covariance(4, c, q)
            path_density = multivariate_normal.logpdf(path, np.full(4, mu), cov)
            initial_density = norm.logpdf(path[0], mu, np.sqrt((c + 1) * q))
            transition_densities = norm.logpdf(np.diff(path), 0, np.sqrt(q))
            check("scalar_path_covariance_factorization", path_density, initial_density + transition_densities.sum())
            check("scalar_initial_schur_complement", cov[3, 3] - cov[3, :3] @ np.linalg.solve(cov[:3, :3], cov[:3, 3]), q)
            dense_prior += path_density
            old_dense_prior += multivariate_normal.logpdf(path[:3], np.full(3, mu), cov[:3, :3])
            initial_logpdf += initial_density
            transition_logpdfs += transition_densities

            # A conditional density and its joint kernel must differ by a constant in q.
            sum_squares = (path[0] - mu)**2 / (c + 1) + np.diff(path) @ np.diff(path)
            offsets = []
            for candidate in np.exp(np.linspace(-1.4, 1.4, 5)):
                joint = invgamma.logpdf(candidate, a, scale=scale) + multivariate_normal.logpdf(
                    path, np.full(4, mu), rw_covariance(4, c, candidate)
                )
                conditional = invgamma.logpdf(candidate, a + 2, scale=scale + sum_squares / 2)
                offsets.append(joint - conditional)
            check("inverse_gamma_full_conditional", np.array(offsets), offsets[0])

        prior_theta = (invwishart.logpdf(v, nu, psi)
                       + invgamma.logpdf(q_h, a, scale=b_h).sum()
                       + invgamma.logpdf(q_r, a, scale=b_r))
        batch_kernel = prior_theta + dense_prior + sum(likelihood)
        old_kernel = prior_theta + old_dense_prior + sum(likelihood[:3])
        sequential_kernel = old_kernel + transition_logpdfs[-1] + likelihood[-1]
        check("batch_vs_sequential_joint_kernel", batch_kernel, sequential_kernel)

        # Each five-state transition retains V's off-diagonal dependence.
        transition_covariance = block_diag(v, np.diag(q_h), np.array([[q_r]]))
        delta = np.r_[b_path[3] - b_path[2], h_path[3] - h_path[2], r_path[3] - r_path[2]]
        check("joint_transition_block_covariance", multivariate_normal.logpdf(delta, cov=transition_covariance), transition_logpdfs[-1])

        delta_b = np.diff(b_path, axis=0)
        posterior_psi = psi + delta_b.T @ delta_b
        offsets = []
        for _ in range(5):
            candidate = invwishart.rvs(nu, psi, random_state=rng)
            joint = invwishart.logpdf(candidate, nu, psi)
            joint += sum(multivariate_normal.logpdf(d, cov=candidate) for d in delta_b)
            conditional = invwishart.logpdf(candidate, nu + 3, posterior_psi)
            offsets.append(joint - conditional)
        check("inverse_wishart_full_conditional", np.array(offsets), offsets[0])

    results = {
        "status": "passed",
        "synthetic_evaluation_cases": 64,
        "max_absolute_errors": errors,
        "scope": "Algebra checks, not posterior estimation or numerical evidence integration."
    }
    destination = ROOT / "research/generated/sequential_bayes_dsc_n2_checks.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
