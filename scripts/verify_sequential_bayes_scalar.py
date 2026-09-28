#!/usr/bin/env python3
"""Verify the synthetic examples in research/sequential_bayes_scalar.tex.

Uses exact Gaussian conditioning and likelihoods; does not run estimation.
"""

import json
from pathlib import Path

import numpy as np
from scipy.stats import multivariate_normal, norm


ROOT = Path(__file__).resolve().parents[1]


def append_observation(mean, covariance, y_new, q, observation_variance=1.0):
    predicted_mean = np.r_[mean, mean[-1]]
    predicted_covariance = np.block([
        [covariance, covariance[:, -1:]],
        [covariance[-1:, :], np.array([[covariance[-1, -1] + q]])],
    ])
    cross_covariance = predicted_covariance[:, -1]
    predictive_variance = predicted_covariance[-1, -1] + observation_variance
    gain = cross_covariance / predictive_variance
    updated_mean = predicted_mean + gain * (y_new - predicted_mean[-1])
    updated_covariance = predicted_covariance - np.outer(cross_covariance, cross_covariance) / predictive_variance
    likelihood = norm.pdf(y_new, predicted_mean[-1], np.sqrt(predictive_variance))
    return updated_mean, updated_covariance, likelihood


def main():
    y = np.array([1.0, 2.0, 3.0])
    old_mean, old_variance = y[:2].sum() / 3, 1 / 3
    sequential_variance = 1 / (1 / old_variance + 1)
    sequential_mean = sequential_variance * (old_mean / old_variance + y[-1])
    np.testing.assert_allclose([sequential_mean, sequential_variance], [1.5, 0.25])

    components, old_evidences, new_evidences = [], [], []
    for q in [1.0, 4.0]:
        t = np.arange(len(y))
        prior_covariance = np.ones((3, 3)) + q * np.minimum.outer(t, t)
        covariance = np.linalg.inv(np.linalg.inv(prior_covariance[:2, :2]) + np.eye(2))
        mean = covariance @ y[:2]
        sequential_mean, sequential_covariance, predictive = append_observation(mean, covariance, y[-1], q)
        batch_covariance = np.linalg.inv(np.linalg.inv(prior_covariance) + np.eye(3))
        batch_mean = batch_covariance @ y
        np.testing.assert_allclose(sequential_mean, batch_mean, rtol=1e-13, atol=1e-14)
        np.testing.assert_allclose(sequential_covariance, batch_covariance, rtol=1e-13, atol=1e-14)
        old_evidence = 0.5 * multivariate_normal.pdf(y[:2], cov=prior_covariance[:2, :2] + np.eye(2))
        batch_evidence = 0.5 * multivariate_normal.pdf(y, cov=prior_covariance + np.eye(3))
        np.testing.assert_allclose(old_evidence * predictive, batch_evidence, rtol=1e-13)
        if q == 1:
            np.testing.assert_allclose(sequential_mean, np.array([12, 23, 31]) / 13)
            np.testing.assert_allclose(sequential_covariance, np.array([[5, 2, 1], [2, 6, 3], [1, 3, 8]]) / 13)
        else:
            np.testing.assert_allclose(sequential_mean, np.array([21, 61, 89]) / 32)
            np.testing.assert_allclose(sequential_covariance, np.array([[29, 5, 1], [5, 45, 9], [1, 9, 53]]) / 64)
        old_evidences.append(old_evidence)
        new_evidences.append(batch_evidence)
        components.append(dict(q=q, old_mean=mean.tolist(), old_covariance=covariance.tolist(),
                               predictive_density=float(predictive), new_mean=sequential_mean.tolist(),
                               new_covariance=sequential_covariance.tolist()))
    old_weights = np.array(old_evidences) / sum(old_evidences)
    new_weights = np.array(new_evidences) / sum(new_evidences)
    sequential_weights = old_weights * [c["predictive_density"] for c in components]
    sequential_weights /= sequential_weights.sum()
    np.testing.assert_allclose(sequential_weights, new_weights, rtol=1e-13)
    new_mean = new_weights @ np.array([c["new_mean"] for c in components])
    new_covariance = sum(w * (np.array(c["new_covariance"]) +
                             np.outer(np.array(c["new_mean"]) - new_mean, np.array(c["new_mean"]) - new_mean))
                         for w, c in zip(new_weights, components))
    results = dict(status="passed", tolerance=1e-13, observations=y.tolist(),
                   static_posterior=dict(mean=1.5, variance=0.25), components=components,
                   old_parameter_weights=old_weights.tolist(), new_parameter_weights=new_weights.tolist(),
                   old_parameter_mean=float(old_weights @ np.array([1, 4])),
                   new_parameter_mean=float(new_weights @ np.array([1, 4])),
                   mixture_smoothed_mean=new_mean.tolist(), mixture_smoothed_covariance=new_covariance.tolist())
    destination = ROOT / "research/generated/sequential_bayes_scalar_checks.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
