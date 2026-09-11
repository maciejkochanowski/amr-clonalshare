"""The latent-scale restatement of a binary share."""
import math

import numpy as np
import pytest
from scipy.stats import multivariate_normal

from amr_clonalshare.latent import (bivariate_normal_cdf_equal, latent_interval,
                                    latent_share, observed_share)


@pytest.mark.parametrize("rho", [0.05, 0.3, 0.5, 0.9, 0.999])
def test_at_half_prevalence_the_map_is_the_arcsine_law(rho):
    """At a threshold of zero the bivariate orthant probability is
    1/4 + arcsin(rho)/2pi, so the observed share is 2 arcsin(rho)/pi."""
    assert observed_share(rho, 0.5) == pytest.approx(2.0 * math.asin(rho) / math.pi,
                                                     abs=1e-12)


@pytest.mark.parametrize("h,rho", [(-1.2, 0.3), (0.7, 0.8), (1.5, 0.05),
                                   (-0.3, 0.99), (0.0, -0.4)])
def test_the_orthant_probability_matches_scipy(h, rho):
    ours = bivariate_normal_cdf_equal(h, rho)
    ref = multivariate_normal.cdf([h, h], mean=[0.0, 0.0],
                                  cov=[[1.0, rho], [rho, 1.0]])
    assert ours == pytest.approx(ref, abs=1e-10)


@pytest.mark.parametrize("prevalence", [0.02, 0.084, 0.25, 0.5, 0.9])
@pytest.mark.parametrize("observed", [0.01, 0.1, 0.5, 0.935])
def test_the_inverse_map_round_trips(prevalence, observed):
    latent = latent_share(observed, prevalence)
    assert 0.0 < latent < 1.0
    assert observed_share(latent, prevalence) == pytest.approx(observed, abs=1e-8)


def test_the_map_is_monotone_and_the_latent_share_is_the_larger():
    grid = np.linspace(0.01, 0.99, 50)
    images = [observed_share(r, 0.2) for r in grid]
    assert all(b > a for a, b in zip(images, images[1:]))
    assert all(img < r for img, r in zip(images, grid))


def test_the_ends_and_the_undefined_cases():
    assert latent_share(0.0, 0.3) == 0.0
    assert latent_share(-0.05, 0.3) == 0.0
    assert latent_share(1.0, 0.3) == 1.0
    assert math.isnan(latent_share(0.4, 0.0))
    assert math.isnan(latent_share(0.4, 1.0))
    assert math.isnan(latent_share(float("nan"), 0.3))
    lo, hi = latent_interval(-0.1, 0.4, 0.3)
    assert lo == 0.0 and 0.0 < hi < 1.0


def test_the_map_reproduces_a_simulated_threshold_model():
    """Lineages with a normal liability effect, calls by threshold: the
    observed-scale intraclass correlation of the calls is the image of the
    latent share under the map, to Monte Carlo error."""
    rng = np.random.default_rng(20260908)
    n_groups, size, latent, prevalence = 4000, 40, 0.4, 0.15
    a = rng.normal(0.0, math.sqrt(latent), n_groups)
    liability = a[:, None] + rng.normal(0.0, math.sqrt(1.0 - latent),
                                        (n_groups, size))
    from scipy.special import ndtri
    y = (liability > ndtri(1.0 - prevalence)).astype(float)
    means = y.mean(axis=1)
    grand = y.mean()
    msb = size * ((means - grand) ** 2).sum() / (n_groups - 1)
    msw = y.var(axis=1, ddof=1).mean()
    icc = (msb - msw) / (msb + (size - 1) * msw)
    assert icc == pytest.approx(observed_share(latent, prevalence), abs=0.01)
