"""Exact Bernoulli boundaries must not manufacture evidence, and the clipped predictor keeps the process alive."""
from itertools import product

import numpy as np
import pytest

from amr_clonalshare.evalues import _log_bernoulli, e_process, sequential_e_process


@pytest.mark.parametrize("outcome", [0., 1.])
def test_deterministic_null_cannot_earn_evidence_from_singleton_training(outcome):
    # Under p=0 or p=1 this is the only possible sequence, so E must be <=1
    # pointwise, without a Monte Carlo tolerance or clipping allowance.
    ys = [np.full(2, outcome) for _ in range(3)]
    labels = [["a", "b"], ["a", "b"], ["a", "b"]]
    result = sequential_e_process(ys, labels)
    assert all(value <= 1. for value in result.e_value)
    assert all(value <= 0. for value in result.log_e)
    split = e_process(np.full(6, outcome), list("abcdef"), folds=3, repeats=2, seed=0)
    assert split.e_value <= 1.
    assert split.log_e <= 0.


@pytest.mark.parametrize("y,p,expected", [
    ([0., 1.], [0., 1.], 0.),
    ([0., 1.], [1., 0.], -np.inf),
    ([0., 1.], [0.25, 0.75], 2 * np.log(0.75)),
])
def test_bernoulli_likelihood_preserves_certain_and_impossible_events(y, p, expected):
    with np.errstate(all="raise"):
        assert _log_bernoulli(np.array(y), np.array(p)) == expected


def test_a_training_singleton_predicting_the_opposite_call_leaves_a_floor():
    # Each of the two training singletons predicts the opposite held-out
    # call. The predictor is clipped away from 0 and 1, so the evidence is
    # the floor of the clip rather than an exact zero, and no NaN arises.
    with np.errstate(all="raise"):
        result = e_process([0, 1], ["a", "b"], folds=2, repeats=1, seed=0)
    assert result.e_value == pytest.approx(1e-6, rel=1e-9)
    assert result.log_e == pytest.approx(np.log(1e-6), rel=1e-9)
    assert not result.reject_05


def test_zero_and_nonzero_split_likelihoods_have_the_correct_average():
    # Leave-one-out likelihood ratios are the clip floor, 1/2 and 1/2, so
    # their mean is 1/3 up to a third of the floor.
    with np.errstate(all="raise"):
        result = e_process([0, 0, 1], list("abc"), folds=3, repeats=1, seed=0)
    assert result.e_value == pytest.approx(1 / 3 + 1e-6 / 3, abs=1e-12)
    assert result.log_e == pytest.approx(np.log(1 / 3 + 1e-6 / 3), abs=1e-12)


def test_a_first_batch_of_identical_outcomes_does_not_end_the_process():
    # Trained on a batch of two susceptible calls, the predictor gives the
    # two resistant calls of the next batch the clipped floor, not zero, so
    # the running product is penalised and can recover from later batches.
    result = sequential_e_process([[0, 0], [1, 1], [0, 1]], [["a", "b"]] * 3)
    assert result.e_value[0] == 1.
    assert result.e_value[1] == pytest.approx(1e-12, rel=1e-9)
    assert result.log_e[1] == pytest.approx(2 * np.log(1e-6), rel=1e-9)
    assert result.log_e[2] == pytest.approx(result.log_e[1], rel=1e-9)
    assert all(np.isfinite(result.log_e))


def test_large_split_evidence_retains_finite_log_after_natural_scale_overflow():
    labels = np.repeat([0, 1], 1200)
    result = e_process(labels, labels, folds=2, repeats=1, seed=0)
    assert np.isfinite(result.log_e) and result.log_e > 710
    assert result.e_value == np.inf
    assert result.reject_05 and result.reject_01


@pytest.mark.parametrize("bad", [-1., 0.5, 2.])
@pytest.mark.parametrize("mode", ["split", "sequential"])
def test_finite_nonbinary_outcomes_are_rejected(bad, mode):
    with pytest.raises(ValueError, match="binary"):
        if mode == "split":
            e_process([0., bad, 1., 0.], list("abcd"), seed=0)
        else:
            sequential_e_process([[0., bad], [1., 0.]], [["a", "b"]] * 2)


@pytest.mark.parametrize("p", [0., 0.2, 0.5, 0.8, 1.])
def test_exhaustive_small_bernoulli_null_expectations(p):
    # Sum all 2**6 sequences, with their exact Bernoulli null probabilities.
    # This exercises singleton and replicated training designs at each look.
    means = np.zeros(3)
    split_mean = 0.
    for calls in product([0, 1], repeat=6):
        k = sum(calls)
        mass = p**k * (1-p)**(6-k)
        if not mass:
            continue
        sequential = sequential_e_process([calls[:2], calls[2:4], calls[4:]],
                                          [["a", "b"]] * 3)
        means += mass * np.asarray(sequential.e_value)
        split_mean += mass * e_process(calls, list("abcdef"), folds=3,
                                       repeats=1, seed=0).e_value
    assert np.all(means <= 1. + 1e-14)
    assert split_mean <= 1. + 1e-14


@pytest.mark.parametrize("mode", ["split", "sequential"])
def test_binary_validation_retains_missing_outcome_and_untyped_accounting(mode):
    calls = [0., 0., np.nan, np.inf, 1.]
    labels = ["a", "b", "c", "d", None]
    if mode == "split":
        result = e_process(calls, labels, folds=2, repeats=1, seed=0)
    else:
        result = sequential_e_process([calls], [labels])
    assert result.n_dropped_non_finite == 2
    assert result.n_dropped_untyped == 1
