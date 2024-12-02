import pytest
from pypythia.raxmlng_parser import (
    get_patterns_gaps_invariant,
    get_raxmlng_rfdist_results,
)


def test_get_patterns_gaps_invariant(raxmlng, raxmlng_inference_log):
    patterns, gaps, invariant = get_patterns_gaps_invariant(raxmlng_inference_log)

    assert isinstance(patterns, int)
    assert isinstance(gaps, float)
    assert isinstance(invariant, float)

    assert patterns == 50
    assert gaps == pytest.approx(0.0163, abs=0.01)
    assert invariant == pytest.approx(0.8970, abs=0.01)


def test_get_patterns_gaps_invariant_raises_value_error(raxmlng_rfdistance_log):
    with pytest.raises(ValueError):
        get_patterns_gaps_invariant(raxmlng_rfdistance_log)


def test_get_raxmlng_rfdist_results(raxmlng_rfdistance_log):
    num_topos, rel_rfdist, abs_rfdist = get_raxmlng_rfdist_results(
        raxmlng_rfdistance_log
    )

    assert num_topos == 2
    assert rel_rfdist == pytest.approx(1 / 3, abs=0.01)
    assert abs_rfdist == pytest.approx(2)


def test_get_raxmlng_rfdist_results_raises_value_error(raxmlng_inference_log):
    with pytest.raises(ValueError):
        test_get_raxmlng_rfdist_results(raxmlng_inference_log)
