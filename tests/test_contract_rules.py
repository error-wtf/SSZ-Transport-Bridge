"""Transformations are explicit declarations, never inferred from direction tags."""
import math

import pytest

from transport_bridge import ObservableRule, ObservableSet, PropagatedState, check_parity
from transport_bridge.rules import evaluate_rules


@pytest.mark.parametrize(('transform', 'a', 'b'), [
    ('odd', 2.0, -2.0), ('even', 2.0, 2.0), ('reciprocal', 2.0, 0.5),
])
def test_explicit_transform(transform, a, b):
    rule = ObservableRule('x', transform, 1e-12)
    assert evaluate_rules((rule,), {'x': a}, {'x': b}) == {'x': 'PASS'}
    assert evaluate_rules((rule,), {'x': a}, {'x': 3.0}) == {'x': 'FAIL'}


def test_swap_checks_both_channels():
    rule = ObservableRule('a', 'swap', 0.0, partner='b')
    assert evaluate_rules((rule,), {'a': 1, 'b': 2}, {'a': 2, 'b': 1}) == {'a': 'PASS'}
    assert evaluate_rules((rule,), {'a': 1, 'b': 2}, {'a': 9, 'b': 1}) == {'a': 'FAIL'}


def test_custom_and_none():
    rules = (ObservableRule('x', 'custom', 0.0, custom=lambda a, b: a + b == 7),
             ObservableRule('y', 'none', 0.0))
    assert evaluate_rules(rules, {'x': 3}, {'x': 4}) == {'x': 'PASS', 'y': 'OPEN'}


@pytest.mark.parametrize('bad', [math.inf, -math.inf, math.nan])
def test_nonfinite_never_parity_pass(bad):
    obs = ObservableSet(rules=(ObservableRule('x', 'even', 1e-12),))
    out = check_parity(PropagatedState(1, {'x': bad}), PropagatedState(-1, {'x': bad}),
                       obs, 1e-12)
    assert out and not all(out.values())


def test_empty_or_missing_and_zero_reciprocal_fail_closed():
    with pytest.raises(ValueError):
        evaluate_rules((), {}, {})
    assert evaluate_rules((ObservableRule('x', 'odd', 0),), {}, {}) == {'x': 'FAIL'}
    assert evaluate_rules((ObservableRule('x', 'reciprocal', 0),),
                          {'x': 0}, {'x': 0}) == {'x': 'FAIL'}


def test_invalid_rule_rejected():
    for transform in ('automatic', 'ODD'):
        with pytest.raises(ValueError):
            ObservableRule('x', transform, 0)
    with pytest.raises(ValueError):
        ObservableRule('x', 'even', -1)
    with pytest.raises(ValueError):
        ObservableRule('x', 'custom', 0)


def test_direction_tag_does_not_choose_transform():
    obs = ObservableSet(rules=(ObservableRule('x', 'even', 0),))
    assert all(check_parity(PropagatedState(1, {'x': 2}),
                            PropagatedState(-1, {'x': 2}), obs, 0).values())
