"""Backend-independent, explicitly declared observable transformations."""
from __future__ import annotations

import math
from collections.abc import Callable, Mapping
from dataclasses import dataclass


@dataclass(frozen=True)
class ObservableRule:
    name: str
    transform: str
    atol: float
    rtol: float = 0.0
    partner: str | None = None
    custom: Callable[[float, float], bool] | None = None

    def __post_init__(self):
        if self.transform not in {'odd', 'even', 'swap', 'reciprocal', 'custom', 'none'}:
            raise ValueError('unknown observable transformation')
        if not self.name or not all(math.isfinite(t) and t >= 0 for t in (self.atol, self.rtol)):
            raise ValueError('invalid name/tolerance')
        if self.transform == 'swap' and (not self.partner or self.partner == self.name):
            raise ValueError('swap requires a distinct partner')
        if self.transform == 'custom' and not callable(self.custom):
            raise ValueError('custom requires a predicate')


def evaluate_rules(rules: tuple[ObservableRule, ...], plus: Mapping, minus: Mapping) -> dict:
    """PASS/FAIL/OPEN, with finite-only comparisons and no vacuous PASS."""
    if not rules or len({r.name for r in rules}) != len(rules):
        raise ValueError('nonempty unique rules required')
    result = {}
    for rule in rules:
        if rule.transform == 'none':
            result[rule.name] = 'OPEN'
            continue
        try:
            a, b = float(plus[rule.name]), float(minus[rule.name])
            if not (math.isfinite(a) and math.isfinite(b)):
                raise ValueError('nonfinite observable')
            def close(x, y, _rtol=rule.rtol, _atol=rule.atol):
                return math.isfinite(x) and math.isfinite(y) and math.isclose(
                    x, y, rel_tol=_rtol, abs_tol=_atol)

            if rule.transform == 'odd':
                ok = close(a, -b)
            elif rule.transform == 'even':
                ok = close(a, b)
            elif rule.transform == 'reciprocal':
                ok = a != 0 and b != 0 and close(a * b, 1.0)
            elif rule.transform == 'swap':
                ok = close(a, float(minus[rule.partner])) and close(
                    float(plus[rule.partner]), b)
            else:
                ok = rule.custom(a, b) is True
            result[rule.name] = 'PASS' if ok else 'FAIL'
        except (KeyError, ValueError, TypeError, OverflowError, ZeroDivisionError):
            result[rule.name] = 'FAIL'
    return result
