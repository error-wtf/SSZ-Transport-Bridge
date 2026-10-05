"""check_convergence contract: order estimation and monotonicity must
work on analytically known grid problems."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from transport_bridge import check_convergence


def _grid_errors(order: float, base_err: float):
    """Error sequence for successively doubled grids."""
    return [base_err / (2 ** (order * i)) for i in range(5)]


def test_second_order_detected():
    errs = _grid_errors(2.0, 1e-2)
    res = check_convergence(lambda g: errs[g], list(range(5)),
                            min_observed_order=1.5)
    assert res["order_ok"]
    assert res["observed_order"] == pytest.approx(2.0, abs=0.1)
    assert res["monotone_decrease"]


def test_first_order_rejected_at_high_bar():
    errs = _grid_errors(1.0, 1e-2)
    res = check_convergence(lambda g: errs[g], list(range(5)),
                            min_observed_order=1.5)
    assert not res["order_ok"]
    assert res["observed_order"] == pytest.approx(1.0, abs=0.1)


def test_non_monotone_flagged():
    errs = [1e-2, 1e-3, 5e-3, 1e-4, 1e-5]  # bump in the middle
    res = check_convergence(lambda g: errs[g], list(range(5)),
                            min_observed_order=1.5)
    assert not res["monotone_decrease"]


def test_zero_error_grid_handled():
    # exact solution on the finest grid: last diff = 0
    errs = [1e-2, 2.5e-3, 6.25e-4, 1.5e-4, 0.0]
    res = check_convergence(lambda g: errs[g], list(range(5)),
                            min_observed_order=1.5)
    # order estimate involves log(0) -> nan is acceptable, but no crash
    assert "observed_order" in res
