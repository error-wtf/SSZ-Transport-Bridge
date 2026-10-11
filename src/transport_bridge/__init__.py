"""Transport Bridge — generic typed contracts.

CORE RULE (from docs/TRANSPORT_CONTRACT.md lineage):
the bridge shares TYPES and CHECK STRUCTURE between backends, never
formulas.  Each backend owns its physics completely; the bridge owns
only the vocabulary in which forward/convergence/parity/observable/
inverse contracts are stated and checked.
"""
from .contracts import (
    ObservableSet,
    PropagatedState,
    SignalLaw,
    TransportArchitecture,
    TransportOperator,
    check_convergence,
    check_forward,
    check_inverse_round_trip,
    check_parity,
)
from .phase_contract import (
    TIME_UNITS,
    PhaseContractSpec,
    PhaseDecomposition,
    check_phase_decomposition,
    check_two_channel_readout,
    dimension_check,
)
from .rules import ObservableRule

__all__ = [
    "TIME_UNITS",
    "ObservableRule",
    "ObservableSet",
    "PhaseContractSpec",
    "PhaseDecomposition",
    "PropagatedState",
    "SignalLaw",
    "TransportArchitecture",
    "TransportOperator",
    "check_convergence",
    "check_forward",
    "check_inverse_round_trip",
    "check_parity",
    "check_phase_decomposition",
    "check_two_channel_readout",
    "dimension_check",
]
