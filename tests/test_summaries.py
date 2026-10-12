import numpy as np
from hsp90memory import amm, simulators as sim, summaries as sm


def test_ck_divergence_vanishes_for_markov_two_state():
    m = amm.modal(1); theta = np.log([0.2, 0.65])
    X = sim.simulate_amm(m.Q(theta), m.classes, 216000, n_rep=4, seed=11)
    D = sm.ck_profile([X[i] for i in range(4)])
    assert np.nanmax(np.abs(D)) < 3e-4           # estimator noise level at this data volume


def test_ck_divergence_positive_for_modal_gating():
    m = amm.modal(2)
    X = sim.simulate_amm(m.Q(m.init), m.classes, 216000, n_rep=4, seed=12)
    D = sm.ck_profile([X[i] for i in range(4)])
    assert np.nanmax(D[1]) > 1e-3


def test_summary_vector_is_finite_and_named():
    m = amm.modal(2)
    X = sim.simulate_amm(m.Q(m.init), m.classes, 60000, n_rep=4, seed=13)
    v, names = sm.summary_vector([X[i] for i in range(4)])
    assert len(v) == len(names)
    assert np.isfinite(v[[i for i, n in enumerate(names) if not n.startswith(("ck", "fano"))]]).all()


def _ck_reference(s, tau_frames=10, lags=sm.CK_LAGS):
    """Original masked-mean implementation of the CK divergence."""
    Gtau = sm.propagator(s, tau_frames); out = np.full((2, len(lags)), np.nan)
    for i, L in enumerate(lags):
        if L <= tau_frames: continue
        G = sm.propagator(s, L); Gck = Gtau @ sm.propagator(s, L - tau_frames)
        for q0 in (0, 1):
            out[q0, i] = np.sum(G[q0] * np.log(np.clip(G[q0], 1e-12, None) / np.clip(Gck[q0], 1e-12, None)))
    return out


def test_ck_divergence_counts_identical_to_reference():
    """The prefix-sum/dot-product propagators give bit-identical CK divergences (including a trace that never
    visits one class)."""
    m = amm.modal(2)
    X = sim.simulate_amm(m.Q(m.init), m.classes, 30000, n_rep=2, seed=21)
    for s in [X[0], X[1], np.zeros(20000, np.int8), np.r_[np.zeros(15000, np.int8), np.ones(100, np.int8)]]:
        a, b = sm.ck_divergence(s), _ck_reference(s)
        assert np.array_equal(a, b, equal_nan=True)
