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
