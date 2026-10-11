import numpy as np
from hsp90memory import amm, simulators as sim
from hsp90memory.io import dwells


def _dwells_from(X):
    out = []
    for s in X:
        dw, dc = dwells(s); out.append((dw, dc))
    return out


def test_fast_and_slow_likelihood_agree():
    m = amm.modal(2)
    X = sim.simulate_amm(m.Q(m.init), m.classes, 20000, n_rep=2, seed=3)
    dl = _dwells_from(X)
    th = m.init + 0.2
    assert abs(amm.loglik(m, th, dl) - amm.loglik_slow(m, th, dl)) < 1e-6


def test_stationary_distribution():
    Q = np.array([[-1.0, 1.0], [2.0, -2.0]])
    pi = amm.stationary(Q)
    assert np.allclose(pi, [2 / 3, 1 / 3])


def test_two_state_rates_recovered():
    """Markov two-state chain: ML rates must match the generating ones within a few percent on 2 x 10^5 frames."""
    m = amm.modal(1); theta = np.log([0.2, 0.65])
    X = sim.simulate_amm(m.Q(theta), m.classes, 200000, n_rep=2, seed=5)
    res = amm.fit(m, _dwells_from(X), n_restart=1, maxiter=100)
    assert np.allclose(np.exp(res.x), np.exp(theta), rtol=0.05)


def test_modal_ou_sigma_zero_matches_amm_populations():
    """With sigma = 0 the nested simulator is the modal AMM (up to the first-order mode-switch step)."""
    m = amm.modal(2); Q = m.Q(m.init)
    kOC, kCO, SO, SC = sim.modal_params_from_Q(Q)
    X = sim.simulate_modal_ou(kOC, kCO, SO, SC, 0.0, 10.0, 0.0, 100000, n_rep=4, seed=1)
    Y = sim.simulate_amm(Q, m.classes, 100000, n_rep=4, seed=2)
    assert abs(X.mean() - Y.mean()) < 0.03
