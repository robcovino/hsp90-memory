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


def test_loglik_Q_matches_reference_and_terms_sum():
    m = amm.modal(3)
    X = sim.simulate_amm(m.Q(m.init), m.classes, 30000, n_rep=2, seed=7)
    dl = _dwells_from(X); th = m.init + 0.1
    a = amm.loglik(m, th, dl); b = amm.loglik_Q(m.Q(th), m.classes, dl)
    assert abs(a - b) < 1e-6
    assert abs(sum(t.sum() for t in amm.loglik_terms(m.Q(th), m.classes, dl)) - b) < 1e-6


def test_modal_ou_reduces_to_modal_and_eq_is_balanced():
    m = amm.modal(2); X = sim.simulate_amm(m.Q(m.init), m.classes, 20000, n_rep=2, seed=8); dl = _dwells_from(X)
    for couple in ("OC", "eq"):
        mo = amm.modal_ou(2, n_grid=7, couple=couple)
        th = np.r_[m.init, np.log(1e-9), np.log(50.0)]
        assert abs(amm.loglik_Q(mo.Q(th), mo.classes, dl) - amm.loglik(m, m.init, dl)) < 1e-5
    mo = amm.modal_ou(3, n_grid=9, couple="eq"); Q = mo.Q(np.r_[amm.modal(3).init, np.log(0.5), np.log(100.0)])
    F = amm.stationary(Q)[:, None] * Q; np.fill_diagonal(F, 0)
    assert np.abs(F - F.T).max() < 1e-12 * F.max()
    z, R = amm.ou_grid(15); pz = amm.stationary(R)
    assert abs(np.sum(pz * z ** 2) - 1) < 1e-12 and abs(np.sort(np.linalg.eigvals(R).real)[-2] + 1) < 1e-10


def test_grid_ou_amm_matches_ou_simulator():
    """The grid-OU AMM (simulated exactly from its generator) and simulate_modal_ou_fast (sigma_CO = 0) are the same
    model: all summary statistics agree within replicate tolerance."""
    import sys; from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent)); import test_simulators as T
    m4 = amm.modal(4); th = m4.init
    mo = amm.modal_ou(4, n_grid=21)
    A = T._summaries(sim.simulate_amm(mo.Q(np.r_[th, np.log(0.6), np.log(60.0)]), mo.classes, T.NF, n_rep=T.R * T.NT, seed=3))
    B = T._summaries(sim.simulate_modal_ou_fast(*sim.modal_params_from_Q(m4.Q(th)), (0.6, 0.0), 60.0, 0.0, T.NF,
                                                 n_rep=T.R * T.NT, seed=4))
    z = T._z(A, B)
    assert np.abs(z).max() < 4.5 and (np.abs(z) > 2).mean() < 0.1, np.abs(z).max()


def test_ring_detailed_balance_and_driven_variant():
    rng = np.random.default_rng(1)
    for db in (True, False):
        r = amm.ring(4, db=db); Q = r.Q(r.init + rng.normal(0, 0.5, r.n_par))
        F = amm.stationary(Q)[:, None] * Q; np.fill_diagonal(F, 0)
        asym = np.abs(F - F.T).max() / F.max()
        assert (asym < 1e-12) if db else (asym > 1e-6)
    mo = amm.modal_ou(4, n_grid=7, couple="eq", topology="ring")
    Q = mo.Q(np.r_[amm.ring(4).init, np.log(0.4), np.log(100.0)])
    F = amm.stationary(Q)[:, None] * Q; np.fill_diagonal(F, 0)
    assert np.abs(F - F.T).max() < 1e-12 * F.max()


def test_graph_reproduces_ring_and_is_balanced():
    rng = np.random.default_rng(2)
    r = amm.ring(4); th = r.init + rng.normal(0, 0.5, r.n_par)
    g = amm.graph(4, amm.RING4)
    assert g.n_par == r.n_par and np.abs(g.Q(th) - r.Q(th)).max() < 1e-14
    for name in amm.NAMED_GRAPHS:
        m = amm.named_graph(name); Q = m.Q(m.init + rng.normal(0, 0.5, m.n_par))
        F = amm.stationary(Q)[:, None] * Q; np.fill_diagonal(F, 0)
        assert np.abs(F - F.T).max() < 1e-12 * F.max()


def test_split_mode_balance():
    """Equilibrium sub-state chain keeps detailed balance; the driven chain breaks it; OU can be layered on top."""
    g = amm.named_graph("ringF1"); rng = np.random.default_rng(3); th = g.init + rng.normal(0, 0.3, g.n_par)
    for driven in (False, True):
        sm = amm.split_mode(g, 3, 2, driven=driven); Q = sm.Q(np.r_[th, np.log(0.01)])
        F = amm.stationary(Q)[:, None] * Q; np.fill_diagonal(F, 0)
        asym = np.abs(F - F.T).max() / F.max()
        assert (asym > 1e-6) if driven else (asym < 1e-12)
        assert sm.K == g.K + 2 and sm.n_par == g.n_par + 1
    mo = amm.modal_ou(0, n_grid=5, couple="eq", base=amm.split_mode(g, 3, 2))
    Q = mo.Q(np.r_[th, np.log(0.01), np.log(0.3), np.log(200.0)]); F = amm.stationary(Q)[:, None] * Q; np.fill_diagonal(F, 0)
    assert np.abs(F - F.T).max() < 1e-12 * F.max()
