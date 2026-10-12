import json, warnings
from pathlib import Path
import numpy as np
from hsp90memory import simulators as sim, summaries as sm

ROOT = Path(__file__).resolve().parents[1]
NT, NF, R = 4, 36000, 20            # 20 replicate measurements of 4 x 1 h


def _summaries(X):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.array([sm.summary_vector([X[k * NT + j] for j in range(NT)])[0] for k in range(R)])


def _z(A, B):
    ok = np.isfinite(A).all(0) & np.isfinite(B).all(0)
    A, B = A[:, ok], B[:, ok]
    se = np.sqrt((A.var(0, ddof=1) + B.var(0, ddof=1)) / R)
    keep = se > 0
    return (A.mean(0) - B.mean(0))[keep] / se[keep]


def test_modal_ou_numba_matches_reference():
    """Numba and numpy simulators agree on every summary statistic within replicate tolerance (fixed seeds, modal3
    mode structure, state-specific disorder with correlation), and the comparison has power: the numba simulator
    without disorder is rejected."""
    Q = np.array(json.load(open(ROOT / "results/modal3.json"))["Q"])
    p = sim.modal_params_from_Q(Q)
    args = ((0.6, 0.4), 20.0, 0.5, NF)
    ref = _summaries(sim.simulate_modal_ou(*p, *args, n_rep=R * NT, seed=7))
    nb = _summaries(sim.simulate_modal_ou_nb(*p, *args, n_rep=R * NT, seed=8))
    z = _z(ref, nb)
    # ~110 correlated statistics: expect ~5 % beyond 2 sd and a maximum around 3 sd under equality
    assert z.size > 100 and np.abs(z).max() < 4.5 and (np.abs(z) > 2).mean() < 0.1, np.abs(z).max()
    nb0 = _summaries(sim.simulate_modal_ou_nb(*p, 0.0, 20.0, 0.5, NF, n_rep=R * NT, seed=9))
    assert np.abs(_z(ref, nb0)).max() > 6.0


def test_modal_ou_numba_reproducible():
    Q = np.array(json.load(open(ROOT / "results/modal3.json"))["Q"])
    p = sim.modal_params_from_Q(Q)
    a = sim.simulate_modal_ou_nb(*p, 0.5, 10.0, 0.0, 5000, n_rep=3, seed=123)
    b = sim.simulate_modal_ou_nb(*p, 0.5, 10.0, 0.0, 5000, n_rep=3, seed=123)
    assert np.array_equal(a, b) and a.dtype == np.int8 and a.shape == (3, 5000)


def test_modal_ou_fast_identical_to_reference():
    """simulate_modal_ou_fast consumes the reference's random draws in the same order: identical trajectories."""
    for key in ("modal3", "modal4"):
        Q = np.array(json.load(open(ROOT / f"results/{key}.json"))["Q"])
        p = sim.modal_params_from_Q(Q)
        for sigma, tau, rho in [(0.0, 10.0, 0.0), ((0.6, 0.1), 300.0, -0.7), ((0.0, 0.8), 1.0, 1.0)]:
            a = sim.simulate_modal_ou(*p, sigma, tau, rho, 30000, n_rep=3, seed=5)
            b = sim.simulate_modal_ou_fast(*p, sigma, tau, rho, 30000, n_rep=3, seed=5)
            assert np.array_equal(a, b), (key, sigma, tau, rho)
