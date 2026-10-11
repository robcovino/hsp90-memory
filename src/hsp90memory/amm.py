"""Aggregated (hidden) Markov models for the two-level plasmon-ruler traces.

Observation model: each frame is an instantaneous snapshot (no straddle frames, see within_state.py), the
SNR is ~11, so thresholding gives an error-free dwell sequence in units of frames.  The likelihood of a
continuous-time Markov chain with generator Q observed at interval dt through a 2-class aggregation is

    L = pi_c0 . prod_i [ (P_cc)^(n_i-1) P_cc' ] . (P_cc)^(n_N-1) . 1 ,   P = expm(Q dt),

where P_cc is the sub-block of P between states of class c.  Using P = expm(Q dt) (rather than Q dt) makes
sub-frame excursions ("missed events") exact at the sampling resolution.

Model families
  modal(n_modes):  n_modes gating modes in a linear chain; each mode has its own O<->C rates; mode switching
                   from both O and C; detailed balance imposed on every cycle (equilibrium, apo).
  modal_nodb(2):   as modal(2) without the detailed-balance constraint (likelihood-ratio test of DB).
  general(kO,kC):  fully connected kO open + kC closed states, all rates free (likelihood ceiling).
"""
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

from .io import DT, load_dwells  # noqa: F401


# ----------------------------------------------------------------------------- models
class Model:
    """Holds the state->class map and a function theta -> Q."""
    def __init__(self, name, classes, build, n_par, init):
        self.name, self.classes, self.build, self.n_par, self.init = name, np.asarray(classes), build, n_par, init
        self.K = len(classes)

    def Q(self, theta):
        return self.build(np.exp(theta))


def _finalize(Q):
    Q = Q.copy(); np.fill_diagonal(Q, 0); np.fill_diagonal(Q, -Q.sum(1)); return Q


def modal(n_modes, db=True):
    """States ordered (O_0, C_0, O_1, C_1, ...), mode 0 slowest.  Rates per mode: kOC_m, kCO_m.
    Mode switching m<->m+1 from O: a_m (up), b_m (down); from C: c_m (up), d_m (down).
    With detailed balance d_m is fixed by the cycle condition kOC_m * c_m * kCO_{m+1} * b_m = kCO_m * a_m * kOC_{m+1} * d_m."""
    K = 2 * n_modes
    classes = [0, 1] * n_modes
    n_sw = 3 if db else 4
    n_par = 2 * n_modes + n_sw * (n_modes - 1)

    def build(r):
        Q = np.zeros((K, K))
        kOC = r[0:2 * n_modes:2]; kCO = r[1:2 * n_modes:2]
        sw = r[2 * n_modes:].reshape(n_modes - 1, n_sw)
        for m in range(n_modes):
            Q[2 * m, 2 * m + 1] = kOC[m]; Q[2 * m + 1, 2 * m] = kCO[m]
        for m in range(n_modes - 1):
            a, b, c = sw[m, 0], sw[m, 1], sw[m, 2]
            d = sw[m, 3] if not db else (kOC[m] * c * kCO[m + 1] * b) / (kCO[m] * a * kOC[m + 1])
            Q[2 * m, 2 * m + 2] = a; Q[2 * m + 2, 2 * m] = b          # O_m <-> O_m+1
            Q[2 * m + 1, 2 * m + 3] = c; Q[2 * m + 3, 2 * m + 1] = d  # C_m <-> C_m+1
        return _finalize(Q)

    init = []
    for m in range(n_modes):                       # slow mode ~ (0.06, 0.25)/s, fast ~ (0.5, 1)/s
        f = m / max(n_modes - 1, 1)
        init += [np.log(0.06 * (8 ** f)), np.log(0.25 * (4 ** f))]
    for m in range(n_modes - 1):
        init += [np.log(1 / 900.), np.log(1 / 500.), np.log(1 / 900.)] + ([np.log(1 / 500.)] if not db else [])
    return Model(f"modal{n_modes}" + ("" if db else "_noDB"), classes, build, n_par, np.array(init))


def general(kO, kC):
    K = kO + kC; classes = [0] * kO + [1] * kC
    idx = [(i, j) for i in range(K) for j in range(K) if i != j]

    def build(r):
        Q = np.zeros((K, K))
        for (i, j), v in zip(idx, r): Q[i, j] = v
        return _finalize(Q)

    rng = np.random.default_rng(0)
    init = np.log(0.1) + rng.normal(0, 1.0, len(idx))
    return Model(f"general{kO}o{kC}c", classes, build, len(idx), init)


# ----------------------------------------------------------------------------- likelihood
def stationary(Q):
    """Stationary distribution by linear solve of pi Q = 0, sum(pi) = 1 (robust to eig ordering/sign)."""
    K = Q.shape[0]
    A = np.vstack([Q.T, np.ones(K)]); b = np.zeros(K + 1); b[-1] = 1
    pi, *_ = np.linalg.lstsq(A, b, rcond=None)
    pi = np.clip(pi, 0, None); return pi / pi.sum()


def _block_powers(Pcc, nmax):
    """(nmax+1, k, k) array of Pcc^n via eigendecomposition, falling back to cumulative products."""
    k = Pcc.shape[0]
    try:
        lam, V = np.linalg.eig(Pcc); Vi = np.linalg.inv(V)
        if np.linalg.cond(V) > 1e8: raise np.linalg.LinAlgError
        n = np.arange(nmax + 1)[:, None]
        lp = lam[None, :] ** n                                   # (nmax+1, k)
        return np.real(np.einsum('ij,nj,jl->nil', V, lp, Vi))
    except np.linalg.LinAlgError:
        out = np.empty((nmax + 1, k, k)); out[0] = np.eye(k)
        for n in range(1, nmax + 1): out[n] = out[n - 1] @ Pcc
        return out


def loglik(model, theta, dwells_list, block=32):
    """Vectorised version of loglik_slow: per-dwell matrices are built in batch and multiplied pairwise
    within blocks of `block` consecutive dwells; only the block products are applied sequentially."""
    Q = model.Q(theta)
    if not np.all(np.isfinite(Q)): return -np.inf
    P = expm(Q * DT); pi = stationary(Q)
    K = model.K; cls = model.classes
    idx = {c: np.flatnonzero(cls == c) for c in (0, 1)}
    total = 0.0
    for (dw, dc) in dwells_list:
        N = len(dw); nmax = int(dw.max())
        pw = {c: _block_powers(P[np.ix_(idx[c], idx[c])], nmax) for c in (0, 1)}
        M = np.zeros((N - 1, K, K))
        for c in (0, 1):
            sel = np.flatnonzero(dc[:-1] == c)
            if sel.size == 0: continue
            d = 1 - c
            Pcx = P[np.ix_(idx[c], idx[d])]
            M[np.ix_(sel, idx[c], idx[d])] = pw[c][dw[sel] - 1] @ Pcx
        pad = (-(N - 1)) % block
        if pad: M = np.concatenate([M, np.broadcast_to(np.eye(K), (pad, K, K))], 0)
        M = M.reshape(-1, block, K, K)
        while M.shape[1] > 1:
            M = M[:, 0::2] @ M[:, 1::2]
        M = M[:, 0]
        c0 = dc[0]; alpha = np.zeros(K); alpha[idx[c0]] = pi[idx[c0]]; alpha /= alpha.sum()
        ll = 0.0
        for Mb in M:
            alpha = alpha @ Mb; sm = alpha.sum()
            if not (sm > 0): return -np.inf
            ll += np.log(sm); alpha /= sm
        cN = dc[-1]; last = np.zeros(K); last[idx[cN]] = pw[cN][dw[-1] - 1] @ np.ones(len(idx[cN]))
        sm = alpha @ last
        if not (sm > 0): return -np.inf
        total += ll + np.log(sm)
    return total


def loglik_slow(model, theta, dwells_list, return_parts=False):
    Q = model.Q(theta)
    if not np.all(np.isfinite(Q)): return -np.inf
    P = expm(Q * DT); pi = stationary(Q)
    cls = model.classes
    idx = {c: np.flatnonzero(cls == c) for c in (0, 1)}
    Pcc = {c: P[np.ix_(idx[c], idx[c])] for c in (0, 1)}
    Pcx = {(c, d): P[np.ix_(idx[c], idx[d])] for c in (0, 1) for d in (0, 1) if c != d}
    total = 0.0
    for (dw, dc) in dwells_list:
        nmax = int(dw.max())
        # powers of the within-class blocks, computed by cumulative products (cheap: n_max x small matrices)
        pw = {}
        for c in (0, 1):
            k = len(idx[c]); arr = np.empty((nmax + 1, k, k)); arr[0] = np.eye(k)
            for n in range(1, nmax + 1): arr[n] = arr[n - 1] @ Pcc[c]
            pw[c] = arr
        c0 = dc[0]; alpha = pi[idx[c0]]; alpha = alpha / alpha.sum()
        ll = 0.0
        N = len(dw)
        for i in range(N):
            c, n = dc[i], dw[i]
            alpha = alpha @ pw[c][n - 1]
            if i < N - 1:
                alpha = alpha @ Pcx[(c, dc[i + 1])]
            sm = alpha.sum()
            if not (sm > 0): return -np.inf
            ll += np.log(sm); alpha = alpha / sm
        total += ll
    return total


def fit(model, dwells_list, n_restart=3, seed=0, maxiter=400, verbose=False):
    rng = np.random.default_rng(seed)
    best = None
    for r in range(n_restart):
        x0 = model.init + (rng.normal(0, 0.7, model.n_par) if r else 0)
        f = lambda th: -loglik(model, th, dwells_list)
        res = minimize(f, x0, method="L-BFGS-B", options=dict(maxiter=maxiter))
        if verbose: print(f"  restart {r}: ll={-res.fun:.2f}, nit={res.nit}", flush=True)
        if best is None or res.fun < best.fun: best = res
    return best


def bic(ll, n_par, n_dwells):
    return -2 * ll + n_par * np.log(n_dwells)


# ----------------------------------------------------------------------------- simulation
def simulate(model, theta, n_frames, seed=0):
    rng = np.random.default_rng(seed)
    Q = model.Q(theta); P = expm(Q * DT); pi = stationary(Q)
    C = np.cumsum(P, 1)
    x = rng.choice(model.K, p=pi); out = np.empty(n_frames, int)
    u = rng.random(n_frames)
    for t in range(n_frames):
        out[t] = x
        x = int(np.searchsorted(C[x], u[t]))
        if x >= model.K: x = model.K - 1
    return model.classes[out]
