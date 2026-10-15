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
import numba
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


def ring(n_modes, db=True):
    """n_modes gating modes on a cycle (mode m <-> m+1 mod n_modes); states ordered (O_0, C_0, O_1, C_1, ...).
    Parameters: kOC_m, kCO_m for every mode, then for every edge m -> m+1: a_m (O_m -> O_m+1), b_m (O_m+1 -> O_m),
    c_m (C_m -> C_m+1); d_m (C_m+1 -> C_m) is fixed by detailed balance on the square (O_m, O_m+1, C_m+1, C_m).
    With db=True the last b (O_0 -> O_{n-1}) is fixed by detailed balance around the ring of open states
    (prod a = prod b); the ring of closed states then balances as well. db=False leaves it free (a driven cycle)."""
    K = 2 * n_modes; classes = [0, 1] * n_modes
    n_par = 2 * n_modes + 3 * n_modes - (1 if db else 0)

    def build(r):
        Q = np.zeros((K, K))
        kOC = r[0:2 * n_modes:2]; kCO = r[1:2 * n_modes:2]
        sw = list(r[2 * n_modes:])
        a, b, c = np.zeros(n_modes), np.zeros(n_modes), np.zeros(n_modes)
        for m in range(n_modes):
            a[m] = sw.pop(0)
            b[m] = sw.pop(0) if (m < n_modes - 1 or not db) else np.nan
            c[m] = sw.pop(0)
        if db:
            b[-1] = np.prod(a) / np.prod(b[:-1])
        for m in range(n_modes):
            Q[2 * m, 2 * m + 1] = kOC[m]; Q[2 * m + 1, 2 * m] = kCO[m]
        for m in range(n_modes):
            n = (m + 1) % n_modes
            d = kOC[m] * c[m] * kCO[n] * b[m] / (kCO[m] * a[m] * kOC[n])
            Q[2 * m, 2 * n] = a[m]; Q[2 * n, 2 * m] = b[m]
            Q[2 * m + 1, 2 * n + 1] = c[m]; Q[2 * n + 1, 2 * m + 1] = d
        return _finalize(Q)

    init = []
    for m in range(n_modes):
        f = m / max(n_modes - 1, 1); init += [np.log(0.06 * (8 ** f)), np.log(0.25 * (4 ** f))]
    for m in range(n_modes):
        init += [np.log(1 / 900.)] + ([np.log(1 / 500.)] if (m < n_modes - 1 or not db) else []) + [np.log(1 / 900.)]
    return Model(f"ring{n_modes}" + ("" if db else "_noDB"), classes, build, n_par, np.array(init))


def graph(n_modes, edges, db=True, name=None):
    """Gating modes on an arbitrary graph. edges: list of (i, j) mode pairs. Per edge (in the listed order): a (O_i ->
    O_j), b (O_j -> O_i), c (C_i -> C_j); d (C_j -> C_i) fixed by detailed balance on the square. The spanning tree is
    built from the edges in the listed order; for every other edge (one per independent cycle) b is fixed by detailed
    balance around the cycle of open states (db=True) or free (db=False). Parameters: kOC_m, kCO_m per mode, then per
    edge a, [b], c. graph(n, [(0,1),(1,2),...,(n-1,0)]) has the parametrisation of ring(n)."""
    K = 2 * n_modes; classes = [0, 1] * n_modes
    parent = list(range(n_modes))
    def find(x):
        while parent[x] != x: x = parent[x]
        return x
    tree = []
    for (i, j) in edges:
        ri, rj = find(i), find(j)
        tree.append(ri != rj)
        if ri != rj: parent[ri] = rj
    n_par = 2 * n_modes + sum(3 if (t or not db) else 2 for t in tree)

    def build(r):
        Q = np.zeros((K, K)); kOC = r[0:2 * n_modes:2]; kCO = r[1:2 * n_modes:2]; sw = list(r[2 * n_modes:])
        A, B, Cc = [], [], []
        for t in tree:
            A.append(sw.pop(0)); B.append(sw.pop(0) if (t or not db) else np.nan); Cc.append(sw.pop(0))
        if db:                                               # open-state potentials on the tree, then cycle closures
            phi = {0: 0.0}; todo = [k for k, t in enumerate(tree) if t]
            while todo:
                for k in list(todo):
                    i, j = edges[k]
                    if i in phi: phi[j] = phi[i] + np.log(A[k] / B[k]); todo.remove(k)
                    elif j in phi: phi[i] = phi[j] - np.log(A[k] / B[k]); todo.remove(k)
            for k, t in enumerate(tree):
                if not t:
                    i, j = edges[k]; B[k] = A[k] * np.exp(phi[i] - phi[j])
        for m in range(n_modes):
            Q[2 * m, 2 * m + 1] = kOC[m]; Q[2 * m + 1, 2 * m] = kCO[m]
        for k, (i, j) in enumerate(edges):
            d = kOC[i] * Cc[k] * kCO[j] * B[k] / (kCO[i] * A[k] * kOC[j])
            Q[2 * i, 2 * j] = A[k]; Q[2 * j, 2 * i] = B[k]; Q[2 * i + 1, 2 * j + 1] = Cc[k]; Q[2 * j + 1, 2 * i + 1] = d
        return _finalize(Q)

    init = np.r_[ring(n_modes).init[:2 * n_modes],
                 np.concatenate([[np.log(1 / 900.)] + ([np.log(1 / 500.)] if (t or not db) else []) + [np.log(1 / 900.)]
                                 for t in tree])]
    m = Model(name or f"graph{n_modes}", classes, build, n_par, init); m.edges, m.tree = edges, tree
    return m


RING4 = [(0, 1), (1, 2), (2, 3), (3, 0)]
NAMED_GRAPHS = {"ringL1": (5, RING4 + [(3, 4)]), "ringL2": (6, RING4 + [(3, 4), (4, 5)]),
                "ringS2": (6, RING4 + [(3, 4), (3, 5)]),
                "ringF1": (5, RING4 + [(0, 4)]),       # ring + a mode attached to the frozen mode
                "full5": (5, RING4 + [(0, 4), (0, 2), (1, 3), (1, 4), (2, 4), (3, 4)])}   # all 10 pairs (ceiling)


def named_graph(name, db=True):
    """Discrete competitors on the ring skeleton (frozen 0, slow 1, open-biased 2, fast 3): a kOC ladder branching off
    the fast mode, as a chain (ringL1, ringL2) or a star (ringS2)."""
    n, e = NAMED_GRAPHS[name]; return graph(n, e, db=db, name=name)


def ou_grid(n_grid=15, z_max=3.5):
    """Unit-variance Ornstein-Uhlenbeck coordinate on a grid, square-root approximation of the Fokker-Planck
    operator for U(z) = z^2/2. Returns (z, R1): z rescaled so that the discrete stationary distribution has unit
    variance, and the generator R1 rescaled so that its slowest relaxation time is exactly 1 (multiply by 1/tau)."""
    z = np.linspace(-z_max, z_max, n_grid); h = z[1] - z[0]; U = z ** 2 / 2
    R = np.zeros((n_grid, n_grid))
    for i in range(n_grid - 1):
        R[i, i + 1] = np.exp(-(U[i + 1] - U[i]) / 2) / h ** 2
        R[i + 1, i] = np.exp(-(U[i] - U[i + 1]) / 2) / h ** 2
    R = _finalize(R)
    pz = stationary(R); z = z / np.sqrt(np.sum(pz * z ** 2))
    ev = np.sort(np.linalg.eigvals(R).real)[::-1]
    return z, R / (-ev[1])


def split_mode(base, mode, k=2, driven=False):
    """Replace gating mode `mode` of `base` (states O_m = 2m, C_m = 2m+1) by k identical sub-states (same kOC, kCO),
    appended as new modes M, M+1, ... Extra parameter: log h (hopping rate between successive sub-states).
    driven=False (equilibrium): every connection to other modes stays on the first copy; copies hop symmetrically
      (rate h both ways, in O and in C), so detailed balance is kept. The mode lifetime becomes a mixture of
      exponentials with a heavier tail (copies act as traps); it cannot become more regular in equilibrium.
    driven=True: the mode is entered at the first copy and left from the last; copies are traversed forwards only
      (rate h, in O and in C). The lifetime becomes Erlang-like (more regular). Breaks detailed balance."""
    M = base.K // 2; nb = base.n_par; K = 2 * (M + k - 1)
    copies = [mode] + [M + j for j in range(k - 1)]
    classes = list(base.classes) + [0, 1] * (k - 1)

    def build(theta):
        Q0 = base.Q(theta[:nb]).copy(); np.fill_diagonal(Q0, 0); h = np.exp(theta[nb])
        Q = np.zeros((K, K)); Q[:base.K, :base.K] = Q0
        for c in copies[1:]:
            Q[2 * c, 2 * c + 1] = Q0[2 * mode, 2 * mode + 1]; Q[2 * c + 1, 2 * c] = Q0[2 * mode + 1, 2 * mode]
        if driven:                                    # move exits to the last copy
            last = copies[-1]
            for n in range(M):
                if n == mode: continue
                for x in (0, 1):
                    Q[2 * last + x, 2 * n + x] = Q0[2 * mode + x, 2 * n + x]; Q[2 * mode + x, 2 * n + x] = 0.0
        for a_, b_ in zip(copies[:-1], copies[1:]):
            for x in (0, 1):
                Q[2 * a_ + x, 2 * b_ + x] = h
                if not driven: Q[2 * b_ + x, 2 * a_ + x] = h
        return _finalize(Q)

    m = Model(f"{base.name}_split{mode}x{k}{'drv' if driven else 'eq'}", classes, None, nb + 1, np.r_[base.init, np.log(1 / 300.)])
    m.Q = build; m.base_model, m.split = base, (mode, k, driven)
    return m


def modal_ou(n_modes, n_grid=15, z_max=3.5, couple="OC", topology="chain", base=None):
    """modal(n_modes) with the open->closed rate (couple='OC'), or both rates (couple='both'), modulated by one
    unit-variance OU coordinate z on a grid (the joint chain is a finite AMM; the gauge is fixed by unit variance and
    relaxation time tau_w):  kOC_m(z) = kOC_m exp(sigma_OC z),  kCO_m(z) = kCO_m exp(sigma_CO z).
    Mode switching does not depend on z and z does not depend on (mode, conformation): the same model as
    simulators.simulate_modal_ou with rho = +-1 (or sigma_CO = 0). Not in detailed balance for sigma > 0.
    couple='eq': equilibrium version. Same rates kOC_m exp(sigma_OC z), but z relaxes in U_O(z) = z^2/2 while open
    and in U_C(z) = U_O(z) - sigma_OC z while closed (the closed state tilts the slow coordinate), so that every
    cycle of the joint chain is balanced (detailed balance holds if it holds for modal(n_modes)). Parameters as 'OC'.
    Parameters: the modal(n_modes) log-rates, then log sigma_OC, log tau_w (couple='OC'), or log sigma_OC,
    sigma_CO (signed), log tau_w (couple='both'). States ordered (z_0: O_0, C_0, O_1, ...), (z_1: ...), ..."""
    if base is not None:
        topology = "custom"; n_modes = base.K // 2
    elif topology in NAMED_GRAPHS:
        base = named_graph(topology); n_modes = NAMED_GRAPHS[topology][0]
    else:
        base = {"chain": lambda: modal(n_modes), "ring": lambda: ring(n_modes),
                "ring_noDB": lambda: ring(n_modes, db=False)}[topology]()
    K2 = base.K; nb = base.n_par
    z, R1 = ou_grid(n_grid, z_max)
    classes = list(base.classes) * n_grid
    O = np.arange(0, K2, 2); C = O + 1

    def build_from(theta):
        Q0 = base.Q(theta[:nb]).copy(); np.fill_diagonal(Q0, 0)
        if couple in ("OC", "eq"):
            sOC, sCO, tau = np.exp(theta[nb]), 0.0, np.exp(theta[nb + 1])
        else:
            sOC, sCO, tau = np.exp(theta[nb]), theta[nb + 1], np.exp(theta[nb + 2])
        if couple == "eq":
            S = np.sqrt(np.clip(R1 * R1.T, 0, None))                     # symmetric part of the SQRA generator
            UO = -np.log(stationary(R1)); UC = UO - sOC * z
            RC = S * np.exp(-(UC[None, :] - UC[:, None]) / 2); np.fill_diagonal(RC, 0); RC = _finalize(RC)
            PO = np.zeros((K2, K2)); PO[O, O] = 1; PC = np.zeros((K2, K2)); PC[C, C] = 1
            Q = np.kron(R1 / tau, PO) + np.kron(RC / tau, PC)
        else:
            Q = np.kron(R1 / tau, np.eye(K2))
        for i in range(n_grid):
            Qi = Q0.copy()
            Qi[O, C] *= np.exp(sOC * z[i]); Qi[C, O] *= np.exp(sCO * z[i])
            Q[i * K2:(i + 1) * K2, i * K2:(i + 1) * K2] += Qi
        return _finalize(Q)

    n_extra = 3 if couple == "both" else 2
    init = np.r_[base.init, [np.log(0.3), 0.0, np.log(100.0)] if couple == "both" else [np.log(0.3), np.log(100.0)]]
    suffix = {"OC": "", "both": "2", "eq": "eq"}[couple]
    if topology == "custom":
        mname = f"{base.name}_ou{suffix}_g{n_grid}"
    elif topology in NAMED_GRAPHS:
        mname = f"{topology}_ou{suffix}_g{n_grid}"
    else:
        pre = "modal" if topology == "chain" else topology.replace("_noDB", "")
        mname = f"{pre}{n_modes}{'_noDB' if topology == 'ring_noDB' else ''}_ou{suffix}_g{n_grid}"
    m = Model(mname, classes, None, nb + n_extra, init)
    m.Q = build_from
    m.base, m.n_grid, m.z, m.couple = base, n_grid, z, couple
    return m


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


@numba.njit(cache=True)
def _forward_eig(dw, dc, alpha0, lam0, V0, W0, u0, lam1, V1, W1, u1):
    """Forward recursion with Pcc^n = V diag(lam^n) V^-1 for both classes; W_c = V_c^-1 P_cd, u_c = V_c^-1 1."""
    terms = np.empty(dw.size)
    return _forward_eig_terms(dw, dc, alpha0, lam0, V0, W0, u0, lam1, V1, W1, u1, terms)


@numba.njit(cache=True)
def _forward_eig_terms(dw, dc, alpha0, lam0, V0, W0, u0, lam1, V1, W1, u1, terms):
    """As _forward_eig; also writes log p(dwell i | dwells < i) into terms (the last term includes censoring)."""
    N = dw.size; ll = 0.0
    alpha = alpha0.astype(np.complex128)
    for i in range(N):
        c = dc[i]; n = dw[i] - 1
        if c == 0:
            lam, V, W, u = lam0, V0, W0, u0
        else:
            lam, V, W, u = lam1, V1, W1, u1
        k = lam.size
        beta = np.zeros(k, np.complex128)
        for j in range(k):
            acc = 0j
            for l in range(k):
                acc += alpha[l] * V[l, j]
            beta[j] = acc * lam[j] ** n
        if i < N - 1:
            kd = W.shape[1]; new = np.zeros(kd, np.complex128); s = 0.0
            for j in range(kd):
                acc = 0j
                for l in range(k):
                    acc += beta[l] * W[l, j]
                new[j] = acc.real; s += acc.real
            if not (s > 0):
                return -np.inf
            terms[i] = np.log(s); ll += terms[i]; alpha = new / s
        else:
            acc = 0j
            for l in range(k):
                acc += beta[l] * u[l]
            if not (acc.real > 0):
                return -np.inf
            terms[i] = np.log(acc.real); ll += terms[i]
    return ll


def _eig_class(P, ic, id_):
    Pcc = P[np.ix_(ic, ic)]
    lam, V = np.linalg.eig(Pcc)
    if np.linalg.cond(V) > 1e10:
        raise np.linalg.LinAlgError("ill-conditioned eigenvectors")
    Vi = np.linalg.inv(V)
    return (lam.astype(np.complex128), np.ascontiguousarray(V, np.complex128),
            np.ascontiguousarray(Vi @ P[np.ix_(ic, id_)], np.complex128),
            np.ascontiguousarray(Vi @ np.ones(len(ic)), np.complex128))


def loglik_Q(Q, classes, dwells_list):
    """Exact log-likelihood of the dwell sequences for generator Q (any size), same as `loglik` but with an
    O(k^2)-per-dwell forward pass (numba) on the eigendecomposition of the within-class blocks of P = expm(Q dt).
    Falls back to `loglik_slow`-style cumulative powers if an eigenvector basis is ill-conditioned."""
    Q = np.asarray(Q, float)
    if not np.all(np.isfinite(Q)): return -np.inf
    classes = np.asarray(classes)
    P = expm(Q * DT); pi = stationary(Q)
    idx = [np.flatnonzero(classes == c) for c in (0, 1)]
    try:
        e0 = _eig_class(P, idx[0], idx[1]); e1 = _eig_class(P, idx[1], idx[0])
    except np.linalg.LinAlgError:
        return _loglik_powers(P, pi, idx, dwells_list)
    total = 0.0
    for dw, dc in dwells_list:
        c0 = dc[0]; a0 = pi[idx[c0]]; s = a0.sum()
        if not (s > 0): return -np.inf
        ll = _forward_eig(np.asarray(dw, np.int64), np.asarray(dc, np.int64), a0 / s, *e0, *e1)
        if not np.isfinite(ll): return -np.inf
        total += ll
    return total


def loglik_terms(Q, classes, dwells_list):
    """Per-dwell predictive log-likelihood terms log p(dwell_i | earlier dwells of the trace), one array per trace;
    the first term of a trace conditions on its first-frame class. Sums equal loglik_Q."""
    Q = np.asarray(Q, float); classes = np.asarray(classes)
    P = expm(Q * DT); pi = stationary(Q)
    idx = [np.flatnonzero(classes == c) for c in (0, 1)]
    e0 = _eig_class(P, idx[0], idx[1]); e1 = _eig_class(P, idx[1], idx[0])
    out = []
    for dw, dc in dwells_list:
        a0 = pi[idx[dc[0]]]; t = np.empty(len(dw))
        _forward_eig_terms(np.asarray(dw, np.int64), np.asarray(dc, np.int64), a0 / a0.sum(), *e0, *e1, t)
        out.append(t)
    return out


def _loglik_powers(P, pi, idx, dwells_list):
    total = 0.0
    for dw, dc in dwells_list:
        nmax = int(dw.max()); pw = {}
        for c in (0, 1):
            Pcc = P[np.ix_(idx[c], idx[c])]; k = len(idx[c])
            arr = np.empty((nmax + 1, k, k)); arr[0] = np.eye(k)
            for n in range(1, nmax + 1): arr[n] = arr[n - 1] @ Pcc
            pw[c] = arr
        alpha = pi[idx[dc[0]]]; alpha = alpha / alpha.sum(); ll = 0.0
        for i in range(len(dw)):
            c = dc[i]; alpha = alpha @ pw[c][dw[i] - 1]
            if i < len(dw) - 1: alpha = alpha @ P[np.ix_(idx[c], idx[1 - c])]
            sm = alpha.sum()
            if not (sm > 0): return -np.inf
            ll += np.log(sm); alpha = alpha / sm
        total += ll
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
