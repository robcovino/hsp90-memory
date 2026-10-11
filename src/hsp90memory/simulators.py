"""Simulators producing 10-Hz two-level snapshot series, vectorised over replicates.

simulate_amm        : continuous-time Markov chain with generator Q observed at DT through a class map (exact, expm).
simulate_modal_ou   : nested model for the identifiability study. Gating modes follow a discrete chain (rates taken
                      from a fitted modal AMM, switching allowed from both conformations); within each mode the
                      open<->closed log-rates are modulated by Ornstein-Uhlenbeck processes:
                          ln k_OC = ln k_OC(m) + sigma * z1,   ln k_CO = ln k_CO(m) + sigma * z2,
                      z1, z2 unit-variance OU with relaxation time tau_w and correlation rho.
                      sigma = 0 recovers the discrete model; sigma > 0 is 'continuous disorder within modes'.
                      The conformational update uses the exact two-state propagator over one frame given the current
                      rates (so multiple sub-frame transitions are included, as in the AMM likelihood).
"""
import numpy as np
from scipy.linalg import expm
from .io import DT


def stationary(Q):
    K = Q.shape[0]
    A = np.vstack([Q.T, np.ones(K)]); b = np.zeros(K + 1); b[-1] = 1
    pi, *_ = np.linalg.lstsq(A, b, rcond=None)
    pi = np.clip(pi, 0, None); return pi / pi.sum()


def simulate_amm(Q, classes, n_frames, n_rep=1, seed=0):
    """Returns an (n_rep, n_frames) int8 array of observed classes."""
    rng = np.random.default_rng(seed)
    Q = np.asarray(Q); classes = np.asarray(classes); K = Q.shape[0]
    P = expm(Q * DT); C = np.cumsum(P, 1); C[:, -1] = 1.0
    x = rng.choice(K, size=n_rep, p=stationary(Q))
    out = np.empty((n_rep, n_frames), np.int8)
    U = rng.random((n_frames, n_rep))
    for t in range(n_frames):
        out[:, t] = classes[x]
        x = (U[t][:, None] > C[x]).sum(1)
    return out


def modal_params_from_Q(Q):
    """Extract per-mode conformational rates and the mode-switch generators from a fitted modal AMM
    (states ordered O_0, C_0, O_1, C_1, ...)."""
    Q = np.asarray(Q); M = Q.shape[0] // 2
    kOC = np.array([Q[2 * m, 2 * m + 1] for m in range(M)])
    kCO = np.array([Q[2 * m + 1, 2 * m] for m in range(M)])
    SO = np.array([[Q[2 * m, 2 * n] if m != n else 0.0 for n in range(M)] for m in range(M)])       # mode switch from O
    SC = np.array([[Q[2 * m + 1, 2 * n + 1] if m != n else 0.0 for n in range(M)] for m in range(M)])  # from C
    return kOC, kCO, SO, SC


def simulate_modal_ou(kOC, kCO, SO, SC, sigma, tau_w, rho, n_frames, n_rep=1, seed=0, pi_mode=None):
    """See module docstring. Returns (n_rep, n_frames) int8 observed conformations (0 open, 1 closed)."""
    rng = np.random.default_rng(seed)
    kOC = np.asarray(kOC); kCO = np.asarray(kCO); M = kOC.size
    # per-frame mode-switch probabilities (first order; rates are << 1/DT)
    PO = SO * DT; PC = SC * DT
    CO_ = np.cumsum(PO, 1); CC_ = np.cumsum(PC, 1)      # cumulative rows have a flat step at the diagonal
    if pi_mode is None:
        # crude: stationary distribution of the mode chain averaged over conformations
        Qm = 0.5 * (SO + SC); np.fill_diagonal(Qm, -Qm.sum(1)); pi_mode = stationary(Qm)
    m = rng.choice(M, size=n_rep, p=pi_mode)
    c = (rng.random(n_rep) < kOC[m] / (kOC[m] + kCO[m])).astype(np.int8)
    z1 = rng.standard_normal(n_rep); z2p = rng.standard_normal(n_rep)
    a = np.exp(-DT / tau_w) if tau_w > 0 else 0.0; sq = np.sqrt(1 - a * a)
    out = np.empty((n_rep, n_frames), np.int8)
    U = rng.random((n_frames, 3, n_rep)); G = rng.standard_normal((n_frames, 2, n_rep))
    for t in range(n_frames):
        out[:, t] = c
        # OU modulators
        if sigma > 0:
            z1 = a * z1 + sq * G[t, 0]; z2p = a * z2p + sq * G[t, 1]
            z2 = rho * z1 + np.sqrt(1 - rho * rho) * z2p
            k1 = kOC[m] * np.exp(sigma * z1); k2 = kCO[m] * np.exp(sigma * z2)
        else:
            k1 = kOC[m]; k2 = kCO[m]
        # mode switching, decided on the conformation at the start of the frame
        u = U[t, 1]
        cum = np.where(c[:, None] == 0, CO_[m], CC_[m])           # (n_rep, M) cumulative switch probabilities
        sw = u < cum[:, -1]
        tgt = (u[:, None] > cum).sum(1)                             # the flat step at the diagonal is never selected
        m = np.where(sw, np.minimum(tgt, M - 1), m)
        # exact two-state propagator over one frame at the current rates
        kt = k1 + k2; pC = k1 / kt; relax = 1 - np.exp(-kt * DT)
        p_flip = np.where(c == 0, pC * relax, (1 - pC) * relax)
        c = np.where(U[t, 0] < p_flip, 1 - c, c).astype(np.int8)
    return out
