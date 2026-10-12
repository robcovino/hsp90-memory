"""Summary statistics of a set of two-level 10-Hz traces (a 'measurement' = list of 0/1 arrays, e.g. 4 x 6 h).

All functions take a list of int arrays (0 = open, 1 = closed).  `summary_vector` concatenates everything into one
fixed-length vector with named entries, for classifiers / embedding-free SBI.
"""
import numpy as np
from .io import dwells, DT

LAGS = np.array([1, 2, 3, 5, 7, 10, 15, 20, 30, 50, 70, 100])
QUANTILES = np.array([0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999])
RATE_EDGES = np.arange(0, 1.1, 0.1)
CK_LAGS = np.unique(np.round(np.logspace(0, np.log10(6000), 25)).astype(int))   # frames: 0.1 s .. 600 s
LOG_EDGES = np.array([0, 0.2, 0.5, 1, 2, 5, 10, 30, 1e4])                         # s, for the 2D dwell histogram


def _rank(a):
    return np.argsort(np.argsort(a))


def spearman(a, b):
    ra, rb = _rank(a), _rank(b)
    return np.corrcoef(ra, rb)[0, 1]


def dwell_pairs(series_list):
    """Open dwells O_i and the closed dwell C_i that follows each, pooled over traces (censored ends dropped)."""
    O, C = [], []
    for s in series_list:
        d, st = dwells(s); d, st = d[1:-1] * DT, st[1:-1]
        if d.size < 4: continue
        if st[0] == 1: d, st = d[1:], st[1:]
        n = min((st == 0).sum(), (st == 1).sum())
        O.append(d[st == 0][:n]); C.append(d[st == 1][:n])
    return np.concatenate(O), np.concatenate(C)


def lag_correlations(O, C):
    oo = np.array([spearman(O[:-l], O[l:]) for l in LAGS])
    cc = np.array([spearman(C[:-l], C[l:]) for l in LAGS])
    return oo, cc, spearman(O, C), spearman(C[:-1], O[1:])


def dwell_quantiles(O, C):
    return np.log(np.quantile(O, QUANTILES)), np.log(np.quantile(C, QUANTILES))


def local_rate_hist(series_list, window_frames):
    r = np.concatenate([np.r_[np.abs(np.diff(s)), 0][: (s.size // window_frames) * window_frames]
                        .reshape(-1, window_frames).sum(1) / (window_frames * DT) for s in series_list])
    h, _ = np.histogram(np.clip(r, 0, RATE_EDGES[-1] - 1e-9), bins=RATE_EDGES)
    return h / h.sum()


def fano_within_modes(series_list, window_frames=600, mode_window=1200, thr=0.3):
    """Overdispersion of transition counts in 1-min windows, separately within slow and fast 2-min episodes.
    Markov-within-mode gives Fano ~ O(1); OU-modulated rates give Fano growing with sigma^2 tau_w."""
    out = []
    for m in (0, 1):
        counts = []
        for s in series_list:
            tr = np.r_[np.abs(np.diff(s)), 0]
            n = (s.size // mode_window) * mode_window
            rate = tr[:n].reshape(-1, mode_window).sum(1) / (mode_window * DT)
            mode = np.repeat((rate > thr).astype(int), mode_window)
            c = tr[:n].reshape(-1, window_frames).sum(1)
            mw = mode[:n].reshape(-1, window_frames).mean(1)
            counts.append(c[(mw > 0.99) if m else (mw < 0.01)])
        counts = np.concatenate(counts)
        out.append(counts.var() / counts.mean() if counts.size > 10 and counts.mean() > 0 else np.nan)
    return np.array(out)


def propagator(s, L):
    a, b = s[:-L], s[L:]; G = np.zeros((2, 2))
    for c0 in (0, 1):
        m = a == c0
        if m.sum() == 0: G[c0] = [0.5, 0.5]; continue
        G[c0, 1] = b[m].mean(); G[c0, 0] = 1 - G[c0, 1]
    return G


class _Counts:
    """Two-point counts of a 0/1 series at arbitrary lags from prefix sums and lagged dot products; gives the same
    propagators as `propagator` (exact integer counts, identical floating-point division)."""
    def __init__(self, s):
        self.s = np.ascontiguousarray(s, dtype=np.float64); self.n = self.s.size   # float dot: BLAS, exact < 2^53
        self.cs = np.r_[0, np.cumsum(np.asarray(s, dtype=np.int64))]

    def G(self, L):
        n, cs = self.n, self.cs
        n1 = int(cs[n - L]); n0 = (n - L) - n1                 # a = s[:-L]: frames in class 1 / 0
        b1 = int(cs[n] - cs[L])                              # b = s[L:]: frames in class 1
        n11 = int(np.dot(self.s[:-L], self.s[L:])); n01 = b1 - n11
        G = np.zeros((2, 2))
        for c0, (nn, k1) in enumerate(((n0, n01), (n1, n11))):
            if nn == 0: G[c0] = [0.5, 0.5]; continue
            G[c0, 1] = k1 / nn; G[c0, 0] = 1 - G[c0, 1]
        return G


def ck_divergence(s, tau_frames=10, lags=CK_LAGS):
    """Vollmar et al. (2024) Chapman-Kolmogorov divergence D_{tau,q0}(t), both q0, for one trace."""
    cnt = _Counts(s)
    Gtau = cnt.G(tau_frames); out = np.full((2, len(lags)), np.nan)
    for i, L in enumerate(lags):
        if L <= tau_frames: continue
        G = cnt.G(L); Gck = Gtau @ cnt.G(L - tau_frames)
        for q0 in (0, 1):
            out[q0, i] = np.sum(G[q0] * np.log(np.clip(G[q0], 1e-12, None) / np.clip(Gck[q0], 1e-12, None)))
    return out


def ck_profile(series_list, tau_frames=10):
    return np.nanmean([ck_divergence(s, tau_frames) for s in series_list], 0)


def dwell2d_spectrum(O, C, n_sv=4):
    """Singular values of the (open, following closed) log-binned joint histogram, normalised.
    For an aggregated Markov model the rank is bounded by the number of gateway states."""
    H, _, _ = np.histogram2d(O, C, bins=[LOG_EDGES, LOG_EDGES]); H = H / H.sum()
    sv = np.linalg.svd(H, compute_uv=False)
    return sv[:n_sv] / sv[0]


def acf(series_list, lags_frames=(1, 10, 50, 100, 300, 600, 1000, 2000, 3000)):
    out = []
    for s in series_list:
        x = s.astype(float) - s.mean(); v = (x * x).mean()
        out.append([np.mean(x[:-l] * x[l:]) / v for l in lags_frames])
    return np.mean(out, 0)


def summary_vector(series_list):
    """Fixed-length named summary vector."""
    O, C = dwell_pairs(series_list)
    oo, cc, oc, co = lag_correlations(O, C)
    qO, qC = dwell_quantiles(O, C)
    h1 = local_rate_hist(series_list, 600); h2 = local_rate_hist(series_list, 1200)
    fano = fano_within_modes(series_list)
    ck = ck_profile(series_list)
    sv = dwell2d_spectrum(O, C)
    ac = acf(series_list)
    closed_frac = np.mean([s.mean() for s in series_list])
    parts = [("closed_frac", np.array([closed_frac])), ("logqO", qO), ("logqC", qC),
             ("oo", oo), ("cc", cc), ("oc", np.array([oc, co])), ("rate1min", h1), ("rate2min", h2),
             ("fano", fano), ("ck_open", ck[0]), ("ck_closed", ck[1]), ("sv2d", sv), ("acf", ac)]
    names = [f"{k}_{i}" for k, v in parts for i in range(len(v))]
    return np.concatenate([v for _, v in parts]), names
