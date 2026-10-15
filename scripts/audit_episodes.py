"""Audit: where in time does the OU gain (ringF1 + OU eq. vs ringF1) come from, and what do those episodes look like?
1. Per-dwell predictive log-lik differences (amm.loglik_terms), summed in 30-min blocks per trace.
2. For every block: dwell counts, mean open/closed dwell, number of 1-frame dwells, local closing rate
   (O->C transitions per open second), and the raw-intensity levels and widths of the open and closed peaks (I/<I>),
   to look for instrumental anomalies (level shifts, broadened peaks) in the blocks that carry the gain.
3. Smoothed hidden coordinate E[z | all data](t) on a 1-s grid (frame-level forward-backward on the grid-15 joint
   chain), its autocorrelation, and the zero-lag cross-correlation of 2-min averages between traces.
Output: results/audit/episodes.json, results/audit/episodes_z.npz.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json
from pathlib import Path
import numpy as np, numba
from scipy.linalg import expm
from hsp90memory import amm, io

OUT = Path("results/audit"); BLOCK = 18000          # frames per 30-min block
th1 = np.array(json.load(open("results/exact/fit_ringF1.json"))["theta"])
th2 = np.array(json.load(open("results/exact/fit_ringF1_oueq_global.json"))["theta"])
m1 = amm.named_graph("ringF1"); m2 = amm.modal_ou(5, n_grid=15, couple="eq", topology="ringF1")
Q1, Q2 = m1.Q(th1), m2.Q(th2)


@numba.njit(cache=True)
def fb_z(s, P, cls, pi, zstate, every):
    """Scaled forward-backward on the frame sequence; returns E[z | data] at frames 0, every, 2*every, ..."""
    T = s.size; K = P.shape[0]; n_out = (T + every - 1) // every
    A = np.empty((n_out, K)); a = np.zeros(K)
    for i in range(K):
        if cls[i] == s[0]: a[i] = pi[i]
    a /= a.sum(); A[0] = a
    for t in range(1, T):
        b = np.zeros(K)
        for j in range(K):
            if cls[j] != s[t]: continue
            acc = 0.0
            for i in range(K):
                if a[i] != 0.0: acc += a[i] * P[i, j]
            b[j] = acc
        a = b / b.sum()
        if t % every == 0: A[t // every] = a
    out = np.empty(n_out); be = np.ones(K)
    for i in range(K):
        if cls[i] != s[T - 1]: be[i] = 0.0
    for t in range(T - 1, -1, -1):
        if t % every == 0:
            w = A[t // every] * be; w /= w.sum(); out[t // every] = (w * zstate).sum()
        if t == 0: break
        nb = np.zeros(K)
        for i in range(K):
            if cls[i] != s[t - 1]: continue
            acc = 0.0
            for j in range(K):
                if be[j] != 0.0: acc += P[i, j] * be[j]
            nb[i] = acc
        be = nb / nb.sum()
    return out


def peak_stats(I, s):
    o, c = I[s == 0], I[s == 1]
    return dict(I_open=float(np.median(o)), I_closed=float(np.median(c)),
                w_open=float(1.4826 * np.median(np.abs(o - np.median(o)))) if o.size else np.nan,
                w_closed=float(1.4826 * np.median(np.abs(c - np.median(c)))) if c.size else np.nan)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    traces = [io.load_trace(f"data/{k}Thermal.txt") for k in "1234"]
    S = [io.idealise(I) for I in traces]; dl = [io.dwells(s) for s in S]
    t1 = amm.loglik_terms(Q1, m1.classes, dl); t2 = amm.loglik_terms(Q2, m2.classes, dl)
    blocks = []
    for k, ((dw, dc), a, b, I, s) in enumerate(zip(dl, t1, t2, traces, S)):
        end = np.cumsum(dw); start = end - dw; blk = np.minimum((end - 1) // BLOCK, 11); d = b - a
        tr_peaks = peak_stats(I, s)
        for j in range(12):
            sel = blk == j; fr = slice(j * BLOCK, (j + 1) * BLOCK)
            o = sel & (dc == 0); c = sel & (dc == 1)
            n_oc = int(np.sum(np.diff(s[fr].astype(int)) == 1))
            ps = peak_stats(I[fr], s[fr])
            blocks.append(dict(trace=k + 1, block=j, t0_h=j * 0.5, gain=float(d[sel].sum()), n_dwells=int(sel.sum()),
                               mean_open_s=float(dw[o].mean() * io.DT), mean_closed_s=float(dw[c].mean() * io.DT),
                               n1_open=int(np.sum(o & (dw == 1))), n1_closed=int(np.sum(c & (dw == 1))),
                               k_close=float(n_oc / (np.sum(s[fr] == 0) * io.DT)),
                               open_frac=float(np.mean(s[fr] == 0)),
                               dI_open=ps["I_open"] - tr_peaks["I_open"], dI_closed=ps["I_closed"] - tr_peaks["I_closed"],
                               w_open=ps["w_open"], w_closed=ps["w_closed"]))
    gains = np.array([b["gain"] for b in blocks]); order = np.argsort(-gains)
    print(f"total gain {gains.sum():.1f}; top-5 share {gains[order[:5]].sum() / gains.sum():.2f}; "
          f"blocks with gain < 0: {(gains < 0).sum()} of {gains.size}")
    keys = ["trace", "t0_h", "gain", "n_dwells", "mean_open_s", "mean_closed_s", "k_close", "open_frac", "n1_open",
            "n1_closed", "dI_open", "dI_closed", "w_open", "w_closed"]
    print(" ".join(f"{k:>9s}" for k in keys))
    for i in list(order[:8]) + list(order[-4:]):
        print(" ".join(f"{blocks[i][k]:9.4g}" for k in keys))
    W = np.array([[b["w_open"], b["w_closed"], b["dI_open"], b["dI_closed"], b["k_close"], b["mean_open_s"]] for b in blocks])
    print("medians over all blocks (w_open, w_closed, dI_open, dI_closed, k_close, mean_open):", np.round(np.median(W, 0), 4))
    # smoothed z
    P2 = expm(Q2 * io.DT); pi2 = amm.stationary(Q2); cls2 = np.asarray(m2.classes, np.int64)
    zst = np.repeat(m2.z, m2.base.K); Z = []
    for s in S:
        Z.append(fb_z(s.astype(np.int64), P2, cls2, pi2, zst, 10))
    Z = np.array(Z)                                            # (4, 21600), 1-s grid
    acf = []
    for z in Z:
        x = z - z.mean(); acf.append([float(np.dot(x[:-L or None], x[L:]) / np.dot(x, x)) for L in (0, 60, 120, 270, 600, 1200)])
    Zb = Z.reshape(4, -1, 120).mean(2)                         # 2-min averages
    cc = np.corrcoef(Zb)
    print("E[z] per trace: mean", np.round(Z.mean(1), 3), "sd", np.round(Z.std(1), 3))
    print("ACF of E[z] at 0, 1, 2, 4.5, 10, 20 min:", np.round(acf, 2).tolist())
    print("zero-lag correlation of 2-min E[z] between traces:\n", np.round(cc, 2))
    # trend: correlation of E[z] with time, and first-half vs second-half means
    tt = np.arange(Z.shape[1])
    print("corr(E[z], t):", np.round([np.corrcoef(z, tt)[0, 1] for z in Z], 2),
          " half means:", np.round([[z[:z.size // 2].mean(), z[z.size // 2:].mean()] for z in Z], 2).tolist())
    json.dump(dict(blocks=blocks, total_gain=float(gains.sum()), top5_share=float(gains[order[:5]].sum() / gains.sum()),
                   z_acf=acf, z_cross_corr_2min=cc.tolist(), z_mean=Z.mean(1).tolist(), z_sd=Z.std(1).tolist()),
              open(OUT / "episodes.json", "w"), indent=1)
    np.savez_compressed(OUT / "episodes_z.npz", Ez=Z, dt_s=1.0)


if __name__ == "__main__":
    main()
