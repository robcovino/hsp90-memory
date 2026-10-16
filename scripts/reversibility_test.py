"""Model-free test of time-reversal symmetry (detailed balance) on the open/closed dwell sequence.

In equilibrium the probability of a dwell sequence equals that of the reversed sequence. For an alternating two-state
record this means: the joint law of (open dwell, the closed dwell k places AFTER it) equals the joint law of (open
dwell, the closed dwell k places BEFORE it). Statistics (interior dwells, pooled over traces):
  - D_k = Spearman(O_i, C_after_k) - Spearman(O_i, C_before_k), k = 0..K (k = 0: the adjacent closed dwells);
  - J = symmetrised KL divergence between the log-binned 2D histograms of (O, next C) and (O, previous C).
Null distribution: the same statistics on data sets simulated from the final model (ringF1 + equilibrium OU), which
satisfies detailed balance exactly, 200 posterior draws (scripts/ppc_paper_figure.py simulator).
Output: results/reversibility/reversibility.json
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import json, sys
from multiprocessing import get_context
from pathlib import Path
import numpy as np
from scipy.linalg import expm
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ppc_paper_figure import _sim, model_of            # noqa: E402
from hsp90memory import amm, io, summaries as sm      # noqa: E402

K = 10
EDGES = np.array([0, 0.2, 0.5, 1, 2, 5, 10, 30, 1e4])


def pairs(series):
    """Per trace: open dwells O_j and closed dwells C_j with C_j immediately after O_j (alternation), interior only."""
    out = []
    for s in series:
        d, c = io.dwells(s); d, c = d[1:-1] * io.DT, c[1:-1]
        if c[0] == 1: d, c = d[1:], c[1:]
        n = min((c == 0).sum(), (c == 1).sum()); out.append((d[c == 0][:n], d[c == 1][:n]))
    return out


def stats(series):
    P = pairs(series); D = []
    for k in range(K + 1):
        oa, ca, ob, cb = [], [], [], []
        for O, C in P:
            n = len(O)
            oa.append(O[:n - k]); ca.append(C[k:])             # C_j+k is k places after the closed dwell right after O_j
            ob.append(O[k + 1:]); cb.append(C[:n - k - 1])     # C_j-1-k is k places before the closed dwell right before O_j
        D.append(sm.spearman(np.concatenate(oa), np.concatenate(ca)) - sm.spearman(np.concatenate(ob), np.concatenate(cb)))
    O = np.concatenate([p[0][1:] for p in P]); Ca = np.concatenate([p[1][1:] for p in P]); Cb = np.concatenate([p[1][:-1] for p in P])
    Ha, _, _ = np.histogram2d(O, Ca, bins=[EDGES, EDGES]); Hb, _, _ = np.histogram2d(O, Cb, bins=[EDGES, EDGES])
    pa = (Ha + 0.5) / (Ha + 0.5).sum(); pb = (Hb + 0.5) / (Hb + 0.5).sum()
    J = float(np.sum((pa - pb) * np.log(pa / pb)))
    return np.array(D), J, (Ha - Hb)


def one(args):
    tag, phi, seed = args
    m, to_theta = model_of(tag); Q = m.Q(to_theta(phi))
    P = expm(Q * io.DT); C = np.cumsum(P, 1); C[:, -1] = 1.0; pi = amm.stationary(Q); cls = np.asarray(m.classes)
    rng = np.random.default_rng(seed)
    series = [cls[_sim(C, rng.choice(pi.size, p=pi), rng.random(216000))].astype(np.int8) for _ in range(4)]
    D, J, _ = stats(series); return D, J


def main():
    data = [io.idealise(io.load_trace(f"data/{k}Thermal.txt")) for k in "1234"]
    D0, J0, dH = stats(data)
    per_trace = [stats([s])[0][:3].tolist() for s in data]
    phi = np.load("results/exact/thin/bayes_ringF1_oueq_thin.npz")["phi"]
    rng = np.random.default_rng(7); idx = rng.choice(len(phi), 200, replace=False)
    with get_context("fork").Pool(2) as pool:
        R = pool.map(one, [("ringF1_oueq", phi[i], 5000 + j) for j, i in enumerate(idx)], chunksize=4)
    Ds = np.array([r[0] for r in R]); Js = np.array([r[1] for r in R])
    z = (D0 - Ds.mean(0)) / Ds.std(0)
    pD = [float(np.mean(np.abs(Ds[:, k] - Ds[:, k].mean()) >= abs(D0[k] - Ds[:, k].mean()))) for k in range(K + 1)]
    chi = float(np.sum(z ** 2)); chis = np.array([np.sum(((d - Ds.mean(0)) / Ds.std(0)) ** 2) for d in Ds])
    res = dict(K=K, D_data=D0.tolist(), D_null_mean=Ds.mean(0).tolist(), D_null_sd=Ds.std(0).tolist(), z=z.tolist(),
               p_two_sided=pD, sum_z2=chi, p_sum_z2=float(np.mean(chis >= chi)), J_data=J0,
               J_null_median=float(np.median(Js)), p_J=float(np.mean(Js >= J0)), per_trace_D012=per_trace,
               hist_after_minus_before=dH.tolist(), edges_s=EDGES.tolist(), n_null=len(R),
               null_model="ringF1_oueq posterior (detailed balance), results/exact/thin/bayes_ringF1_oueq_thin.npz")
    Path("results/reversibility").mkdir(parents=True, exist_ok=True)
    json.dump(res, open("results/reversibility/reversibility.json", "w"), indent=1)
    print("D_k data :", np.round(D0, 3)); print("null sd  :", np.round(Ds.std(0), 3)); print("z        :", np.round(z, 1))
    print(f"sum z^2 = {chi:.1f} (p = {res['p_sum_z2']:.3f});  J = {J0:.4f}, null median {res['J_null_median']:.4f}, p = {res['p_J']:.3f}")
    print("per trace D_0..2:", np.round(per_trace, 3).tolist())
    print("after - before counts (rows O bins, cols C bins):\n", dH.astype(int))


if __name__ == "__main__":
    main()
