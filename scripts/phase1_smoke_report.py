"""Phase 1 smoke test report: distribution of the summary statistics over the prior-predictive simulations in
results/phase1/batch_*.npz, timings, non-finite entries, and where the real data fall (prior-predictive percentile).

Usage: python scripts/phase1_smoke_report.py [results/phase1]   ->  <dir>/smoke_report.json, figures/phase1_smoke.png
"""
import json, sys, warnings
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from hsp90memory import io, summaries as sm

d = Path(sys.argv[1] if len(sys.argv) > 1 else "results/phase1")
files = sorted(d.glob("batch_*.npz"))
X, P, T, rss = [], [], [], []
for f in files:
    with np.load(f) as z:
        X.append(z["x"]); P.append(z["params"]); T.append(z["sim_time"]); rss.append(float(z["peak_rss_mb"]))
        names = list(z["names"]); pnames = list(z["param_names"])
X, P, T = np.concatenate(X), np.concatenate(P), np.concatenate(T)

with warnings.catch_warnings():
    warnings.simplefilter("ignore", RuntimeWarning)
    obs, _ = sm.summary_vector([io.idealise(io.load_trace(f"data/{k}Thermal.txt")) for k in (1, 2, 3, 4)])

fin = np.isfinite(X)
q = np.nanquantile(X, [0.01, 0.05, 0.5, 0.95, 0.99], axis=0)
pct = np.array([np.mean(X[fin[:, j], j] < obs[j]) if np.isfinite(obs[j]) and fin[:, j].any() else np.nan
                for j in range(len(names))])
per_stat = {n: {"q01": q[0, j], "q05": q[1, j], "median": q[2, j], "q95": q[3, j], "q99": q[4, j],
                "frac_nonfinite": float(1 - fin[:, j].mean()), "data": float(obs[j]), "data_prior_pct": pct[j]}
            for j, n in enumerate(names)}
structural = [n for j, n in enumerate(names) if fin[:, j].mean() == 0]
partial = {n: float(1 - fin[:, j].mean()) for j, n in enumerate(names) if 0 < fin[:, j].mean() < 1}
outside = {n: pct[j] for j, n in enumerate(names) if np.isfinite(pct[j]) and (pct[j] < 0.005 or pct[j] > 0.995)}
rep = {"n_sim": int(len(X)), "n_batches": len(files), "sim_time_mean_s": float(T.mean()),
       "sim_time_q": np.quantile(T, [0.05, 0.5, 0.95]).tolist(), "peak_rss_mb_max": max(rss),
       "params_range": {n: [float(P[:, i].min()), float(P[:, i].max())] for i, n in enumerate(pnames)},
       "structural_nan": structural, "partial_nonfinite": partial, "data_outside_99pct_prior_predictive": outside,
       "per_stat": per_stat}
json.dump(rep, open(d / "smoke_report.json", "w"), indent=1, default=float)

groups = ["closed_frac", "logqO", "logqC", "oo", "cc", "oc", "rate1min", "rate2min", "fano", "ck_open", "ck_closed",
          "sv2d", "acf"]
fig, axs = plt.subplots(len(groups), 1, figsize=(11, 2.0 * len(groups)))
for ax, g in zip(axs, groups):
    idx = [j for j, n in enumerate(names) if n.rsplit("_", 1)[0] == g and fin[:, j].any()]
    data = [X[fin[:, j], j] for j in idx]
    ax.violinplot(data, showextrema=False, quantiles=[[0.05, 0.5, 0.95]] * len(data))
    ax.plot(range(1, len(idx) + 1), obs[idx], "o", color="C3", ms=4, label="data")
    ax.set_xticks(range(1, len(idx) + 1), [names[j].rsplit("_", 1)[1] for j in idx], fontsize=7)
    ax.set_ylabel(g, fontsize=8)
    if g.startswith("ck"): ax.set_yscale("symlog", linthresh=1e-4)
axs[0].legend(fontsize=8)
fig.suptitle(f"Phase 1 prior predictive, {len(X)} simulations (violins: 5/50/95 %), red: data", fontsize=10)
fig.tight_layout(); Path("figures").mkdir(exist_ok=True); fig.savefig("figures/phase1_smoke.png", dpi=110)

print(f"{len(X)} sims in {len(files)} batches; {T.mean():.3f} s/sim (5/50/95 %: {np.round(rep['sim_time_q'], 3)}); "
      f"peak RSS/worker {max(rss):.0f} MB")
print("structural NaN:", structural)
print("partially non-finite:", partial)
print("data outside central 99 % of prior predictive:", {k: round(v, 4) for k, v in outside.items()})
print(f"{'stat':14s} {'q05':>9s} {'median':>9s} {'q95':>9s} {'data':>9s} {'pct':>6s}")
for j, n in enumerate(names):
    s = per_stat[n]
    if n in structural: continue
    print(f"{n:14s} {s['q05']:9.4g} {s['median']:9.4g} {s['q95']:9.4g} {s['data']:9.4g} {s['data_prior_pct']:6.3f}")
