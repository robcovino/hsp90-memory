"""Phase 1: posterior-predictive checks on the real traces.

modal-OU: n draws from the NPE posterior for the data (npe_<tag>_posterior_data.npz), one simulated measurement
(4 x 6 h) each. Reference: modal3 and modal4 at their ML generators (results/modal3.json, modal4.json), n measurements
each with the exact AMM simulator (plug-in predictive). For every model: summaries.summary_vector, dwell survival
curves on a fixed time grid; for the data the same. Per statistic: predictive percentile of the data and
z = (data - mean) / sd.

Usage: python scripts/phase1_ppc.py --tag v1 [--n 200] [--workers 4]
Outputs: results/phase1/ppc_<tag>.json, results/phase1/ppc_<tag>.npz, figures/phase1_ppc_<tag>.png
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, json, sys, time, warnings
from multiprocessing import get_context
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from hsp90memory import amm, io, simulators as sim, summaries as sm

NF, NT = 216000, 4
T_GRID = np.logspace(-1, 3.3, 44)     # s
D = Path("results/phase1")


def survival_grid(series_list, cls):
    d = np.concatenate([(lambda dw, dc: dw[1:-1][dc[1:-1] == cls])(*io.dwells(s)) for s in series_list]) * io.DT
    return np.array([(d > t - 1e-9).mean() for t in T_GRID])


def describe(Y):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        v, names = sm.summary_vector(Y)
    return v, survival_grid(Y, 0), survival_grid(Y, 1)


def job_ou(args):
    i, p = args
    model = amm.modal(3); theta, (sOC, sCO, l10tau, rho) = p[:12], p[12:]
    kOC, kCO, SO, SC = sim.modal_params_from_Q(model.Q(theta))
    Y = sim.simulate_modal_ou_nb(kOC, kCO, SO, SC, (sOC, sCO), 10 ** l10tau, rho, NF, n_rep=NT, seed=900000 + i)
    return describe(list(Y))


def job_amm(args):
    key, i = args
    r = json.load(open(f"results/{key}.json"))
    Y = sim.simulate_amm(np.array(r["Q"]), r["classes"], NF, n_rep=NT, seed=700000 + 1000 * len(key) + i)
    return describe(list(Y))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--dir", default="results/phase1", help="where npe_<tag>_posterior_data.npz is and outputs go")
    ap.add_argument("--figdir", default="figures")
    a = ap.parse_args()
    D = Path(a.dir)
    out_js = D / f"ppc_{a.tag}.json"
    if out_js.exists():
        sys.exit(f"{out_js} exists")
    t0 = time.time()
    z = np.load(D / f"npe_{a.tag}_posterior_data.npz")
    rng = np.random.default_rng(a.seed)
    draws = z["samples"][rng.choice(len(z["samples"]), a.n, replace=False)]
    data = [io.idealise(io.load_trace(f"data/{k}Thermal.txt")) for k in (1, 2, 3, 4)]
    obs = describe(data)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        names = sm.summary_vector(data)[1]
    res = {}
    with get_context("fork").Pool(a.workers) as pool:
        res["modal-OU"] = pool.map(job_ou, list(enumerate(draws)), chunksize=4)
        print(f"modal-OU done {time.time() - t0:.0f} s", flush=True)
        for key in ("modal3", "modal4"):
            res[key] = pool.map(job_amm, [(key, i) for i in range(a.n)], chunksize=4)
            print(f"{key} done {time.time() - t0:.0f} s", flush=True)

    report = {"tag": a.tag, "n": a.n, "time_s": None, "models": {}}
    arrays = {"names": np.array(names), "t_grid": T_GRID, "obs_v": obs[0], "obs_sO": obs[1], "obs_sC": obs[2]}
    groups = sorted({n.rsplit("_", 1)[0] for n in names})
    for key, R in res.items():
        V = np.array([r[0] for r in R]); SO = np.array([r[1] for r in R]); SCc = np.array([r[2] for r in R])
        arrays[f"{key}_v"] = V; arrays[f"{key}_sO"] = SO; arrays[f"{key}_sC"] = SCc
        ok = np.isfinite(obs[0]) & np.isfinite(V).all(0) & (np.nanstd(V, 0) > 0)
        zz = np.full(len(names), np.nan); pc = np.full(len(names), np.nan)
        zz[ok] = (obs[0][ok] - V[:, ok].mean(0)) / V[:, ok].std(0)
        pc[ok] = (V[:, ok] < obs[0][ok]).mean(0)
        okv = ok.sum()
        outside95 = int(np.sum((pc[ok] < 0.025) | (pc[ok] > 0.975)))
        report["models"][key] = {
            "n_stats": int(okv), "outside_95": outside95, "frac_outside_95": outside95 / okv,
            "n_abs_z_gt_3": int(np.sum(np.abs(zz[ok]) > 3)), "sum_z2": float(np.sum(zz[ok] ** 2)),
            "per_group_rms_z": {g: float(np.sqrt(np.nanmean(zz[[i for i, n in enumerate(names)
                                                               if n.rsplit('_', 1)[0] == g]] ** 2))) for g in groups},
            "z": dict(zip(names, zz.tolist())), "pct": dict(zip(names, pc.tolist()))}
    report["time_s"] = time.time() - t0
    json.dump(report, open(out_js, "w"), indent=1, default=float)
    np.savez(D / f"ppc_{a.tag}.npz", **arrays)

    # figure
    cols = {"modal-OU": "C3", "modal3": "C0", "modal4": "C2"}
    idx = lambda g: [i for i, n in enumerate(names) if n.rsplit("_", 1)[0] == g]
    fig, axs = plt.subplots(2, 4, figsize=(19, 8.5)); axs = axs.ravel()

    def band(ax, x, Y, c, label, log=False):
        lo, md, hi = np.nanpercentile(Y, [5, 50, 95], 0)
        ax.fill_between(x, lo, hi, color=c, alpha=0.25, lw=0); ax.plot(x, md, color=c, lw=1.2, label=label)

    for key in res:
        c = cols[key]
        band(axs[0], T_GRID, arrays[f"{key}_sO"], c, key); band(axs[1], T_GRID, arrays[f"{key}_sC"], c, key)
        V = arrays[f"{key}_v"]
        band(axs[2], sm.LAGS, V[:, idx("oo")], c, key); band(axs[2], sm.LAGS, V[:, idx("cc")], c, None)
        rc = 0.5 * (sm.RATE_EDGES[1:] + sm.RATE_EDGES[:-1])
        band(axs[3], rc, V[:, idx("rate2min")], c, key)
        tck = sm.CK_LAGS * io.DT
        band(axs[4], tck, V[:, idx("ck_open")], c, key); band(axs[5], tck, V[:, idx("ck_closed")], c, key)
        lagsf = np.array([1, 10, 50, 100, 300, 600, 1000, 2000, 3000]) * io.DT
        band(axs[6], lagsf, V[:, idx("acf")], c, key)
        f = V[:, idx("fano")]
        for j in (0, 1):
            ff = f[:, j][np.isfinite(f[:, j])]
            axs[7].boxplot(ff, positions=[j * 4 + list(cols).index(key)], widths=0.7, showfliers=False,
                           patch_artist=True, boxprops=dict(facecolor=c, alpha=0.4))
    ov = obs[0]
    axs[0].plot(T_GRID, obs[1], "k", lw=2, label="data"); axs[1].plot(T_GRID, obs[2], "k", lw=2, label="data")
    axs[2].plot(sm.LAGS, ov[idx("oo")], "ko-", ms=3, label="data O-O"); axs[2].plot(sm.LAGS, ov[idx("cc")], "k^:", ms=3, label="data C-C")
    axs[3].plot(0.5 * (sm.RATE_EDGES[1:] + sm.RATE_EDGES[:-1]), ov[idx("rate2min")], "ko-", ms=3, label="data")
    axs[4].plot(sm.CK_LAGS * io.DT, ov[idx("ck_open")], "ko-", ms=3); axs[5].plot(sm.CK_LAGS * io.DT, ov[idx("ck_closed")], "ko-", ms=3)
    axs[6].plot(np.array([1, 10, 50, 100, 300, 600, 1000, 2000, 3000]) * io.DT, ov[idx("acf")], "ko-", ms=3)
    for j in (0, 1): axs[7].plot([j * 4 - 0.5, j * 4 + 2.5], [ov[idx("fano")][j]] * 2, "k", lw=2)
    for ax in axs[:2]: ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-5, 1.1); ax.set_xlabel("t (s)")
    axs[0].set_title("open dwell survival"); axs[1].set_title("closed dwell survival")
    axs[2].set_xscale("log"); axs[2].set_title("dwell-lag Spearman (O-O, C-C)"); axs[2].set_xlabel("lag (cycles)")
    axs[3].set_title("2-min local-rate histogram"); axs[3].set_xlabel("transitions / s")
    for ax, t in ((axs[4], "CK divergence, q0 = open"), (axs[5], "CK divergence, q0 = closed")):
        ax.set_xscale("log"); ax.set_yscale("symlog", linthresh=1e-4); ax.set_title(t); ax.set_xlabel("t (s)")
    axs[6].set_xscale("log"); axs[6].set_title("ACF of the two-level signal"); axs[6].set_xlabel("lag (s)")
    axs[7].set_xticks([1, 5], ["Fano, slow episodes", "Fano, fast episodes"]); axs[7].set_title("within-episode Fano")
    axs[0].legend(fontsize=8)
    m = report["models"]
    fig.suptitle(f"Phase 1 posterior-predictive checks (bands 5-95 %, n = {a.n} each). Data outside 95 %: "
                 + ", ".join(f"{k} {v['outside_95']}/{v['n_stats']}" for k, v in m.items()), fontsize=10)
    fig.tight_layout(); fig.savefig(Path(a.figdir) / f"phase1_ppc_{a.tag}.png", dpi=110)
    for k, v in m.items():
        print(k, {kk: v[kk] for kk in ("n_stats", "outside_95", "n_abs_z_gt_3", "sum_z2")})
        print("  rms z by group:", {g: round(x, 2) for g, x in v["per_group_rms_z"].items()})
    print(f"total {report['time_s']:.0f} s")


if __name__ == "__main__":
    main()
