"""Phase 1: posterior-predictive checks for several model classes, each from its own NPE posterior (symmetric
comparison), optionally with plug-in references at the ML generators.

For each class: n draws from <dir>/npe_<tag>_posterior_data.npz, one simulated measurement (4 x 6 h) each, with the
class's own simulator and settings (<dir>/settings.json). Statistics: summaries.summary_vector and dwell survival
curves on a fixed grid. For the data the same. Per statistic and class: predictive percentile of the data and
z = (data - mean) / sd; counts outside the central 95 % band; sum of z^2.

Class specs: label:<dir>:<npe tag> (NPE posterior, simulated with the class's SBI simulator), or
label:exact:<bayes tag> (exact-likelihood MCMC posterior results/exact/bayes_<tag>_samples.npz; each draw simulated
exactly from its generator with simulators.simulate_amm; grid-OU models on a 15-point grid).

Usage: python scripts/phase1_ppc_classes.py --name all_v2 \\
          --classes modal3-OU:results/phase1_m3ou_v2:m3ou_v2 modal4:results/phase1_m4:m4 modal5x:exact:modal5 ... \\
          [--plugin modal3 modal4] [--n 200] [--workers 4]
Outputs: results/phase1_classes/ppc_<name>.{json,npz}, figures/phase1_ppc_<name>.png
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
sys.path.insert(0, str(Path(__file__).resolve().parent))
from phase1_simulate import simulate_one                                     # noqa: E402
from hsp90memory import io, simulators as sim, summaries as sm               # noqa: E402

T_GRID = np.logspace(-1, 3.3, 44)
OUT = Path("results/phase1_classes")


def survival_grid(series_list, cls):
    d = np.concatenate([(lambda dw, dc: dw[1:-1][dc[1:-1] == cls])(*io.dwells(s)) for s in series_list]) * io.DT
    return np.array([(d > t - 1e-9).mean() for t in T_GRID])


def describe(Y):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        v, _ = sm.summary_vector(Y)
    return v, survival_grid(Y, 0), survival_grid(Y, 1)


def job_class(args):
    S, p, seed = args
    return describe(list(simulate_one(S, p, seed)))


def job_exact(args):
    name, couple, phi, seed = args
    import bayes_exact as be
    T = be.Target(name, 15, couple)
    Q = T.model.Q(T.to_theta(phi))
    return describe(list(sim.simulate_amm(Q, T.model.classes, 216000, n_rep=4, seed=seed)))


def job_plugin(args):
    key, i = args
    r = json.load(open(f"results/{key}.json"))
    return describe(list(sim.simulate_amm(np.array(r["Q"]), r["classes"], 216000, n_rep=4, seed=700000 + 97 * len(key) + i)))


def to_sim_params(samples, names):
    P = samples.copy()
    if "log10_tau_w" in names:
        i = names.index("log10_tau_w"); P[:, i] = 10 ** P[:, i]
    return P


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--classes", nargs="+", required=True, help="label:dir:tag")
    ap.add_argument("--plugin", nargs="*", default=[])
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--outdir", default=str(OUT))
    ap.add_argument("--figdir", default="figures")
    a = ap.parse_args()
    od = Path(a.outdir); od.mkdir(parents=True, exist_ok=True)
    out_js = od / f"ppc_{a.name}.json"
    if out_js.exists():
        sys.exit(f"{out_js} exists")
    t0 = time.time()
    data = [io.idealise(io.load_trace(f"data/{k}Thermal.txt")) for k in (1, 2, 3, 4)]
    obs = describe(data)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        names = sm.summary_vector(data)[1]
    rng = np.random.default_rng(a.seed)
    res, meta = {}, {}
    with get_context("fork").Pool(a.workers) as pool:
        for spec in a.classes:
            label, d, tag = spec.split(":")
            if d == "exact":
                js = json.load(open(f"results/exact/bayes_{tag}.json")); z = np.load(f"results/exact/bayes_{tag}_samples.npz")
                ph = z["phi"].reshape(-1, z["phi"].shape[-1]); ph = ph[rng.choice(len(ph), a.n, replace=False)]
                res[label] = pool.map(job_exact, [(js["args"]["name"], js["args"]["couple"], ph[i], 800000 + i)
                                                  for i in range(a.n)], chunksize=4)
                meta[label] = {"bayes": tag, "kind": "exact-likelihood posterior predictive", "simulator": "simulate_amm"}
                print(f"{label} done {time.time() - t0:.0f} s", flush=True)
                continue
            S = json.load(open(Path(d) / "settings.json"))
            z = np.load(Path(d) / f"npe_{tag}_posterior_data.npz")
            pn = [str(x) for x in z["names"]]
            draws = to_sim_params(z["samples"][rng.choice(len(z["samples"]), a.n, replace=False)], pn)
            res[label] = pool.map(job_class, [(S, draws[i], 900000 + i) for i in range(a.n)], chunksize=4)
            meta[label] = {"dir": d, "tag": tag, "kind": "posterior predictive", "simulator": S["simulator"]}
            print(f"{label} done {time.time() - t0:.0f} s", flush=True)
        for key in a.plugin:
            res[f"{key} (ML)"] = pool.map(job_plugin, [(key, i) for i in range(a.n)], chunksize=4)
            meta[f"{key} (ML)"] = {"kind": "plug-in at ML generator", "simulator": "simulate_amm"}
            print(f"{key} plug-in done {time.time() - t0:.0f} s", flush=True)

    report = {"name": a.name, "n": a.n, "classes": a.classes, "plugin": a.plugin, "meta": meta, "models": {}}
    arrays = {"names": np.array(names), "t_grid": T_GRID, "obs_v": obs[0], "obs_sO": obs[1], "obs_sC": obs[2]}
    groups = sorted({n.rsplit("_", 1)[0] for n in names})
    for key, R in res.items():
        V = np.array([r[0] for r in R]); SO = np.array([r[1] for r in R]); SCc = np.array([r[2] for r in R])
        arrays[f"{key}_v"] = V; arrays[f"{key}_sO"] = SO; arrays[f"{key}_sC"] = SCc
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            ok = np.isfinite(obs[0]) & (np.isfinite(V).mean(0) > 0.9) & (np.nanstd(V, 0) > 0)
            zz = np.full(len(names), np.nan); pc = np.full(len(names), np.nan)
            zz[ok] = (obs[0][ok] - np.nanmean(V[:, ok], 0)) / np.nanstd(V[:, ok], 0)
            pc[ok] = np.array([np.mean(V[np.isfinite(V[:, j]), j] < obs[0][j]) for j in np.flatnonzero(ok)])
        outside = int(np.sum((pc[ok] < 0.025) | (pc[ok] > 0.975)))
        report["models"][key] = {
            "n_stats": int(ok.sum()), "outside_95": outside, "n_abs_z_gt_3": int(np.sum(np.abs(zz[ok]) > 3)),
            "sum_z2": float(np.sum(zz[ok] ** 2)),
            "per_group_rms_z": {g: float(np.sqrt(np.nanmean(zz[[i for i, n in enumerate(names)
                                                               if n.rsplit('_', 1)[0] == g]] ** 2))) for g in groups},
            "z": dict(zip(names, zz.tolist())), "pct": dict(zip(names, pc.tolist()))}
    report["time_s"] = time.time() - t0
    json.dump(report, open(out_js, "w"), indent=1, default=float)
    np.savez(od / f"ppc_{a.name}.npz", **arrays)

    keys = list(res); cols = {k: f"C{i}" for i, k in enumerate(keys)}
    idx = lambda g: [i for i, n in enumerate(names) if n.rsplit("_", 1)[0] == g]
    fig, axs = plt.subplots(2, 4, figsize=(19, 8.5)); axs = axs.ravel()
    lagsf = np.array([1, 10, 50, 100, 300, 600, 1000, 2000, 3000]) * io.DT
    rc = 0.5 * (sm.RATE_EDGES[1:] + sm.RATE_EDGES[:-1]); tck = sm.CK_LAGS * io.DT

    def band(ax, x, Y, c, label):
        lo, md, hi = np.nanpercentile(Y, [5, 50, 95], 0)
        ax.fill_between(x, lo, hi, color=c, alpha=0.18, lw=0); ax.plot(x, md, color=c, lw=1.2, label=label)

    for k in keys:
        c = cols[k]; V = arrays[f"{k}_v"]
        band(axs[0], T_GRID, arrays[f"{k}_sO"], c, k); band(axs[1], T_GRID, arrays[f"{k}_sC"], c, k)
        band(axs[2], sm.LAGS, V[:, idx("oo")], c, k); band(axs[2], sm.LAGS, V[:, idx("cc")], c, None)
        band(axs[3], rc, V[:, idx("rate2min")], c, k)
        band(axs[4], tck, V[:, idx("ck_open")], c, k); band(axs[5], tck, V[:, idx("ck_closed")], c, k)
        band(axs[6], lagsf, V[:, idx("acf")], c, k)
        f = V[:, idx("fano")]
        for j in (0, 1):
            ff = f[:, j][np.isfinite(f[:, j])]
            axs[7].boxplot(ff, positions=[j * (len(keys) + 1) + keys.index(k)], widths=0.7, showfliers=False,
                           patch_artist=True, boxprops=dict(facecolor=c, alpha=0.4))
    ov = obs[0]
    axs[0].plot(T_GRID, obs[1], "k", lw=2, label="data"); axs[1].plot(T_GRID, obs[2], "k", lw=2, label="data")
    axs[2].plot(sm.LAGS, ov[idx("oo")], "ko-", ms=3); axs[2].plot(sm.LAGS, ov[idx("cc")], "k^:", ms=3)
    axs[3].plot(rc, ov[idx("rate2min")], "ko-", ms=3)
    axs[4].plot(tck, ov[idx("ck_open")], "ko-", ms=3); axs[5].plot(tck, ov[idx("ck_closed")], "ko-", ms=3)
    axs[6].plot(lagsf, ov[idx("acf")], "ko-", ms=3)
    for j in (0, 1):
        axs[7].plot([j * (len(keys) + 1) - 0.5, j * (len(keys) + 1) + len(keys) - 0.5], [ov[idx("fano")][j]] * 2, "k", lw=2)
    for ax in axs[:2]: ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-5, 1.1); ax.set_xlabel("t (s)")
    axs[0].set_title("open dwell survival"); axs[1].set_title("closed dwell survival")
    axs[2].set_xscale("log"); axs[2].set_title("dwell-lag Spearman (O-O solid, C-C)"); axs[2].set_xlabel("lag (cycles)")
    axs[3].set_title("2-min local-rate histogram"); axs[3].set_xlabel("transitions / s")
    for ax, t in ((axs[4], "CK divergence, q0 = open"), (axs[5], "CK divergence, q0 = closed")):
        ax.set_xscale("log"); ax.set_yscale("symlog", linthresh=1e-4); ax.set_title(t); ax.set_xlabel("t (s)")
    axs[6].set_xscale("log"); axs[6].set_title("ACF of the two-level signal"); axs[6].set_xlabel("lag (s)")
    nk = len(keys); axs[7].set_xticks([(nk - 1) / 2, nk + 1 + (nk - 1) / 2], ["Fano, slow episodes", "Fano, fast episodes"])
    axs[7].set_title("within-episode Fano"); axs[0].legend(fontsize=7)
    m = report["models"]
    fig.suptitle(f"Posterior-predictive checks, n = {a.n} each (bands 5-95 %). Data outside 95 %: "
                 + ", ".join(f"{k} {v['outside_95']}/{v['n_stats']}" for k, v in m.items()), fontsize=9)
    fig.tight_layout(); fig.savefig(Path(a.figdir) / f"phase1_ppc_{a.name}.png", dpi=110)
    for k, v in m.items():
        print(k, {kk: v[kk] for kk in ("n_stats", "outside_95", "n_abs_z_gt_3", "sum_z2")})
        print("  rms z by group:", {g: round(x, 2) for g, x in v["per_group_rms_z"].items()})
    print(f"total {report['time_s']:.0f} s")


if __name__ == "__main__":
    main()
