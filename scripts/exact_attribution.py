"""Where does one model beat another? Per-dwell predictive log-likelihood terms log p(dwell_i | earlier dwells)
(amm.loglik_terms) at the ML parameters of each model, differences aggregated by trace, by state of the dwell, by
dwell length, by the length of the preceding open dwell, and cumulatively along each trace.

Usage: python scripts/exact_attribution.py --models modal4 modal5 modal6 modal4_ou [--ref modal5] [--grid 15]
Outputs: results/exact/attribution.json, figures/exact_attribution.png
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import argparse, json, re
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from hsp90memory import amm, io

EDGES = np.array([0, 0.5, 1, 2, 5, 10, 20, 50, 100, 1e5])   # s


def model_Q(name, grid):
    m = re.match(r"modal(\d+)(_ou(eq|2)?F?)?$", name); K = int(m.group(1))
    if m.group(2):
        couple = {None: "OC", "eq": "eq", "2": "both"}[m.group(3)]
        mo = amm.modal_ou(K, n_grid=grid, couple=couple)
        th = np.array(json.load(open(f"results/exact/fit_{name}.json"))["theta"]); return mo.Q(th), mo.classes, mo.n_par
    p = next(q for q in (Path(f"results/exact/fit_{name}_global.json"), Path(f"results/{name}.json"),
                         Path(f"results/exact/fit_{name}.json")) if q.exists())
    mo = amm.modal(K); th = np.array(json.load(open(p))["theta"]); return mo.Q(th), mo.classes, mo.n_par


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", nargs="+", default=["modal4", "modal5", "modal6", "modal4_ou"])
    ap.add_argument("--ref", default="modal5")
    ap.add_argument("--grid", type=int, default=15)
    ap.add_argument("--out", default="results/exact/attribution.json")
    a = ap.parse_args()
    if Path(a.out).exists():
        raise SystemExit(f"{a.out} exists")
    dl = [amm.load_dwells(f"data/{k}Thermal.txt")[1:] for k in "1234"]
    terms, npar = {}, {}
    for m in a.models:
        Q, cls, n = model_Q(m, a.grid); terms[m] = amm.loglik_terms(Q, cls, dl); npar[m] = n
        print(m, "ll", sum(t.sum() for t in terms[m]), flush=True)
    out = {"ref": a.ref, "grid": a.grid, "models": {}}
    fig, axs = plt.subplots(1, 4, figsize=(20, 4.5))
    col = {m: f"C{i}" for i, m in enumerate([x for x in a.models if x != a.ref])}
    for m in a.models:
        if m == a.ref: continue
        D = [t - r for t, r in zip(terms[m], terms[a.ref])]
        res = {"n_par_diff": npar[m] - npar[a.ref], "total": float(sum(d.sum() for d in D)),
               "per_trace": [float(d.sum()) for d in D]}
        allD, allL, allC, prevO = [], [], [], []
        for (dw, dc), d in zip(dl, D):
            dur = dw * io.DT
            allD.append(d); allL.append(dur); allC.append(dc)
            po = np.full(len(dw), np.nan)
            for i in range(1, len(dw)):
                if dc[i - 1] == 0: po[i] = dur[i - 1]
                elif i >= 2: po[i] = dur[i - 2]
            prevO.append(po)
        allD, allL, allC, prevO = map(np.concatenate, (allD, allL, allC, prevO))
        res["by_state"] = {"open": float(allD[allC == 0].sum()), "closed": float(allD[allC == 1].sum())}
        by_len = {}
        for c, lab in ((0, "open"), (1, "closed")):
            b = np.digitize(allL, EDGES) - 1
            by_len[lab] = [[float(EDGES[j]), float(EDGES[j + 1]), int(np.sum((b == j) & (allC == c))),
                            float(allD[(b == j) & (allC == c)].sum())] for j in range(len(EDGES) - 1)]
        res["by_length"] = by_len
        b = np.digitize(prevO, EDGES) - 1
        res["by_previous_open"] = [[float(EDGES[j]), float(EDGES[j + 1]), int(np.sum(b == j)),
                                    float(allD[b == j].sum())] for j in range(len(EDGES) - 1)]
        out["models"][m] = res
        ctr = 0.5 * (EDGES[1:] + EDGES[:-1]); ctr[-1] = 200
        for k, lab in enumerate(("open", "closed")):
            v = np.array(by_len[lab]); axs[k].semilogx(np.sqrt(np.maximum(v[:, 0], 0.1) * v[:, 1].clip(max=500)), v[:, 3], "o-", color=col[m], label=m)
        v = np.array(res["by_previous_open"]); axs[2].semilogx(np.sqrt(np.maximum(v[:, 0], 0.1) * v[:, 1].clip(max=500)), v[:, 3], "o-", color=col[m], label=m)
        for i, d in enumerate(D):
            axs[3].plot(np.cumsum(dl[i][0]) * io.DT / 3600, np.cumsum(d), color=col[m], lw=1,
                        label=m if i == 0 else None)
        print(m, "vs", a.ref, json.dumps({k: res[k] for k in ("n_par_diff", "total", "per_trace", "by_state")}), flush=True)
    for ax, t in zip(axs, ("open dwells: sum of dlog p by duration", "closed dwells", "by preceding open dwell",
                           "cumulative dlog L along each trace")):
        ax.set_title(t, fontsize=9); ax.axhline(0, color="k", lw=0.5)
    axs[0].set_xlabel("dwell (s)"); axs[1].set_xlabel("dwell (s)"); axs[2].set_xlabel("preceding open dwell (s)")
    axs[3].set_xlabel("time in trace (h)"); axs[0].legend(fontsize=8); axs[3].legend(fontsize=8)
    fig.suptitle(f"per-dwell log-likelihood differences relative to {a.ref} (ML parameters)", fontsize=10)
    fig.tight_layout(); fig.savefig(Path("figures") / (Path(a.out).stem.replace("attribution", "exact_attribution") + ".png"), dpi=110)
    json.dump(out, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
