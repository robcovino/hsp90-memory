"""Model-comparison table of the exact-likelihood programme, assembled from the result files (no recomputation).

Per model: number of parameters, best maximum log-likelihood (best stored fit; OU models on a 15-point grid), BIC
(n = number of dwells), log evidence (preferred Bayes run; OU models corrected to a 21-point grid), chain spread and
max R-hat of that run, held-out log-likelihood of the latest leave-one-trace-out CV that includes the model, and a flag
for known problems (e.g. multimodal posterior). Output: results/exact/summary_table.{json,md}.

Usage: python scripts/exact_summary.py
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fits import best_fit_file                                    # noqa: E402

E = Path("results/exact")
N_DW = 25571
MODELS = [  # name, description, preferred Bayes tags (first existing wins), flag
    ("modal3", "chain of 3 modes", ["modal3"], ""),
    ("modal4", "chain of 4 modes", ["modal4_v3", "modal4"], ""),
    ("modal5", "chain of 5 modes", ["modal5_v3", "modal5_g"], ""),
    ("modal6", "chain of 6 modes", ["modal6_v3", "modal6_g"], ""),
    ("modal7", "chain of 7 modes", ["modal7_v3"], "multimodal posterior: log Z of the dominant sampled mode"),
    ("modal3_ou", "chain 3 + OU on kOC", [], ""),
    ("modal4_ou", "chain 4 + OU on kOC (non-eq.)", ["modal4_ou"], ""),
    ("modal4_oueq", "chain 4 + OU on kOC (eq.)", ["modal4_oueq_v3", "modal4_oueq"], ""),
    ("modal4_ou2", "chain 4 + OU on kOC and kCO", [], ""),
    ("modal5_oueqF", "chain 5 + OU (eq.)", ["modal5_oueqF"], "log Z from an inferior peak (ML -98034.7 at the time)"),
    ("ring4", "ring of 4 modes", ["ring4"], ""),
    ("ring4_noDB", "ring of 4, cycle not balanced", [], ""),
    ("ringL1", "ring 4 + 1 discrete kOC ladder mode", ["ringL1"], ""),
    ("ringL2", "ring 4 + 2-step discrete kOC ladder", ["ringL2"], ""),
    ("ringS2", "ring 4 + 2 discrete kOC branches", ["ringS2"], ""),
    ("ringF1", "ring 4 + second frozen state", ["ringF1"], ""),
    ("ring4_oueq", "ring 4 + OU on kOC (eq.)", ["ring4_oueq"], ""),
    ("ring4_noDB_oueq", "ring 4 + OU, cycle not balanced", [], ""),
    ("ringF1_oueq", "ring 4 + second frozen state + OU (eq.)", ["ringF1_oueq"], ""),
    ("full5_oueq", "all 10 connections of 5 modes + OU (eq.)", [], "likelihood ceiling for 5 modes"),
]
CVS = ["cv_v6", "cv_v5", "cv_v4", "cv_v3", "cv_v2"]


def main():
    cvs = {c: json.load(open(E / f"{c}.json")) for c in CVS if (E / f"{c}.json").exists()}
    rows = []
    for name, desc, tags, flag in MODELS:
        try:
            f = best_fit_file(name); d = json.load(open(f))
        except (StopIteration, ValueError):
            continue
        ll, k = d["ll"], d["n_par"]
        row = {"model": name, "description": desc, "n_par": k, "max_loglik": ll, "fit_file": str(f),
               "bic": -2 * ll + k * np.log(N_DW), "flag": flag}
        for t in tags:
            p = E / f"bayes_{t}.json"
            if p.exists():
                b = json.load(open(p))
                row.update({"logZ": b["logZ_corrected"].get("21", b["logZ"]), "logZ_run": t,
                            "logZ_chain_spread": float(np.ptp(b["logZ_chain"])), "max_rhat": float(max(b["rhat"]))})
                break
        for c, js in cvs.items():
            if name in js["summary"]:
                row.update({"cv_heldout": js["summary"][name]["ll_test_total"], "cv_run": c}); break
        rows.append(row)
    best_z = max(r.get("logZ", -np.inf) for r in rows)
    for r in rows:
        if "logZ" in r: r["dlogZ_vs_best"] = r["logZ"] - best_z
    json.dump(rows, open(E / "summary_table.json", "w"), indent=1)
    lines = ["| model | description | k | max log L | BIC | log Z | Δlog Z | held-out log L (CV) | note |",
             "|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        z = f"{r['logZ']:.1f}" if "logZ" in r else "—"; dz = f"{r['dlogZ_vs_best']:.1f}" if "logZ" in r else "—"
        cv = f"{r['cv_heldout']:.1f} ({r['cv_run']})" if "cv_heldout" in r else "—"
        lines.append(f"| {r['model']} | {r['description']} | {r['n_par']} | {r['max_loglik']:.1f} | {r['bic']:.1f} | {z} | {dz} | {cv} | {r['flag']} |")
    open(E / "summary_table.md", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
