"""Thinned copies of the MCMC posterior samples for version control (the full results/exact/bayes_*_samples.npz are
not committed). Keeps ~--n draws per run spread evenly over chains and iterations: phi (sampling parameters: log-rates,
sigma, log tau), log-likelihood, chain index. Output: results/exact/thin/bayes_<tag>_thin.npz.

Usage: python scripts/thin_samples.py [--n 2000]
"""
import argparse, json
from pathlib import Path
import numpy as np

ap = argparse.ArgumentParser(); ap.add_argument("--n", type=int, default=2000); a = ap.parse_args()
out = Path("results/exact/thin"); out.mkdir(exist_ok=True)
for f in sorted(Path("results/exact").glob("bayes_*_samples.npz")):
    tag = f.name[len("bayes_"):-len("_samples.npz")]; dst = out / f"bayes_{tag}_thin.npz"
    if dst.exists(): continue
    z = np.load(f); C, N, d = z["phi"].shape; k = max(1, (C * N) // a.n)
    idx = np.arange(0, N, k)
    js = json.load(open(f"results/exact/bayes_{tag}.json"))
    np.savez_compressed(dst, phi=z["phi"][:, idx].reshape(-1, d), loglik=z["loglik"][:, idx].ravel(),
                        chain=np.repeat(np.arange(C), len(idx)), thin=k, model=js["model"], args=json.dumps(js["args"]))
    print(tag, dst.stat().st_size // 1024, "kB")
