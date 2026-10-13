"""Locate the best stored ML fit of a model: among results/<name>.json, results/exact/fit_<name>.json and
results/exact/fit_<name>_global*.json, the one with the highest log-likelihood (all at the same grid for OU models)."""
import json
from pathlib import Path


def best_fit_file(name):
    cands = [Path(f"results/{name}.json"), Path(f"results/exact/fit_{name}.json")] + \
        sorted(Path("results/exact").glob(f"fit_{name}_global*.json"))
    cands = [p for p in cands if p.exists()]
    return max(cands, key=lambda p: json.load(open(p))["ll"])
