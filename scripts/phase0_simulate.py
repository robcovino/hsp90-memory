"""Phase 0, step 1: simulate replicate measurements (4 x 6 h at 10 Hz) from the nested modal-OU model over a grid of
within-mode disorder parameters and store their summary vectors.

Usage: python scripts/phase0_simulate.py <settings_json> <out_npz> [n_rep]
settings_json: list of {"name":..., "sigma":..., "tau_w":..., "rho":...} or {"name":"amm"} for the exact AMM reference.
"""
import sys, json, time, numpy as np
from hsp90memory import simulators as sim, summaries as sm

settings = json.load(open(sys.argv[1])); out = sys.argv[2]; R = int(sys.argv[3]) if len(sys.argv) > 3 else 50
NF = 216000; NT = 4
r = json.load(open("results/modal3.json")); Q = np.array(r["Q"]); classes = r["classes"]
kOC, kCO, SO, SC = sim.modal_params_from_Q(Q)
store = {}
for i, st in enumerate(settings):
    t0 = time.time(); seed = 1000 + 17 * i
    if st["name"] == "amm":
        X = sim.simulate_amm(Q, classes, NF, n_rep=R * NT, seed=seed)
    else:
        X = sim.simulate_modal_ou(kOC, kCO, SO, SC, st["sigma"], st["tau_w"], st["rho"], NF, n_rep=R * NT, seed=seed)
    V = []
    for k in range(R):
        v, names = sm.summary_vector([X[k * NT + j] for j in range(NT)]); V.append(v)
    store[st["name"]] = np.array(V)
    print(f"{st['name']}: {R} replicates, {time.time() - t0:.0f} s", flush=True)
np.savez(out, names=np.array(names), **store)
print("wrote", out)
