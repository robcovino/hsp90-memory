"""Phase 0, step 3 (the fair test): for a continuous-disorder setting, simulate one measurement (4 x 6 h), fit the
discrete modal AMM to it by maximum likelihood, simulate replicates from the refitted discrete model, and store their
summaries.  If a classifier still separates refitted-discrete from continuous replicates, the continuum has a signature
no discrete model of that size absorbs.

Usage: python scripts/phase0_refit.py <settings_json> <out_npz> <n_modes> [n_rep]
"""
import sys, json, time, numpy as np
from hsp90memory import amm, simulators as sim, summaries as sm
from hsp90memory.io import dwells

settings = json.load(open(sys.argv[1])); out = sys.argv[2]; n_modes = int(sys.argv[3])
R = int(sys.argv[4]) if len(sys.argv) > 4 else 50
NF = 216000; NT = 4
r = json.load(open("results/modal3.json")); Q3 = np.array(r["Q"])
kOC, kCO, SO, SC = sim.modal_params_from_Q(Q3)
store = {}; fits = {}
for i, st in enumerate(settings):
    t0 = time.time(); seed = 5000 + 31 * i
    X = sim.simulate_modal_ou(kOC, kCO, SO, SC, st["sigma"], st["tau_w"], st["rho"], NF, n_rep=NT, seed=seed)
    dl = [dwells(X[j]) for j in range(NT)]
    m = amm.modal(n_modes)
    # start from the generating discrete rates (best case for the discrete model)
    r0 = json.load(open(f"results/modal{n_modes}.json")); m.init = np.array(r0["theta"])
    res = amm.fit(m, dl, n_restart=2, maxiter=300)
    Qf = m.Q(res.x)
    Y = sim.simulate_amm(Qf, m.classes, NF, n_rep=R * NT, seed=seed + 1)
    V = []
    for k in range(R):
        v, names = sm.summary_vector([Y[k * NT + j] for j in range(NT)]); V.append(v)
    store[f"refit{n_modes}_" + st["name"]] = np.array(V)
    fits[st["name"]] = dict(theta=list(res.x), ll=float(-res.fun))
    print(f"{st['name']}: refit modal{n_modes} ll={-res.fun:.1f}, {time.time() - t0:.0f} s", flush=True)
np.savez(out, names=np.array(names), **store)
json.dump(fits, open(out.replace(".npz", "_fits.json"), "w"), indent=1)
print("wrote", out)
