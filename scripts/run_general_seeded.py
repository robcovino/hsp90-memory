"""Fit general (fully connected) kO-open/kC-closed models seeded from the modal2 solution, to get the
likelihood ceiling for a given number of states without the random-init local-optimum problem."""
import sys, json, time, numpy as np
from hsp90memory import amm
kO, kC = int(sys.argv[1]), int(sys.argv[2]); key = f"general{kO}o{kC}c_seeded"
dl = []
for k in "1234":
    s, dw, dc = amm.load_dwells(f"data/{k}Thermal.txt"); dl.append((dw, dc))
m = amm.general(kO, kC)
# embed modal2: O_s, C_s, O_f, C_f -> general states O1=O_s, O2=O_f, C1=C_s, C2=C_f ; extras weakly connected
r2 = json.load(open("results/modal2.json")); Q2 = np.array(r2["Q"]); K = kO + kC
Q0 = np.full((K, K), 0.01); map2 = {0: 0, 2: 1, 1: kO, 3: kO + 1}
for i in range(4):
    for j in range(4):
        if i != j: Q0[map2[i], map2[j]] = max(Q2[i, j], 1e-4)
idx = [(i, j) for i in range(K) for j in range(K) if i != j]
rng = np.random.default_rng(1); best = None; t0 = time.time()
for rs in range(3):
    th0 = np.log(np.array([Q0[i, j] for i, j in idx])) + (rng.normal(0, 0.5, len(idx)) if rs else 0)
    from scipy.optimize import minimize
    res = minimize(lambda th: -amm.loglik(m, th, dl), th0, method="L-BFGS-B", options=dict(maxiter=400))
    print(f"restart {rs}: ll={-res.fun:.2f} nit={res.nit}", flush=True)
    if best is None or res.fun < best.fun: best = res
n = sum(len(x[0]) for x in dl); Q = m.Q(best.x)
out = dict(model=key, n_par=m.n_par, theta=list(best.x), rates=list(np.exp(best.x)), ll=float(-best.fun), n_dwells=n,
           bic=float(amm.bic(-best.fun, m.n_par, n)), time=time.time() - t0, classes=list(map(int, m.classes)), Q=Q.tolist(), pi=amm.stationary(Q).tolist())
json.dump(out, open(f"results/{key}.json", "w"), indent=1)
print(key, "ll", out["ll"], "BIC", out["bic"]); print(np.round(Q, 4)); print("pi", np.round(out["pi"], 3))
