"""Fit a battery of aggregated Markov models to the four pooled traces (and modal2 per trace).
Usage: python3 run_fits.py <model_key> ; writes results/<model_key>.json"""
import sys, json, time, numpy as np
from hsp90memory import amm

key = sys.argv[1]
dl = []
for k in "1234":
    s, dw, dc = amm.load_dwells(f"data/{k}Thermal.txt"); dl.append((dw, dc))

models = {
    "modal1": amm.modal(1), "modal2": amm.modal(2), "modal2_noDB": amm.modal(2, db=False), "modal3": amm.modal(3),
    "modal4": amm.modal(4),
    "general2o2c": amm.general(2, 2), "general3o2c": amm.general(3, 2), "general3o3c": amm.general(3, 3),
}
t0 = time.time()
if key.startswith("pertrace_"):
    mk = key.split("_", 1)[1]; m = models[mk]; out = {}
    for i, k in enumerate("1234"):
        r = amm.fit(m, [dl[i]], n_restart=2, maxiter=300)
        out[k] = dict(theta=list(r.x), rates=list(np.exp(r.x)), ll=float(-r.fun), n_dwells=int(len(dl[i][0])),
                      bic=float(amm.bic(-r.fun, m.n_par, len(dl[i][0]))))
        print(k, "ll", -r.fun, "rates", np.round(np.exp(r.x), 5), flush=True)
    json.dump(dict(model=mk, n_par=m.n_par, per_trace=out, time=time.time() - t0), open(f"results/{key}.json", "w"), indent=1)
else:
    m = models[key]
    r = amm.fit(m, dl, n_restart=2, maxiter=300, verbose=True)
    n = sum(len(x[0]) for x in dl)
    out = dict(model=key, n_par=m.n_par, theta=list(r.x), rates=list(np.exp(r.x)), ll=float(-r.fun), n_dwells=n,
               bic=float(amm.bic(-r.fun, m.n_par, n)), time=time.time() - t0, classes=list(map(int, m.classes)))
    Q = m.Q(r.x); out["Q"] = Q.tolist(); out["pi"] = amm.stationary(Q).tolist()
    json.dump(out, open(f"results/{key}.json", "w"), indent=1)
    print(key, "ll", -r.fun, "BIC", out["bic"], "time", out["time"], flush=True)
    print("Q (1/s):\n", np.round(Q, 5)); print("pi:", np.round(out["pi"], 4))
