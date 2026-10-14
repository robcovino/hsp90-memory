"""Maximum-likelihood fit of a discrete modal-K AMM or a modal-K AMM with grid-OU modulation (amm.modal_ou), exact
likelihood (amm.loglik_Q) on the four pooled traces.

OU models are optimised on a coarse grid (--grid-opt) from the modal-K ML rates and the given (sigma, tau_w) starts,
then refined on --grid-fine from the best coarse optimum, and the log-likelihood is evaluated on --grid-eval for the
grid extrapolation. All log-likelihoods condition on the class of the first frame of each trace (as amm.loglik).
Output: results/exact/fit_<name>.json (settings embedded); log to results/exact/fit_<name>.log.

Usage: python scripts/fit_exact.py modal4_ou [--couple OC|both|eq] [--starts 0.25,170 0.5,30 0.15,1000]
       python scripts/fit_exact.py modal6            (discrete, seeded from modal5 split like fit_modal5.py)
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse, json, re, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from hsp90memory import amm

OUT = Path("results/exact")


def data():
    return [amm.load_dwells(f"data/{k}Thermal.txt")[1:] for k in "1234"]


def negll(model, dl):
    def f(th):
        ll = amm.loglik_Q(model.Q(th), model.classes, dl)
        return -ll if np.isfinite(ll) else 1e12
    return f


def opt(model, dl, x0, maxiter, bounds=None, log=print):
    t0 = time.time()
    r = minimize(negll(model, dl), x0, method="L-BFGS-B", bounds=bounds, options=dict(maxiter=maxiter))
    log(f"  {model.name}: ll={-r.fun:.3f} nit={r.nit} nfev={r.nfev} {time.time() - t0:.0f} s ({r.message})")
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("name")
    ap.add_argument("--couple", default="OC")
    ap.add_argument("--starts", nargs="*", default=["0.25,170", "0.5,30", "0.15,1000"])
    ap.add_argument("--grid-opt", type=int, default=9)
    ap.add_argument("--grid-fine", type=int, default=15)
    ap.add_argument("--grid-eval", type=int, nargs="*", default=[21, 31])
    ap.add_argument("--maxiter", type=int, default=400)
    ap.add_argument("--init-json", default=None, help="JSON list of {'start': name, 'theta': [...]} (full theta, OU "
                    "models: log-rates, log sigma, log tau); replaces the default starts")
    ap.add_argument("--tag", default=None, help="output tag (default from name and coupling)")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    tag = a.tag or a.name + ({"OC": "", "both": "2", "eq": "eq"}[a.couple] if "_ou" in a.name else "")
    out = OUT / f"fit_{tag}.json"
    if out.exists():
        sys.exit(f"{out} exists")
    logf = open(OUT / f"fit_{tag}.log", "a")
    def log(msg):
        print(msg, flush=True); logf.write(msg + "\n"); logf.flush()
    dl = data(); n_dw = sum(len(x[0]) for x in dl); t0 = time.time()
    gm = re.match(r"([A-Za-z]+\d)(_ou)?$", a.name)
    if gm and gm.group(1) in amm.NAMED_GRAPHS:
        topo = gm.group(1); K = amm.NAMED_GRAPHS[topo][0]; is_ou = bool(gm.group(2))
    else:
        mm = re.match(r"(modal|ring)(\d+)(_noDB)?(_ou)?$", a.name); K = int(mm.group(2)); is_ou = bool(mm.group(4))
        topo = "chain" if mm.group(1) == "modal" else ("ring_noDB" if mm.group(3) else "ring")
    res = {"name": tag, "n_dwells": n_dw, "args": vars(a)}
    if is_ou:
        th0 = np.array(json.load(open(f"results/modal{K}.json"))["theta"]) if topo == "chain" else None
        nb = amm.modal_ou(K, n_grid=3, couple=a.couple, topology=topo).base.n_par
        bounds = [(-20, 5)] * nb + ([(np.log(1e-3), np.log(3.0))] +
                                   ([(-3.0, 3.0)] if a.couple == "both" else []) + [(0.0, np.log(1e5))])
        mc = amm.modal_ou(K, n_grid=a.grid_opt, couple=a.couple, topology=topo)
        starts = []
        if a.init_json:
            for d in json.load(open(a.init_json)):
                log(f"start {d['start']} (grid {a.grid_opt})")
                r = opt(mc, dl, np.array(d["theta"]), a.maxiter, bounds, log)
                starts.append({"start": d["start"], "ll": -r.fun, "theta": r.x.tolist()})
        for s in ([] if a.init_json else a.starts):
            sg, tau = map(float, s.split(","))
            x0 = np.r_[th0, np.log(sg), 0.0] if a.couple == "both" else np.r_[th0, np.log(sg)]
            x0 = np.r_[x0, np.log(tau)]
            log(f"start sigma={sg} tau={tau} (grid {a.grid_opt})")
            r = opt(mc, dl, x0, a.maxiter, bounds, log); starts.append({"start": s, "ll": -r.fun, "theta": r.x.tolist()})
        best = max(starts, key=lambda d: d["ll"])
        mf = amm.modal_ou(K, n_grid=a.grid_fine, couple=a.couple, topology=topo)
        log(f"refine on grid {a.grid_fine} from start {best['start']}")
        r = opt(mf, dl, np.array(best["theta"]), a.maxiter, bounds, log)
        th = r.x; ll = -r.fun
        ev = {a.grid_fine: ll}
        for g in a.grid_eval:
            mg = amm.modal_ou(K, n_grid=g, couple=a.couple, topology=topo)
            ev[g] = float(amm.loglik_Q(mg.Q(th), mg.classes, dl)); log(f"  grid {g}: ll={ev[g]:.3f}")
        # Richardson-type extrapolation assuming error ~ 1/n^2 from the two finest grids
        g1, g2 = sorted(ev)[-2:]; ll_inf = ev[g2] + (ev[g2] - ev[g1]) * g1 ** 2 / (g2 ** 2 - g1 ** 2)
        n_par = mf.n_par
        res.update({"model": mf.name, "K_modes": K, "topology": topo, "couple": a.couple, "n_par": n_par, "starts": starts,
                    "grid_fine": a.grid_fine, "theta": th.tolist(), "ll": ll, "ll_grid": {str(k): v for k, v in ev.items()},
                    "ll_extrapolated": float(ll_inf), "bic": float(amm.bic(ll, n_par, n_dw)),
                    "bic_extrapolated": float(amm.bic(ll_inf, n_par, n_dw)),
                    "sigma_OC": float(np.exp(th[nb])), "tau_w": float(np.exp(th[-1])),
                    "sigma_CO": float(th[nb + 1]) if a.couple == "both" else 0.0})
    else:
        m = (amm.named_graph(topo) if topo in amm.NAMED_GRAPHS else
             amm.modal(K) if topo == "chain" else amm.ring(K, db=(topo == "ring"))); Kp = K - 1
        if not a.init_json:
            prev = np.array(json.load(open(f"results/modal{K - 1}.json"))["theta"])
            modes = [prev[2 * i:2 * i + 2] for i in range(Kp)]; coup = [prev[2 * Kp + 3 * i:2 * Kp + 3 * i + 3] for i in range(Kp - 1)]
        starts = []
        if a.init_json:
            for d in json.load(open(a.init_json)):
                log(f"start {d['start']}")
                r = opt(m, dl, np.array(d["theta"]), a.maxiter, None, log); starts.append({"start": d["start"], "ll": -r.fun, "theta": r.x.tolist()})
        for j in ([] if a.init_json else range(Kp)):
            f = np.log(1.5)
            x0 = np.concatenate(modes[:j] + [modes[j] + [f, -f], modes[j] + [-f, f]] + modes[j + 1:]
                                + coup[:j] + [np.full(3, np.log(1 / 300.))] + coup[j:])
            log(f"start split{j}")
            r = opt(m, dl, x0, a.maxiter, None, log); starts.append({"start": f"split{j}", "ll": -r.fun, "theta": r.x.tolist()})
        best = max(starts, key=lambda d: d["ll"]); th = np.array(best["theta"]); ll = best["ll"]
        Q = m.Q(th)
        res.update({"model": m.name, "K_modes": K, "n_par": m.n_par, "starts": starts, "theta": th.tolist(), "ll": ll,
                    "bic": float(amm.bic(ll, m.n_par, n_dw)), "Q": Q.tolist(), "pi": amm.stationary(Q).tolist(),
                    "classes": list(map(int, m.classes))})
    res["time_s"] = time.time() - t0
    json.dump(res, open(out, "w"), indent=1)
    log(f"done {tag}: ll={res['ll']:.3f} BIC={res['bic']:.1f} time {res['time_s']:.0f} s")


if __name__ == "__main__":
    main()
