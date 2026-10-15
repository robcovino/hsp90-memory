"""Audit: continuous Gaussian modulator vs a two-level (possibly asymmetric) switch of the closing rates.
ringF1 + a hidden two-level variable z in {z0, z1} with populations (1-p, p), standardised to mean 0, variance 1
(z1 = sqrt((1-p)/p), z0 = -sqrt(p/(1-p))), relaxation time tau_w; every closing rate (O->C) is multiplied by
exp(sigma z). Equilibrium coupling as amm.modal_ou(couple='eq'): z relaxes in U_O = -ln pi_z while open and in
U_O - sigma z while closed (symmetric part of the generator kept), so detailed balance holds. p = 1/2 is the
symmetric switch (= modal_ou n_grid=2); p != 1/2 gives rare episodes of faster (p < 1/2) or slower (p > 1/2) closing.
Parameters: 24 ringF1 log-rates, log sigma, log tau_w, logit p (27).
Output: results/audit/switch2.json, log results/audit/switch2.log.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
from hsp90memory import amm

OUT = Path("results/audit")
base = amm.named_graph("ringF1"); K2 = base.K; nb = base.n_par
O = np.arange(0, K2, 2); C = O + 1


def Q_switch(theta):
    Q0 = base.Q(theta[:nb]).copy(); np.fill_diagonal(Q0, 0)
    sig, tau, p = np.exp(theta[nb]), np.exp(theta[nb + 1]), expit(theta[nb + 2])
    z = np.array([-np.sqrt(p / (1 - p)), np.sqrt((1 - p) / p)])
    R1 = np.array([[-p, p], [1 - p, -(1 - p)]]) / tau                 # relaxation rate 1/tau, pi = (1-p, p)
    S = np.sqrt(R1 * R1.T); UO = -np.log([1 - p, p]); UC = UO - sig * z
    RC = S * np.exp(-(UC[None, :] - UC[:, None]) / 2); np.fill_diagonal(RC, 0); np.fill_diagonal(RC, -RC.sum(1))
    PO = np.zeros((K2, K2)); PO[O, O] = 1; PC = np.zeros((K2, K2)); PC[C, C] = 1
    Q = np.kron(R1, PO) + np.kron(RC, PC)
    for i in range(2):
        Qi = Q0.copy(); Qi[O, C] *= np.exp(sig * z[i]); Q[i * K2:(i + 1) * K2, i * K2:(i + 1) * K2] += Qi
    np.fill_diagonal(Q, 0); np.fill_diagonal(Q, -Q.sum(1)); return Q


def main():
    OUT.mkdir(parents=True, exist_ok=True); logf = open(OUT / "switch2.log", "a")
    def log(s): print(s, flush=True); logf.write(s + "\n"); logf.flush()
    dl = [amm.load_dwells(f"data/{k}Thermal.txt")[1:] for k in "1234"]
    classes = list(base.classes) * 2
    # check: p = 1/2 reproduces modal_ou(n_grid=2, eq)
    th15 = np.array(json.load(open("results/exact/fit_ringF1_oueq_global.json"))["theta"])
    m2 = amm.modal_ou(5, n_grid=2, couple="eq", topology="ringF1")
    a = amm.loglik_Q(m2.Q(th15), m2.classes, dl); b = amm.loglik_Q(Q_switch(np.r_[th15, 0.0]), classes, dl)
    log(f"check p=1/2 vs modal_ou g2 at the g15 optimum: {a:.4f} {b:.4f}")
    f = lambda th: (lambda v: -v if np.isfinite(v) else 1e12)(amm.loglik_Q(Q_switch(th), classes, dl))
    bounds = [(-20, 5)] * nb + [(np.log(1e-3), np.log(5.0)), (0.0, np.log(1e5)), (-6.0, 6.0)]
    g2 = json.load(open(OUT / "grid_switch.json"))["grids"]["2"]["best"]["theta"] if (OUT / "grid_switch.json").exists() else th15
    g2 = np.array(g2); out = []
    for p in (0.5, 0.2, 0.8, 0.07, 0.93):
        # keep the rate contrast roughly fixed: sigma*(z1 - z0) = sigma / sqrt(p(1-p)) equal to the symmetric fit's 2 sigma
        s0 = np.exp(g2[nb]) * 2 * np.sqrt(p * (1 - p))
        x0 = np.r_[g2[:nb], np.log(s0), g2[nb + 1], np.log(p / (1 - p))]
        t0 = time.time(); r = minimize(f, x0, method="L-BFGS-B", bounds=bounds, options=dict(maxiter=400))
        sig, tau, pp = np.exp(r.x[nb]), np.exp(r.x[nb + 1]), expit(r.x[nb + 2])
        z = np.array([-np.sqrt(pp / (1 - pp)), np.sqrt((1 - pp) / pp)])
        out.append(dict(start_p=p, ll=float(-r.fun), sigma=float(sig), tau_w=float(tau), p=float(pp),
                        kOC_ratio=float(np.exp(sig * (z[1] - z[0]))), nit=int(r.nit), theta=r.x.tolist(), time_s=time.time() - t0))
        log(f"start p={p}: ll={-r.fun:.3f} sigma={sig:.3f} tau={tau:.1f} p={pp:.3f} kOC ratio high/low={out[-1]['kOC_ratio']:.2f} "
            f"nit={r.nit} {time.time()-t0:.0f}s")
        json.dump(dict(model="ringF1 + two-level switch (eq coupling)", n_par=nb + 3, starts=out,
                       best=max(out, key=lambda d: d["ll"]), ref={"ringF1": -97996.319, "ringF1_oueq_g15": -97918.856}),
                  open(OUT / "switch2.json", "w"), indent=1)


if __name__ == "__main__":
    main()
