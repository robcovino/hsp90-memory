"""Sensitivity of the summary-based model odds (phase1_classify v2) to the prior of the modal4-OU class.
Restricting a class prior to a region A multiplies its evidence by P(A | x) / P(A); P(A | x) from the modal4-OU NPE
posterior for the data, P(A) from the (uniform / log-uniform) prior. Output: results/phase1_classes/prior_sensitivity_v2.json
"""
import json
import numpy as np

c = json.load(open("results/phase1_classes/classify_v2.json")); p = c["p_data"]
z = np.load("results/phase1_m4ou/npe_m4ou_posterior_data.npz"); S = z["samples"]; n = list(z["names"])
sOC, sCO, lt = (S[:, n.index(k)] for k in ("sigma_OC", "sigma_CO", "log10_tau_w"))
lmax = np.log10(3000.0)
regions = {"full prior": (1, 1, 1, 3000), "sigma_OC<=0.5": (0.5, 1, 1, 3000),
           "sigma_OC<=0.5, sigma_CO<=0.2": (0.5, 0.2, 1, 3000),
           "sigma_OC<=0.5, sigma_CO<=0.2, tau_w 30-3000 s": (0.5, 0.2, 30, 3000)}
out = {"classifier_p": p, "regions": {}}
disc = p["modal4"] + p["modal5"]
for name, (a1, a2, tl, th) in regions.items():
    PA = a1 * a2 * (np.log10(th) - np.log10(tl)) / lmax
    PAx = float(np.mean((sOC <= a1) & (sCO <= a2) & (lt >= np.log10(tl)) & (lt <= np.log10(th))))
    out["regions"][name] = {"P(A)": PA, "P(A|x)": PAx, "factor": PAx / PA, "odds_discrete_vs_OU": disc / (p["modal4-OU"] * PAx / PA)}
    print(f"{name:48s} P(A)={PA:.3f} P(A|x)={PAx:.3f} odds discrete:OU = {out['regions'][name]['odds_discrete_vs_OU']:.2f}")
json.dump(out, open("results/phase1_classes/prior_sensitivity_v2.json", "w"), indent=1)
