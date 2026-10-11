"""Phase 0, step 2: can the summaries tell within-mode continuous disorder (sigma > 0) from discrete modes (sigma = 0)
at the experimental data volume (4 x 6 h)?  Cross-validated AUC of a classifier (sigma = 0 vs each setting), the
summaries carrying the separation, and the position of the real data relative to each setting.

Usage: python scripts/phase0_analyse.py results/phase0/summaries_A.npz results/phase0/summaries_B.npz
"""
import sys, json, numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import cross_val_score, StratifiedKFold
from hsp90memory import io as hio, summaries as sm

S = {}
for f in sys.argv[1:]:
    z = np.load(f, allow_pickle=True); names = list(z["names"])
    for k in z.files:
        if k != "names": S[k] = z[k]
ref = S["s0"]
data_series = [hio.idealise(hio.load_trace(f"data/{k}Thermal.txt")) for k in "1234"]
vdata, _ = sm.summary_vector(data_series)

def auc_cv(A, B, clf):
    X = np.vstack([A, B]); y = np.r_[np.zeros(len(A)), np.ones(len(B))]
    X = np.nan_to_num(X, nan=0.0)
    cv = StratifiedKFold(5, shuffle=True, random_state=0)
    return cross_val_score(make_pipeline(StandardScaler(), clf), X, y, cv=cv, scoring="roc_auc").mean()

def smd(A, B):
    """standardised mean difference per summary"""
    sd = np.sqrt(0.5 * (np.nanvar(A, 0) + np.nanvar(B, 0))) + 1e-12
    return (np.nanmean(B, 0) - np.nanmean(A, 0)) / sd

rows = []
print(f"{'setting':16s} {'AUC logit':>9s} {'AUC RF':>7s}   top discriminating summaries (standardised mean difference)")
for k in sorted(S, key=lambda s: (s != "amm", s)):
    if k == "s0": continue
    a1 = auc_cv(ref, S[k], LogisticRegression(max_iter=2000, C=0.3))
    a2 = auc_cv(ref, S[k], RandomForestClassifier(200, random_state=0))
    d = smd(ref, S[k]); top = np.argsort(-np.abs(d))[:4]
    print(f"{k:16s} {a1:9.2f} {a2:7.2f}   " + ", ".join(f"{names[i]} {d[i]:+.1f}" for i in top))
    rows.append(dict(setting=k, auc_logit=float(a1), auc_rf=float(a2), top=[(names[i], float(d[i])) for i in top]))

# the fair test: refitted discrete model vs the continuous setting it was fitted to
ref_keys = [k for k in S if k.startswith("refit")]
if ref_keys:
    print("\nRefit test: AUC between replicates of the discrete model refitted to continuous-disorder data and the continuous replicates")
    for k in sorted(ref_keys):
        base = k.split("_", 1)[1]
        if base not in S: continue
        a1 = auc_cv(S[k], S[base], LogisticRegression(max_iter=2000, C=0.3)); a2 = auc_cv(S[k], S[base], RandomForestClassifier(200, random_state=0))
        d = smd(S[k], S[base]); top = np.argsort(-np.abs(d))[:4]
        print(f"  {k:26s} {a1:6.2f} {a2:6.2f}   " + ", ".join(f"{names[i]} {d[i]:+.1f}" for i in top))
        rows.append(dict(setting=k, auc_logit=float(a1), auc_rf=float(a2), top=[(names[i], float(d[i])) for i in top]))

# where do the data sit?  z-score of the data vector against each setting's cloud, summed over the summaries that
# discriminate sigma (|smd| > 1 for the strongest setting), plus a global Mahalanobis-like distance (diagonal).
print("\nData vs settings: diagonal-Mahalanobis distance (all summaries) and the mean z-score on the within-mode summaries")
wm = [i for i, n in enumerate(names) if n.startswith(("fano", "rate1min", "rate2min", "sv2d"))]
for k in sorted(S, key=lambda s: (s != "amm", s)):
    mu = np.nanmean(S[k], 0); sd = np.nanstd(S[k], 0); ok = sd > 1e-6 * (np.abs(mu) + 1e-6)
    z = np.where(ok, (np.nan_to_num(vdata) - mu) / np.where(ok, sd, 1.0), 0.0)
    print(f"  {k:16s} dist {np.sqrt(np.nanmean(z[ok] ** 2)):6.2f}   within-mode z {np.nanmean(z[wm]):+6.2f}  (fano open-slow/fast: data {vdata[names.index('fano_0')]:.2f}/{vdata[names.index('fano_1')]:.2f}, model {mu[names.index('fano_0')]:.2f}/{mu[names.index('fano_1')]:.2f})")
json.dump(dict(rows=rows, data_summary=dict(zip(names, map(float, np.nan_to_num(vdata))))), open("results/phase0/phase0_results.json", "w"), indent=1)
