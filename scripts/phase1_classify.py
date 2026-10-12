"""Phase 1: model posterior over classes by a calibrated classifier on the summary statistics.

Training data: the simulation sets of the classes (equal numbers per class, i.e. equal prior model probabilities),
summaries.summary_vector. Batches listed in --test-batches are held out of every class (the same batches the NPEs
held out), a further 5 % of the training simulations is the validation set (early stopping, temperature scaling).
Preprocessing: drop all-NaN and constant columns; impute partially non-finite columns with the median + indicator;
robust scaling (median, IQR) followed by arcsinh (heavy-tailed CK entries). Classifier: ensemble of MLPs
(3 x 256, ReLU, dropout 0.1), cross-entropy, Adam; softmax temperature fitted on the validation set.
Calibration on the held-out simulations: accuracy, log loss, confusion matrix, reliability per class, ECE.
Applied to the data: p(class | summaries) with the spread over ensemble members.

Caveat (Robert et al. 2011): a model posterior from non-sufficient summaries is the posterior given the summaries,
not given the data; posterior-predictive checks per class (phase1_ppc_classes.py) must accompany it.

Usage: python scripts/phase1_classify.py --name v2 --classes modal3-OU:results/phase1_m3ou_v2 modal4:results/phase1_m4 ...
Outputs: results/phase1_classes/classify_<name>.{json,pt}, figures/phase1_classify_<name>.png
"""
import argparse, json, sys, time, warnings
from pathlib import Path
import numpy as np
import torch
import torch.nn as nn
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from hsp90memory import io, summaries as sm

OUT = Path("results/phase1_classes")


def load_x(d):
    B, X = [], []
    for f in sorted(Path(d).glob("batch_*.npz")):
        with np.load(f) as z:
            B.append(np.full(len(z["x"]), int(z["batch"]))); X.append(z["x"]); names = [str(n) for n in z["names"]]
    return np.concatenate(B), np.concatenate(X), names


class Prep:
    def fit(self, X):
        fin = np.isfinite(X)
        sd = np.array([X[fin[:, j], j].std() if fin[:, j].any() else 0 for j in range(X.shape[1])])
        self.keep = np.flatnonzero(fin.any(0) & (sd > 0))
        self.partial = [j for j in self.keep if not fin[:, j].all()]
        self.med = np.array([np.median(X[fin[:, j], j]) for j in self.keep])
        q75 = np.array([np.percentile(X[fin[:, j], j], 75) for j in self.keep])
        q25 = np.array([np.percentile(X[fin[:, j], j], 25) for j in self.keep])
        self.iqr = np.where(q75 - q25 > 0, q75 - q25, np.array([X[fin[:, j], j].std() for j in self.keep]))
        return self

    def __call__(self, X):
        X = np.atleast_2d(X)
        Z = X[:, self.keep].copy()
        ind = [~np.isfinite(X[:, j]) for j in self.partial]
        bad = ~np.isfinite(Z); Z[bad] = np.broadcast_to(self.med, Z.shape)[bad]
        Z = np.arcsinh((Z - self.med) / self.iqr)
        if ind: Z = np.column_stack([Z] + [b.astype(float) for b in ind])
        return Z.astype(np.float32)


def mlp(d_in, k):
    return nn.Sequential(nn.Linear(d_in, 256), nn.ReLU(), nn.Dropout(0.1), nn.Linear(256, 256), nn.ReLU(),
                         nn.Dropout(0.1), nn.Linear(256, 256), nn.ReLU(), nn.Linear(256, k))


def train_one(Xtr, ytr, Xva, yva, k, seed, max_epochs=200, patience=10):
    torch.manual_seed(seed)
    net = mlp(Xtr.shape[1], k); opt = torch.optim.Adam(net.parameters(), lr=1e-3, weight_decay=1e-5)
    lossf = nn.CrossEntropyLoss(); best, best_state, bad, hist = np.inf, None, 0, []
    Xtr_t, ytr_t, Xva_t, yva_t = map(torch.tensor, (Xtr, ytr, Xva, yva))
    for ep in range(max_epochs):
        net.train(); perm = torch.randperm(len(Xtr_t))
        for i in range(0, len(perm), 512):
            b = perm[i:i + 512]; opt.zero_grad(); l = lossf(net(Xtr_t[b]), ytr_t[b]); l.backward(); opt.step()
        net.eval()
        with torch.no_grad():
            vl = float(lossf(net(Xva_t), yva_t))
        hist.append(vl)
        if vl < best - 1e-4:
            best, bad = vl, 0; best_state = {k_: v.clone() for k_, v in net.state_dict().items()}
        else:
            bad += 1
            if bad >= patience: break
    net.load_state_dict(best_state); net.eval()
    return net, hist


def logits(nets, X):
    with torch.no_grad():
        return torch.stack([n(torch.tensor(X)) for n in nets])          # (n_ens, N, k)


def fit_temperature(lg, y):
    """One temperature for the ensemble-averaged log-probabilities, by grid search on validation log loss."""
    ts = np.exp(np.linspace(np.log(0.3), np.log(5), 200)); yt = torch.tensor(y)
    lp = torch.log_softmax(lg, -1).mean(0)
    nll = [float(nn.functional.cross_entropy(lp / t, yt)) for t in ts]
    return float(ts[int(np.argmin(nll))]), float(min(nll))


def probs(nets, X, T):
    lp = torch.log_softmax(logits(nets, X), -1)                           # members
    return torch.softmax(lp.mean(0) / T, -1).numpy(), torch.softmax(lp / T, -1).numpy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--classes", nargs="+", required=True, help="label:dir")
    ap.add_argument("--test-batches", type=int, nargs="+", default=[0])
    ap.add_argument("--n-ens", type=int, default=5)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--outdir", default=str(OUT))
    ap.add_argument("--figdir", default="figures")
    a = ap.parse_args()
    torch.set_num_threads(a.threads)
    od = Path(a.outdir); od.mkdir(parents=True, exist_ok=True)
    out_js = od / f"classify_{a.name}.json"
    if out_js.exists():
        sys.exit(f"{out_js} exists")
    t0 = time.time()
    labels = [c.split(":")[0] for c in a.classes]; K = len(labels)
    Bs, Xs, ys = [], [], []
    for k, c in enumerate(a.classes):
        B, X, names = load_x(c.split(":")[1]); Bs.append(B); Xs.append(X); ys.append(np.full(len(X), k))
    n_min = min(len(X) for X in Xs)
    for k in range(K):                                                    # equal counts per class: equal prior
        Bs[k], Xs[k], ys[k] = Bs[k][:n_min], Xs[k][:n_min], ys[k][:n_min]
    B, X, y = np.concatenate(Bs), np.concatenate(Xs), np.concatenate(ys)
    test = np.isin(B, a.test_batches)
    rng = np.random.default_rng(0)
    va = (~test) & (rng.random(len(X)) < 0.05); tr = (~test) & (~va)
    prep = Prep().fit(X[tr])
    Xtr, Xva, Xte = prep(X[tr]), prep(X[va]), prep(X[test])
    print(f"classes {labels}; train {tr.sum()}, val {va.sum()}, test {test.sum()}, dim {Xtr.shape[1]}", flush=True)
    nets, hists = [], []
    for s in range(a.n_ens):
        net, h = train_one(Xtr, y[tr], Xva, y[va], K, seed=s); nets.append(net); hists.append(h)
        print(f"member {s}: {len(h)} epochs, best val loss {min(h):.4f}", flush=True)
    T, nll_va = fit_temperature(logits(nets, Xva), y[va])
    p_te, _ = probs(nets, Xte, T); yte = y[test]
    acc = float(np.mean(p_te.argmax(1) == yte)); ll = float(-np.mean(np.log(p_te[np.arange(len(yte)), yte] + 1e-12)))
    conf = np.zeros((K, K), int)
    for t_, p_ in zip(yte, p_te.argmax(1)): conf[t_, p_] += 1
    edges = np.linspace(0, 1, 11); rel = {}; ece = 0.0
    for k in range(K):
        pk = p_te[:, k]; hit = (yte == k)
        rows = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (pk >= lo) & (pk < hi) if hi < 1 else (pk >= lo)
            if m.sum(): rows.append([float(pk[m].mean()), float(hit[m].mean()), int(m.sum())]); ece += m.sum() * abs(pk[m].mean() - hit[m].mean())
        rel[labels[k]] = rows
    ece /= (len(yte) * K)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        xo = sm.summary_vector([io.idealise(io.load_trace(f"data/{k}Thermal.txt")) for k in (1, 2, 3, 4)])[0]
    p_o, p_members = probs(nets, prep(xo), T)
    # out-of-distribution check: nearest-neighbour distance (in the scaled summary space) of the data to training
    # simulations, compared with that of held-out simulations; a fraction near 1 means the data are atypical
    from scipy.spatial import cKDTree
    tree = cKDTree(Xtr[rng.choice(len(Xtr), min(50000, len(Xtr)), replace=False)])
    d_obs = float(tree.query(prep(xo), k=1)[0][0]); d_te = tree.query(Xte[:2000], k=1)[0]
    rep = {"name": a.name, "classes": a.classes, "labels": labels, "n_per_class": int(n_min),
           "n_train": int(tr.sum()), "n_val": int(va.sum()), "n_test": int(test.sum()), "dim": int(Xtr.shape[1]),
           "n_ens": a.n_ens, "epochs": [len(h) for h in hists], "temperature": T, "val_nll": nll_va,
           "test_accuracy": acc, "test_logloss": ll, "test_confusion": conf.tolist(), "test_ece": float(ece),
           "reliability": rel,
           "p_data": dict(zip(labels, p_o[0].tolist())),
           "p_data_members": [dict(zip(labels, m_[0].tolist())) for m_ in p_members],
           "nn_distance_data": d_obs, "frac_heldout_nearer_than_data": float(np.mean(d_te < d_obs)),
           "time_s": time.time() - t0}
    json.dump(rep, open(out_js, "w"), indent=1)
    torch.save({"nets": nets, "prep": prep, "T": T, "labels": labels}, od / f"classify_{a.name}.pt")

    fig, axs = plt.subplots(1, 3, figsize=(16, 4.8))
    for k in range(K):
        r = np.array(rel[labels[k]]); axs[0].plot(r[:, 0], r[:, 1], "o-", ms=3, label=labels[k])
    axs[0].plot([0, 1], [0, 1], "k:"); axs[0].set_xlabel("predicted p(class)"); axs[0].set_ylabel("empirical frequency")
    axs[0].set_title(f"reliability, held-out (ECE {ece:.3f}, acc {acc:.3f})"); axs[0].legend(fontsize=8)
    axs[1].imshow(conf / conf.sum(1, keepdims=True), cmap="Blues", vmin=0, vmax=1)
    for i in range(K):
        for j in range(K): axs[1].text(j, i, f"{conf[i, j] / conf[i].sum():.2f}", ha="center", va="center", fontsize=8)
    axs[1].set_xticks(range(K), labels, rotation=30, fontsize=8); axs[1].set_yticks(range(K), labels, fontsize=8)
    axs[1].set_xlabel("predicted"); axs[1].set_ylabel("true"); axs[1].set_title("confusion (held-out)")
    pm = np.array([m_[0] for m_ in p_members])
    axs[2].bar(range(K), p_o[0], color=[f"C{k}" for k in range(K)])
    for k in range(K): axs[2].plot([k] * len(pm), pm[:, k], "k.", ms=4)
    axs[2].set_xticks(range(K), labels, rotation=30, fontsize=8); axs[2].set_ylim(0, 1)
    axs[2].set_title("p(class | data summaries), equal prior (dots: ensemble members)")
    fig.tight_layout(); fig.savefig(Path(a.figdir) / f"phase1_classify_{a.name}.png", dpi=110)
    print(json.dumps({k: rep[k] for k in ("test_accuracy", "test_logloss", "test_ece", "temperature", "p_data",
                                          "frac_heldout_nearer_than_data", "time_s")}, indent=1))
    print("confusion:\n", conf)


if __name__ == "__main__":
    main()
