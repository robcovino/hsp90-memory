"""Phase 1: neural posterior estimation for one model class (modal-K, discrete or with within-mode OU disorder),
simulation-based calibration, and the posterior for the real traces.

Inputs: <sims>/batch_*.npz (scripts/phase1_simulate.py). Parameters as trained: the log-rates of modalK
(prior N(theta_hat, ln 2)) and, for OU classes, sigma_OC, sigma_CO ~ U[0, 1], log10 tau_w ~ U[0, log10 tau_max],
rho ~ U[-1, 1].
Summaries: summaries.summary_vector without its structurally undefined entries (CK lags <= tau); entries that are
non-finite in some simulations (Fano factors when an episode type is absent) are imputed with the training median and
flagged by an indicator column.
Held-out batches (--test-batches) are never used for training; SBC (rank statistics, 1000 posterior samples per
held-out simulation) and TARP expected coverage are computed on them.

Usage: python scripts/phase1_npe.py --tag v1 [--sims results/phase1 --outdir results/phase1] [--test-batches 0]
Outputs (<outdir>): npe_<tag>.pt (posterior + preprocessing), npe_<tag>.json (settings, training curve, SBC,
coverage, posterior summary for the data), npe_<tag>_posterior_data.npz; figures/phase1_npe_<tag>_*.png.
"""
import argparse, json, os, pickle, sys, time, warnings
from pathlib import Path
import numpy as np
import torch
from torch.distributions import Independent, Normal, Uniform
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sbi.inference import NPE
from sbi.neural_nets import posterior_nn
from sbi.neural_nets.embedding_nets import FCEmbedding
from sbi.utils import MultipleIndependent, within_support
from sbi.diagnostics import check_sbc, check_tarp
from sbi.diagnostics.tarp import _run_tarp, get_tarp_references
from hsp90memory import io, summaries as sm

D = Path("results/phase1")
OU = ("sigma_OC", "sigma_CO", "log10_tau_w", "rho")


def train_names(settings):
    return [("log10_tau_w" if n == "tau_w" else n) for n in settings["param_names"]]


def load_batches(d=D):
    files = sorted(d.glob("batch_*.npz")); B, X, P = [], [], []
    for f in files:
        with np.load(f) as z:
            B.append(np.full(len(z["x"]), int(z["batch"]))); X.append(z["x"]); P.append(z["params"])
            names = [str(n) for n in z["names"]]; settings = json.loads(str(z["settings"]))
    return np.concatenate(B), np.concatenate(X), np.concatenate(P), names, settings


def to_train_params(P, settings):
    T = P.copy()
    if "tau_w" in settings["param_names"]:
        i = settings["param_names"].index("tau_w"); T[:, i] = np.log10(P[:, i])
    return T


def make_prior(settings):
    th = torch.tensor(settings["theta_hat"], dtype=torch.float32)
    sd = torch.full_like(th, settings["theta_sd"])
    if settings.get("discrete"):
        return Independent(Normal(th, sd), 1)
    return MultipleIndependent([
        Independent(Normal(th, sd), 1),
        Independent(Uniform(torch.zeros(1), torch.ones(1)), 1),
        Independent(Uniform(torch.zeros(1), torch.ones(1)), 1),
        Independent(Uniform(torch.zeros(1), torch.tensor([np.log10(settings["tau_w_loguniform"][1])],
                                                           dtype=torch.float32)), 1),
        Independent(Uniform(-torch.ones(1), torch.ones(1)), 1)])


class Preproc:
    """Drop all-NaN and constant columns (e.g. sv2d_0 = 1 by normalisation); impute partially non-finite columns with the training median and add indicators."""
    def fit(self, X, names):
        fin = np.isfinite(X)
        sd = np.array([X[fin[:, j], j].std() if fin[:, j].any() else 0.0 for j in range(X.shape[1])])
        self.keep = np.flatnonzero(fin.any(0) & (sd > 0))
        self.partial = [j for j in self.keep if not fin[:, j].all()]
        self.med = {j: float(np.median(X[fin[:, j], j])) for j in self.partial}
        self.names = [names[j] for j in self.keep] + [f"isnan_{names[j]}" for j in self.partial]
        return self

    def __call__(self, X):
        X = np.atleast_2d(X)
        out = X[:, self.keep].copy()
        ind = []
        for j in self.partial:
            bad = ~np.isfinite(X[:, j]); out[bad, list(self.keep).index(j)] = self.med[j]; ind.append(bad)
        if ind: out = np.column_stack([out] + [b.astype(float) for b in ind])
        return out


def observed_summary():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return sm.summary_vector([io.idealise(io.load_trace(f"data/{k}Thermal.txt")) for k in (1, 2, 3, 4)])[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--test-batches", type=int, nargs="+", default=[0])
    ap.add_argument("--n-sbc", type=int, default=1000)
    ap.add_argument("--n-post", type=int, default=1000, help="posterior samples per SBC simulation")
    ap.add_argument("--max-epochs", type=int, default=500)
    ap.add_argument("--hidden", type=int, default=128)
    ap.add_argument("--transforms", type=int, default=5)
    ap.add_argument("--embed", type=int, default=64)
    ap.add_argument("--batch", type=int, default=256)
    ap.add_argument("--lr", type=float, default=5e-4)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--sims", default=str(D), help="directory with batch_*.npz")
    ap.add_argument("--outdir", default=str(D))
    ap.add_argument("--figdir", default="figures")
    a = ap.parse_args()
    od, fd = Path(a.outdir), Path(a.figdir); od.mkdir(parents=True, exist_ok=True); fd.mkdir(exist_ok=True)
    torch.set_num_threads(a.threads); torch.manual_seed(a.seed); np.random.seed(a.seed)
    out_pt, out_js = od / f"npe_{a.tag}.pt", od / f"npe_{a.tag}.json"
    if out_js.exists():
        sys.exit(f"{out_js} exists; choose a new --tag (results are never overwritten)")
    t_start = time.time()

    B, X, P, names, settings = load_batches(Path(a.sims))
    test = np.isin(B, a.test_batches); train = ~test
    T = to_train_params(P, settings)
    TN = train_names(settings); NP = len(TN); DIS = [i for i, n in enumerate(TN) if n in OU]
    pp = Preproc().fit(X[train], names)
    xtr = torch.tensor(pp(X[train]), dtype=torch.float32); ttr = torch.tensor(T[train], dtype=torch.float32)
    xte = torch.tensor(pp(X[test]), dtype=torch.float32)[: a.n_sbc]; tte = torch.tensor(T[test], dtype=torch.float32)[: a.n_sbc]
    print(f"train {len(xtr)} sims, test {len(xte)}, summary dim {xtr.shape[1]} "
          f"(dropped {[n for n in names if n not in pp.names]}, imputed {len(pp.partial)})", flush=True)

    prior = make_prior(settings)
    if out_pt.exists():
        # a run whose evaluation did not finish: reuse the trained estimator (same data, split and preprocessing)
        saved = torch.load(out_pt, weights_only=False)
        posterior = saved["posterior"]; s = saved.get("summary", {}); t_train = saved.get("train_time_s", np.nan)
        assert pickle.loads(saved["preproc"]).names == pp.names
        print(f"loaded trained posterior from {out_pt}", flush=True)
    else:
        emb = FCEmbedding(input_dim=xtr.shape[1], output_dim=a.embed, num_layers=3, num_hiddens=256)
        de = posterior_nn("nsf", hidden_features=a.hidden, num_transforms=a.transforms, embedding_net=emb)
        inf = NPE(prior=prior, density_estimator=de, show_progress_bars=False)
        inf.append_simulations(ttr, xtr)
        t0 = time.time()
        est = inf.train(training_batch_size=a.batch, learning_rate=a.lr, stop_after_epochs=20,
                        max_num_epochs=a.max_epochs)
        t_train = time.time() - t0
        s = {k: [float(v) for v in vv] for k, vv in inf._summary.items() if isinstance(vv, list)}
        print(f"trained in {t_train:.0f} s, {s['epochs_trained'][-1]:.0f} epochs, best validation loss (-log q) "
              f"{s['best_validation_loss'][-1]:.3f}", flush=True)
        posterior = inf.build_posterior(est)
        torch.save({"posterior": posterior, "preproc": pickle.dumps(pp), "train_names": TN,
                    "summary": s, "train_time_s": t_train}, out_pt)

    # held-out log probability of the true parameters (flow density, without the leakage normalisation)
    with torch.no_grad():
        lp_test = np.concatenate([posterior.log_prob_batched(tte[i:i + 100].unsqueeze(0), x=xte[i:i + 100],
                                                             norm_posterior=False).squeeze(0).numpy()
                                  for i in range(0, len(xte), 100)])

    # SBC and TARP on held-out simulations; posterior samples drawn in chunks to bound memory
    t0 = time.time()
    S = torch.cat([posterior.sample_batched((a.n_post,), x=xte[i:i + 50], show_progress_bars=False)
                   for i in range(0, len(xte), 50)], 1)                        # (n_post, n_sbc, n_par)
    ranks = (S < tte.unsqueeze(0)).sum(0).float()                               # marginal ranks
    stats = check_sbc(ranks, tte, S[0], num_posterior_samples=a.n_post)
    ecp, alpha = _run_tarp(S, tte, get_tarp_references(tte), z_score_theta=True)
    atc, ks_p = check_tarp(ecp, alpha)
    t_sbc = time.time() - t0
    # per-parameter marginal coverage of central credible intervals, and posterior contraction
    levels = np.array([0.5, 0.8, 0.9, 0.95])
    r = ranks.numpy() / a.n_post
    cov = {TN[i]: [float(np.mean(np.abs(r[:, i] - 0.5) <= l / 2)) for l in levels] for i in range(NP)}
    Sn = S.numpy()
    post_sd = Sn.std(0); post_err = np.abs(np.median(Sn, 0) - tte.numpy())
    prior_sd = T[train].std(0)
    contraction = 1 - (np.mean(post_sd, 0) / prior_sd) ** 2
    np.savez(od / f"npe_{a.tag}_sbc.npz", ranks=ranks.numpy(), theta=tte.numpy(), post_median=np.median(Sn, 0),
             post_sd=post_sd, names=np.array(TN))
    del S, Sn

    # posterior for the data
    xo = torch.tensor(pp(observed_summary()), dtype=torch.float32)
    so = posterior.sample((20000,), x=xo, show_progress_bars=False).numpy()
    with torch.no_grad():
        raw = posterior.posterior_estimator.sample((20000,), condition=xo.reshape(1, -1)).reshape(-1, NP)
    leak = float(1 - within_support(prior, raw).float().mean())
    q = np.quantile(so, [0.05, 0.25, 0.5, 0.75, 0.95], 0)
    post_data = {n: {"q05": q[0, i], "q25": q[1, i], "median": q[2, i], "q75": q[3, i], "q95": q[4, i],
                     "mean": float(so[:, i].mean()), "sd": float(so[:, i].std())} for i, n in enumerate(TN)}
    post_data["leakage_outside_prior"] = leak
    if DIS:
        i1, i2 = TN.index("sigma_OC"), TN.index("sigma_CO")
        post_data["P(sigma_OC<0.1)"] = float(np.mean(so[:, i1] < 0.1))
        post_data["P(sigma_CO<0.1)"] = float(np.mean(so[:, i2] < 0.1))
        post_data["P(max sigma<0.1)"] = float(np.mean(np.maximum(so[:, i1], so[:, i2]) < 0.1))
        post_data["corr_disorder"] = np.corrcoef(so[:, DIS].T).tolist()
    np.savez(od / f"npe_{a.tag}_posterior_data.npz", samples=so, names=np.array(TN), x_obs=xo.numpy())

    rep = {"tag": a.tag, "args": vars(a), "simulation_settings": settings, "n_train": int(train.sum()),
           "n_test": int(len(xte)), "test_batches": a.test_batches, "summary_dim": int(xtr.shape[1]),
           "summary_names": pp.names, "train_time_s": t_train, "sbc_time_s": t_sbc,
           "total_time_s": time.time() - t_start,
           "training": s,
           "heldout_logprob_mean": float(np.mean(lp_test)),
           "sbc": {"ks_pvals": dict(zip(TN, stats["ks_pvals"].numpy().tolist())),
                   "c2st_ranks": dict(zip(TN, stats["c2st_ranks"].numpy().ravel().tolist())),
                   "c2st_dap": dict(zip(TN, stats["c2st_dap"].numpy().ravel().tolist()))},
           "tarp": {"atc": float(atc), "ks_pval": float(ks_p)},
           "coverage_levels": levels.tolist(), "marginal_coverage": cov,
           "contraction": dict(zip(TN, contraction.tolist())),
           "mean_abs_error_heldout": dict(zip(TN, np.mean(post_err, 0).tolist())),
           "heldout_post_sd_mean": dict(zip(TN, np.mean(post_sd, 0).tolist())),
           "prior_sd": dict(zip(TN, prior_sd.tolist())),
           "posterior_data": post_data}
    json.dump(rep, open(out_js, "w"), indent=1, default=float)

    # figures: SBC rank histograms, coverage, posterior for the data
    nc = 5; nr = int(np.ceil(NP / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(16, 2.5 * nr)); nb = 20
    for ax in axs.ravel()[NP:]: ax.axis("off")
    lo, hi = np.percentile(np.random.binomial(len(r), 1 / nb, (5000, nb)), [0.5, 99.5])
    for i, ax in enumerate(axs.ravel()[:NP]):
        ax.hist(ranks[:, i].numpy(), bins=nb, range=(0, a.n_post), color="C0")
        ax.axhspan(lo, hi, color="0.85", zorder=0)
        ax.set_title(f"{TN[i]}  KS p={stats['ks_pvals'][i]:.2f}", fontsize=8); ax.set_yticks([])
    fig.suptitle(f"SBC ranks, {len(r)} held-out simulations (grey: 99 % band for uniform)"); fig.tight_layout()
    fig.savefig(fd / f"phase1_npe_{a.tag}_sbc.png", dpi=110)

    fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))
    axs[0].plot(alpha.numpy(), ecp.numpy(), "C0"); axs[0].plot([0, 1], [0, 1], "k:")
    axs[0].set_xlabel("credibility level"); axs[0].set_ylabel("expected coverage")
    axs[0].set_title(f"TARP (joint): ATC={float(atc):.3f}, KS p={float(ks_p):.2f}")
    for i in (DIS or range(NP)):
        axs[1].plot(levels, cov[TN[i]], "o-", label=TN[i])
    axs[1].plot([0.4, 1], [0.4, 1], "k:"); axs[1].legend(fontsize=8); axs[1].set_xlabel("central interval level")
    axs[1].set_ylabel("empirical coverage"); axs[1].set_title("marginal coverage" + (", disorder parameters" if DIS else ""))
    fig.tight_layout(); fig.savefig(fd / f"phase1_npe_{a.tag}_coverage.png", dpi=110)

    fig, axs = plt.subplots(nr, nc, figsize=(16, 2.5 * nr)); pr = prior.sample((20000,)).numpy()
    for ax in axs.ravel()[NP:]: ax.axis("off")
    for i, ax in enumerate(axs.ravel()[:NP]):
        ax.hist(pr[:, i], 60, density=True, color="0.8", label="prior")
        ax.hist(so[:, i], 60, density=True, color="C3", alpha=0.7, label="posterior (data)")
        ax.set_title(TN[i], fontsize=8); ax.set_yticks([])
    axs.ravel()[0].legend(fontsize=7); fig.tight_layout(); fig.savefig(fd / f"phase1_npe_{a.tag}_posterior.png", dpi=110)

    print(json.dumps({k: rep[k] for k in ("n_train", "summary_dim", "train_time_s", "sbc_time_s", "tarp",
                                          "heldout_logprob_mean")}, default=float))
    print("SBC KS p:", {k: round(v, 3) for k, v in rep["sbc"]["ks_pvals"].items()})
    print("contraction:", {k: round(v, 3) for k, v in rep["contraction"].items()})
    print("posterior (data):", {n: [round(post_data[n][k], 3) for k in ("q05", "median", "q95")] for n in TN})
    print({k: v for k, v in post_data.items() if k.startswith("P(")})


if __name__ == "__main__":
    main()
