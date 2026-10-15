"""Figure for the extended abstract: kinetic model (a, b) and posterior-predictive checks of the best model (c-g).

(a) gating-mode network of ringF1 at z = 0 (ML fit): five modes, each an open (O) and a closed (C) state; mode
    population, mean residence time, mean open/closed dwell; the closing rates (O->C, blue) carry exp(sigma z).
(b) free energy of open and closed along the hidden coordinate z (exact marginals of the full model, ML fit);
    inset: kOC(z)/kOC(0).
(c-g) data vs 200 posterior-predictive data sets (4 x 6 h each) of ringF1 + OU (blue: median and central 95 %) and of
    the discrete ringF1 (orange: median; whiskers where shown: central 95 %). Source: results/audit/ppc_paper.npz
    (scripts/ppc_paper_figure.py).
Output: figures/paper_figure.{png,pdf}
"""
import json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from hsp90memory import amm

BLUE, ORANGE = "#2a78d6", "#eb6834"          # best model (ringF1 + OU), discrete skeleton (ringF1); validated pair
INK, MUTED, FAINT, BOX = "#1b1b1b", "#6b6a66", "#bdbcb7", "#f4f4f2"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.linewidth": 0.6, "axes.edgecolor": INK,
                     "xtick.color": INK, "ytick.color": INK, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "xtick.minor.width": 0.4, "ytick.minor.width": 0.4, "xtick.major.size": 2.5,
                     "ytick.major.size": 2.5, "xtick.minor.size": 1.5, "ytick.minor.size": 1.5,
                     "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
                     "mathtext.fontset": "dejavusans"})

# ------------------------------------------------------------------ model (ML fit)
fit = json.load(open("results/exact/fit_ringF1_oueq_global.json")); th = np.array(fit["theta"])
m = amm.modal_ou(5, n_grid=31, couple="eq", topology="ringF1"); nb = m.base.n_par
sig, tau = np.exp(th[nb]), np.exp(th[nb + 1])
Q0 = m.base.Q(th[:nb]); p0 = amm.stationary(Q0)
pi = amm.stationary(m.Q(th)).reshape(m.n_grid, 10); pmode = pi.sum(0)
z = np.load("results/audit/ppc_paper.npz")


def fmt_t(t):
    if t < 10: return f"{t:.1f} s"
    if t < 60: return f"{t:.0f} s"
    return f"{t / 60:.1f} min" if t < 600 else f"{t / 60:.0f} min"


def band(ax, x, S, color, median=True, alpha=0.22, lw=1.0, ls="-"):
    lo, md, hi = np.nanpercentile(S, [2.5, 50, 97.5], axis=0)
    if alpha: ax.fill_between(x, lo, hi, color=color, alpha=alpha, lw=0)
    if median: ax.plot(x, md, color=color, lw=lw, ls=ls)


def panel(ax, letter, title):
    ax.set_title(f"{letter}", loc="left", fontsize=8.5, fontweight="bold", x=-0.02, pad=3)
    ax.text(0.045, 1.0, title, transform=ax.transAxes, fontsize=7, va="bottom", ha="left", color=INK)


def clean(ax):
    ax.spines[["top", "right"]].set_visible(False)


fig = plt.figure(figsize=(7.2, 6.6))
gs = fig.add_gridspec(3, 6, height_ratios=[1.15, 1, 1], hspace=0.44, wspace=1.25, left=0.07, right=0.985,
                      top=0.965, bottom=0.07)

# ------------------------------------------------------------------ (a) network
ax = fig.add_subplot(gs[0, :4]); ax.set_aspect("equal"); ax.axis("off")
NAMES = {0: "H", 1: "S", 2: "B", 3: "F", 4: "Z"}
LONG = {0: "hub", 1: "slow", 2: "open-biased", 3: "fast", 4: "frozen"}
pos = {4: (0.0, 0.0), 0: (2.5, 0.0), 1: (4.9, 1.05), 3: (4.9, -1.05), 2: (7.3, 0.0)}
bw, bh = 1.25, 0.52
WB = dict(fc="white", ec="none", pad=0.25)
for a, b in [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4)]:
    (xa, ya), (xb, yb) = pos[a], pos[b]
    ax.plot([xa, xb], [ya, yb], color=FAINT, lw=1.6, solid_capstyle="round", zorder=1)
for k in range(5):
    x, y = pos[k]
    ax.add_patch(FancyBboxPatch((x - bw / 2, y - bh / 2), bw, bh, boxstyle="round,pad=0.02,rounding_size=0.1",
                                fc=BOX, ec=MUTED, lw=0.6, zorder=2))
    xo, xc = x - 0.34, x + 0.34
    r = 0.06 + 0.12 * np.sqrt(pmode[2 * k] / pmode.max()); rc = 0.06 + 0.12 * np.sqrt(pmode[2 * k + 1] / pmode.max())
    ax.add_patch(plt.Circle((xo, y), r, fc="white", ec=INK, lw=0.8, zorder=4))
    ax.add_patch(plt.Circle((xc, y), rc, fc=INK, ec=INK, lw=0.8, zorder=4))
    if r > 0.1: ax.text(xo, y, "O", fontsize=5.2, ha="center", va="center", zorder=5, color=INK)
    if rc > 0.1: ax.text(xc, y, "C", fontsize=5.2, ha="center", va="center", zorder=5, color="white")
    kOC, kCO = Q0[2 * k, 2 * k + 1], Q0[2 * k + 1, 2 * k]
    if kOC > 1e-5:
        ax.annotate("", xy=(xc - 0.16, y + 0.065), xytext=(xo + 0.16, y + 0.065),
                    arrowprops=dict(arrowstyle="-|>", lw=1.1, color=BLUE, shrinkA=0, shrinkB=0, mutation_scale=6), zorder=6)
        ax.annotate("", xy=(xo + 0.16, y - 0.065), xytext=(xc - 0.16, y - 0.065),
                    arrowprops=dict(arrowstyle="-|>", lw=0.8, color=INK, shrinkA=0, shrinkB=0, mutation_scale=6), zorder=6)
        sub = f"open {fmt_t(1 / kOC)}, closed {fmt_t(1 / kCO)}"
    else:
        sub = "no direct O⇄C"
    out = sum(p0[i] * Q0[i, j] for i in (2 * k, 2 * k + 1) for j in range(10) if j // 2 != k)
    res = (p0[2 * k] + p0[2 * k + 1]) / out
    head = f"{NAMES[k]} ({LONG[k]}) {100 * (pmode[2*k] + pmode[2*k+1]):.0f}%, {fmt_t(res)}"
    ax.text(x, y + bh / 2 + 0.07, head, ha="center", va="bottom", fontsize=5.8, color=INK, zorder=7, bbox=WB)
    ax.text(x, y - bh / 2 - 0.07, sub, ha="center", va="top", fontsize=5.6, color=MUTED, zorder=7, bbox=WB)
ax.annotate("", xy=(0.36, 1.62), xytext=(0.0, 1.62), arrowprops=dict(arrowstyle="-|>", lw=1.1, color=BLUE,
            mutation_scale=6, shrinkA=0, shrinkB=0))
ax.text(0.5, 1.62, r"closing, $k_\mathrm{OC}^{(m)}e^{\sigma z}$", fontsize=6.3, va="center", ha="left", color=BLUE)
ax.annotate("", xy=(0.0, 1.27), xytext=(0.36, 1.27), arrowprops=dict(arrowstyle="-|>", lw=0.8, color=INK,
            mutation_scale=6, shrinkA=0, shrinkB=0))
ax.text(0.5, 1.27, r"opening, $k_\mathrm{CO}^{(m)}$", fontsize=6.3, va="center", ha="left", color=INK)
ax.set_xlim(-0.75, 8.15); ax.set_ylim(-1.55, 1.85)
panel(ax, "a", "Kinetic model: five gating modes, closing rates modulated by z")

# ------------------------------------------------------------------ (b) hidden coordinate
ax = fig.add_subplot(gs[0, 4:]); clean(ax)
zg = m.z; h = zg[1] - zg[0]
pO = pi[:, 0::2].sum(1); pC = pi[:, 1::2].sum(1)
GO = -np.log(pO / h); GC = -np.log(pC / h); g0 = GO.min(); GO -= g0; GC -= g0
zz = np.linspace(-3, 3, 300); fO = np.polyfit(zg, GO, 2); fC = np.polyfit(zg, GC, 2)
ax.plot(zz, np.polyval(fO, zz), color=INK, lw=1.1)
ax.plot(zz, np.polyval(fC, zz), color=INK, lw=1.1, ls=(0, (3, 1.5)))
ax.text(2.05, np.polyval(fO, 2.05) - 0.25, "open", fontsize=6, ha="left", va="top", color=INK)
ax.text(-1.8, 4.05, "closed", fontsize=6, ha="left", va="center", color=INK)
zO, zC = -fO[1] / (2 * fO[0]), -fC[1] / (2 * fC[0])
for zm, f in ((zO, fO), (zC, fC)):
    ax.plot([zm, zm], [-0.45, np.polyval(f, zm)], color=MUTED, lw=0.5, ls=":")
ax.text((zO + zC) / 2, -0.5, f"minima shifted by σ = {sig:.2f}", fontsize=5.8, ha="center", va="top", color=INK)
ax.set_xlim(-3, 3); ax.set_ylim(-1.2, 5.6)
ax.set_xlabel("hidden coordinate z (sd units)", labelpad=1); ax.set_ylabel(r"$G(z)/k_\mathrm{B}T$", labelpad=1)
ax.text(0.3, 5.35, r"$k_\mathrm{OC}(z)=k_\mathrm{OC}\,e^{\sigma z}$" + f"\n×{np.exp(-sig):.2f} (z = −1) to ×{np.exp(sig):.2f} (z = +1)",
        fontsize=5.7, ha="center", va="top", color=BLUE, linespacing=1.35)
panel(ax, "b", f"Slow coordinate, τ$_w$ = {tau / 60:.1f} min")

# ------------------------------------------------------------------ (c) dwell survival
ax = fig.add_subplot(gs[1, 0:2]); clean(ax); t = z["t_grid"]
for cls, mk, lab in (("sO", "o", "open"), ("sC", "s", "closed")):
    S = np.where(z[f"ringF1_oueq_{cls}"] > 0, z[f"ringF1_oueq_{cls}"], np.nan)
    band(ax, t, S, BLUE)
    band(ax, t, np.where(z[f"ringF1_{cls}"] > 0, z[f"ringF1_{cls}"], np.nan), ORANGE, alpha=0, lw=0.8, ls=(0, (2.5, 1.5)))
    d = z[f"data_{cls}"]; ok = d > 0
    ax.plot(t[ok][::2], d[ok][::2], mk, ms=2.6, mfc="white" if cls == "sO" else INK, mec=INK, mew=0.6, ls="none")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.1, 2000); ax.set_ylim(1e-4, 1.3)
ax.text(110, 0.04, "open", fontsize=6); ax.text(1.6, 3e-3, "closed", fontsize=6)
ax.set_xlabel("dwell time t (s)", labelpad=1); ax.set_ylabel("P(dwell > t)", labelpad=1)
panel(ax, "c", "Dwell-time survival")

# ------------------------------------------------------------------ (d) successive open dwells
ax = fig.add_subplot(gs[1, 2:4]); clean(ax); L = z["oo_lags"]
band(ax, L, z["ringF1_oueq_oo"], BLUE); band(ax, L, z["ringF1_oo"], ORANGE, alpha=0, lw=0.8, ls=(0, (2.5, 1.5)))
ax.plot(L, z["data_oo"], "o", ms=2.6, color=INK, ls="none")
ax.axhline(0, color=FAINT, lw=0.5); ax.set_xscale("log"); ax.set_xlim(0.8, 130); ax.set_ylim(-0.05, 0.42)
ax.set_xlabel("lag n (cycles)", labelpad=1); ax.set_ylabel(r"corr($t_{O,i}$, $t_{O,i+n}$)", labelpad=1)
panel(ax, "d", "Memory between open dwells")

# ------------------------------------------------------------------ (e) CK divergence from closed
ax = fig.add_subplot(gs[1, 4:6]); clean(ax); tl = z["ck_lags_s"]; ok = tl > 1.0
band(ax, tl[ok], z["ringF1_oueq_ck"][:, 1, ok], BLUE)
band(ax, tl[ok], z["ringF1_ck"][:, 1, ok], ORANGE, alpha=0, lw=0.8, ls=(0, (2.5, 1.5)))
ax.plot(tl[ok], z["data_ck"][1, ok], "o", ms=2.6, color=INK, ls="none")
ax.set_xscale("log"); ax.set_xlim(1, 700); ax.set_ylim(bottom=0)
ax.set_xlabel("time t (s)", labelpad=1); ax.set_ylabel(r"$D_{\tau,\,\mathrm{C}}(t)$  ($\tau$ = 1 s)", labelpad=1)
panel(ax, "e", "Non-Markovianity (CK test)")

# ------------------------------------------------------------------ (f) local-rate distribution
ax = fig.add_subplot(gs[2, 0:3]); clean(ax); e = z["rate_edges"]; xc = 0.5 * (e[1:] + e[:-1])
for S, col, off in ((z["ringF1_oueq_rate"], BLUE, -0.022), (z["ringF1_rate"], ORANGE, 0.022)):
    lo, md, hi = np.percentile(S, [2.5, 50, 97.5], axis=0)
    ax.errorbar(xc + off, md, yerr=[md - lo, hi - md], fmt="o", ms=2.4, color=col, elinewidth=0.9, capsize=0, mew=0)
ax.plot(xc, z["data_rate"], "o", ms=2.8, color=INK, ls="none", zorder=5)
ax.set_yscale("log"); ax.set_xlim(0, 1.0); ax.set_ylim(3e-3, 0.7)
ax.set_xlabel(r"transitions per second in a 1-min window (s$^{-1}$)", labelpad=1); ax.set_ylabel("fraction of windows", labelpad=1)
panel(ax, "f", "Local transition rate")

# ------------------------------------------------------------------ (g) closing-rate persistence
ax = fig.add_subplot(gs[2, 3:6]); clean(ax); Lm = z["kacf_lags_min"]
band(ax, Lm, z["ringF1_oueq_kacf"], BLUE)
lo, md, hi = np.percentile(z["ringF1_kacf"], [2.5, 50, 97.5], axis=0)
ax.plot(Lm, md, color=ORANGE, lw=0.8, ls=(0, (2.5, 1.5)))
ax.plot(Lm, z["data_kacf"], "o", ms=2.6, color=INK, ls="none")
ax.axhline(0, color=FAINT, lw=0.5); ax.set_xlim(0, 20.5); ax.set_ylim(-0.15, 0.85)
ax.set_xlabel("lag (min)", labelpad=1); ax.set_ylabel("autocorrelation", labelpad=1)
panel(ax, "g", "Persistence of the local closing rate")
ax.annotate("both models too persistent", xy=(7.9, 0.0), xytext=(9.5, 0.42), fontsize=6, color=INK,
            arrowprops=dict(arrowstyle="-", lw=0.5, color=MUTED))

# ------------------------------------------------------------------ legend for c-g
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
handles = [Line2D([], [], marker="o", color=INK, ls="none", ms=3, label="data (4 molecules × 6 h)"),
           (Patch(fc=BLUE, alpha=0.22, lw=0), Line2D([], [], color=BLUE, lw=1.0)),
           Line2D([], [], color=ORANGE, lw=0.8, ls=(0, (2.5, 1.5)))]
labels = ["data (4 molecules × 6 h)", "modes + slow modulation of k$_\\mathrm{OC}$ (median, 95% predictive)",
          "modes only (median)"]
fig.legend(handles, labels, loc="center", bbox_to_anchor=(0.5, 0.655), ncol=3, frameon=False, fontsize=6.3,
           handlelength=2.0, columnspacing=1.4)
for ext in ("png", "pdf"):
    fig.savefig(f"figures/paper_figure.{ext}", dpi=300)
print("ok", f"sigma={sig:.3f} tau={tau:.0f}")
