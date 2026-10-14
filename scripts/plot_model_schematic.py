"""Schematic of the best-fitting model: gating-mode network (ring + second frozen state, 'ringF1') with equilibrium
OU modulation of the open->closed rate. Numbers from the ML fit (results/exact/fit_<tag>_global.json) and the
posterior (results/exact/bayes_<tag>.json). Output: figures/model_schematic.{png,svg}.

Usage: python scripts/plot_model_schematic.py [--tag ringF1_oueq]
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fits import best_fit_file                                                   # noqa: E402
from hsp90memory import amm                                                       # noqa: E402

INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#8a8984", "#d9d8d4", "#fcfcfb"
MODE = {"slow": "#2a78d6", "open-biased": "#eb6834", "fast": "#1baf7a"}          # categorical slots 1-3
FROZEN, HUB = "#6f6e6a", "#a9a8a3"


def fmt_t(s):
    return f"{s:.1f} s" if s < 10 else (f"{s:.0f} s" if s < 120 else f"{s / 60:.1f} min")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", default="ringF1_oueq"); a = ap.parse_args()
    fit = json.load(open(best_fit_file(a.tag))); th = np.array(fit["theta"])
    bay = json.load(open(f"results/exact/bayes_{a.tag}.json"))
    g = amm.named_graph("ringF1"); Q = g.Q(th[:g.n_par]); pi = amm.stationary(Q)
    sig = (bay["phi_q05"][-2], bay["phi_median"][-2], bay["phi_q95"][-2])
    tau = tuple(np.exp([bay["phi_q05"][-1], bay["phi_median"][-1], bay["phi_q95"][-1]]))
    logZ = bay["logZ_corrected"].get("21", bay["logZ"])
    idx = {"hub": 0, "slow": 1, "open-biased": 2, "fast": 3, "X": 4}
    share = {k: pi[2 * i] + pi[2 * i + 1] for k, i in idx.items()}
    dwell = {k: (1 / Q[2 * i, 2 * i + 1], 1 / Q[2 * i + 1, 2 * i]) for k, i in idx.items()}

    fig = plt.figure(figsize=(15, 7.6), facecolor=SURF)
    ax = fig.add_axes([0.01, 0.04, 0.63, 0.90]); ax.set_xlim(-4.3, 2.0); ax.set_ylim(-1.9, 1.9); ax.axis("off")
    pos = {"X": (-3.45, 0.0), "hub": (-1.6, 0.0), "slow": (0.0, 0.95), "open-biased": (1.25, 0.0), "fast": (0.0, -0.95)}
    lab = {("hub", "slow"): (-1.05, 0.62), ("slow", "open-biased"): (0.95, 0.62), ("open-biased", "fast"): (0.95, -0.62),
           ("fast", "hub"): (-1.05, -0.62), ("hub", "X"): (-2.53, 0.42)}

    # edges first (under the boxes): width ~ switching events per day
    for (i, j) in g.edges:
        ni = [k for k, v in idx.items() if v == i][0]; nj = [k for k, v in idx.items() if v == j][0]
        a_, c_ = Q[2 * i, 2 * j], Q[2 * i + 1, 2 * j + 1]
        flux = (pi[2 * i] * a_ + pi[2 * i + 1] * c_) * 86400; via_o = pi[2 * i] * a_ * 86400 / flux
        (x1, y1), (x2, y2) = pos[ni], pos[nj]
        ax.plot([x1, x2], [y1, y2], color=GRID, lw=1.5 + 9 * flux / 260, solid_capstyle="round", zorder=1)
        lx, ly = lab[(ni, nj)]
        how = "mostly from open" if via_o > 0.7 else ("mostly from closed" if via_o < 0.3 else "from open and closed")
        ax.text(lx, ly, f"~{flux:.0f} switches/day\n{how}", ha="center", va="center", fontsize=8.5, color=INK2,
                bbox=dict(boxstyle="round,pad=0.2", fc=SURF, ec="none"), zorder=2)

    def box(name, label, color, w=1.0, h=0.62, frozen=False):
        x, y = pos[name]
        c = np.array(matplotlib.colors.to_rgb(color)); s_ = np.array(matplotlib.colors.to_rgb(SURF))
        tint = 0.18 * c + 0.82 * s_                                      # opaque tint: edges stay behind the box
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                                    fc=tint, ec=color, lw=2, zorder=3))
        ax.text(x, y + h / 2 - 0.1, label, ha="center", va="center", fontsize=11.5, fontweight="bold", color=INK, zorder=5)
        return x, y, w, h

    for name in ("slow", "open-biased", "fast"):
        x, y, w, h = box(name, name.capitalize() + " mode", MODE[name], w=1.0 if name != "open-biased" else 1.2)
        to, tc = dwell[name]
        ax.text(x - 0.27, y - 0.06, "open", ha="center", fontsize=9, color=INK2, zorder=5)
        ax.text(x + 0.27, y - 0.06, "closed", ha="center", fontsize=9, color=INK2, zorder=5)
        ax.add_patch(FancyArrowPatch((x - 0.09, y - 0.02), (x + 0.09, y - 0.02), arrowstyle="<|-|>", mutation_scale=9,
                                     color=MODE[name], lw=1.5, zorder=5))
        ax.text(x - 0.27, y - 0.19, fmt_t(to), ha="center", fontsize=9.5, color=INK, zorder=5)
        ax.text(x + 0.27, y - 0.19, fmt_t(tc), ha="center", fontsize=9.5, color=INK, zorder=5)
        ax.text(x, y - h / 2 - 0.1, f"{100 * share[name]:.0f} % of the time", ha="center", fontsize=9, color=INK2, zorder=5)

    x, y, w, h = box("X", "Frozen state", FROZEN, w=1.25)
    ax.text(x, y - 0.06, "no opening or closing", ha="center", fontsize=9, color=INK2, zorder=5)
    ax.text(x, y - 0.19, "lasts ~4 min (mostly open)", ha="center", fontsize=9.5, color=INK, zorder=5)
    ax.text(x, y - h / 2 - 0.1, f"{100 * share['X']:.0f} % of the time", ha="center", fontsize=9, color=INK2, zorder=5)
    x, y, w, h = box("hub", "Transient hub", HUB, w=0.92, frozen=True)
    ax.text(x, y - 0.06, "frozen, short", ha="center", fontsize=9, color=INK2, zorder=5)
    ax.text(x, y - 0.19, "lasts ≲ 3 s", ha="center", fontsize=9.5, color=INK, zorder=5)
    ax.text(x, y - h / 2 - 0.1, f"{100 * share['hub']:.0f} % of the time", ha="center", fontsize=9, color=INK2, zorder=5)
    ax.text(0.0, 0.0, "cycle\nin equilibrium", ha="center", va="center", fontsize=9, color=MUTED, style="italic")
    ax.text(-4.25, 1.85, "Gating modes and how they connect", fontsize=13.5, fontweight="bold", color=INK, va="top")
    ax.text(-4.25, 1.62, "Inside each box: mean open and closed times in that mode.\n"
            "Line width: how often the protein switches along a connection.", fontsize=9.5, color=INK2, va="top")
    ax.text(-4.25, -1.85, "Detailed balance holds around the cycle (no driven circulation; tested). "
            "Frozen states: the protein stays in one conformation.", fontsize=9, color=INK2, va="bottom")

    # right panel: the slow drift of the closing rate
    ax2 = fig.add_axes([0.68, 0.40, 0.30, 0.42], facecolor=SURF)
    rng = np.random.default_rng(3); dt = 1.0; n = 3600; t = np.arange(n) * dt / 60
    z = np.zeros(n); aa = np.exp(-dt / tau[1])
    for k in range(1, n): z[k] = aa * z[k - 1] + np.sqrt(1 - aa * aa) * rng.standard_normal()
    mult = np.exp(sig[1] * z)
    ax2.axhspan(np.exp(-sig[1]), np.exp(sig[1]), color=GRID, alpha=0.6, lw=0)
    ax2.axhline(1, color=MUTED, lw=0.8, ls=":")
    ax2.plot(t, mult, color=INK, lw=1.6)
    ax2.set_xlim(0, 60); ax2.set_xlabel("time (min)", color=INK2, fontsize=10)
    ax2.set_ylabel("closing-rate multiplier", color=INK2, fontsize=10)
    for sp in ("top", "right"): ax2.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax2.spines[sp].set_color(MUTED)
    ax2.tick_params(colors=INK2, labelsize=9)
    ax2.text(60, np.exp(sig[1]) + 0.01, "±1 sd", ha="right", va="bottom", fontsize=8.5, color=INK2)
    fig.text(0.68, 0.93, "Within every mode, the closing rate drifts", fontsize=13.5, fontweight="bold", color=INK, va="top")
    fig.text(0.68, 0.885, "An example of the drift (one simulated hour):", fontsize=9.5, color=INK2, va="top")
    fig.text(0.68, 0.33,
             f"rate of leaving the open state = mode value × e^(σ z(t))\n"
             f"σ = {sig[1]:.2f} [{sig[0]:.2f}, {sig[2]:.2f}]: about ±{100 * (np.exp(sig[1]) - 1):.0f} % (1 sd)\n"
             f"z drifts slowly: correlation time {tau[1] / 60:.1f} min [{tau[0] / 60:.1f}, {tau[2] / 60:.1f}]\n"
             f"the closed→open rate does not drift\n"
             f"the drift is coupled to the conformation so that\n"
             f"detailed balance holds (equilibrium)", fontsize=9.5, color=INK, va="top", linespacing=1.5)
    fig.text(0.68, 0.06, f"Model: ring of four modes + second frozen state + drift (26 parameters).\n"
             f"Max log L {fit['ll']:.1f}; log evidence {logZ:.1f}. [90 % credible intervals]", fontsize=8.5, color=MUTED)
    Path("figures").mkdir(exist_ok=True)
    fig.savefig("figures/model_schematic.png", dpi=150, facecolor=SURF); fig.savefig("figures/model_schematic.svg", facecolor=SURF)
    print("wrote figures/model_schematic.{png,svg}")


if __name__ == "__main__":
    main()
