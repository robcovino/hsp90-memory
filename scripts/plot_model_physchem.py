"""Chemical-physics rendering of the best model (ringF1 + equilibrium OU on kOC, exact-likelihood ML fit).

(a) kinetic network: 5 gating modes x {open, closed}; node area ~ equilibrium population; edge width ~ log of the
    one-way equilibrium flux pi_i k_ij; labels give mean waiting times 1/k.
(b) disconnectivity graph of the 10-state network (z = 0 slice): state free energies G_i = -ln pi_i and effective
    transition-state free energies G_ij = G_i - ln(k_ij / k0) (detailed balance makes G_ij = G_ji). k0 = 1e6 s^-1 is an
    arbitrary attempt frequency: only differences between barriers are meaningful.
(c) free energy along the hidden slow coordinate z (exact marginals of the full model): open and closed surfaces,
    and the opening-rate modulation kOC(z) = kOC exp(sigma z).
Usage: python scripts/plot_model_physchem.py  -> figures/model_physchem.{png,pdf}
"""
import json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from hsp90memory import amm

K0 = 1e6
NAMES = {0: "H", 1: "S", 2: "B", 3: "F", 4: "Z"}
LONG = {0: "hub", 1: "slow", 2: "open-biased", 3: "fast", 4: "frozen"}
CO, CC = "#1f5fbf", "#d1622b"       # open (low intensity), closed (high intensity)
INK, MUTED = "#222222", "#777777"

d = json.load(open("results/exact/fit_ringF1_oueq_global.json"))
th = np.array(d["theta"]); m = amm.modal_ou(5, n_grid=31, couple="eq", topology="ringF1")
nb = m.base.n_par; sig, tau = np.exp(th[nb]), np.exp(th[nb + 1])
Q0 = m.base.Q(th[:nb]); p0 = amm.stationary(Q0)                    # z = 0 slice (detailed balance)
Q = m.Q(th); pi = amm.stationary(Q).reshape(m.n_grid, 10)
pmode = pi.sum(0)                                                   # exact populations, marginal over z

def fmt_t(t):
    if t < 60: return f"{t:.2g} s" if t < 10 else f"{t:.0f} s"
    if t < 3600: return f"{t/60:.2g} min" if t < 600 else f"{t/60:.0f} min"
    if t < 86400: return f"{t/3600:.2g} h"
    return f"{t/86400:.2g} d"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.linewidth": 0.8,
                     "axes.edgecolor": INK, "xtick.color": INK, "ytick.color": INK})
fig = plt.figure(figsize=(15, 5.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.0, 0.95], wspace=0.28)

# ------------------------------------------------------------------ (a) network
ax = fig.add_subplot(gs[0]); ax.set_aspect("equal"); ax.axis("off")
pos = {0: (0.0, 0.0), 1: (-2.3, 2.3), 2: (0.0, 4.6), 3: (2.3, 2.3), 4: (0.0, -2.9)}   # ring H-S-B-F-H, Z below H
dx = 0.42
def xy(s):                                     # state index -> coordinates (O left, C right)
    x, y = pos[s // 2]; return (x - dx, y) if s % 2 == 0 else (x + dx, y)
flux = {(i, j): p0[i] * Q0[i, j] for i in range(10) for j in range(10) if i != j and Q0[i, j] > 1e-12}
lf = np.log10(np.array(list(flux.values()))); lo, hi = lf.min(), lf.max()
def width(f): return 0.4 + 3.6 * (np.log10(f) - lo) / (hi - lo)
drawn = set()
for (i, j), f in sorted(flux.items(), key=lambda kv: kv[1]):
    if (j, i) in drawn or i // 2 == j // 2: continue
    drawn.add((i, j))
    col = CO if i % 2 == 0 else CC
    (x1, y1), (x2, y2) = xy(i), xy(j)
    ax.plot([x1, x2], [y1, y2], color=col, lw=width(f), alpha=0.55, solid_capstyle="round", zorder=1)
for k in range(5):                                                   # mode boxes and O<->C edges
    x, y = pos[k]
    ax.add_patch(FancyBboxPatch((x - 0.85, y - 0.42), 1.70, 0.84, boxstyle="round,pad=0.02,rounding_size=0.12",
                                fc="#f6f6f4", ec=MUTED, lw=0.8, zorder=2))
    kOC, kCO = Q0[2 * k, 2 * k + 1], Q0[2 * k + 1, 2 * k]
    if kOC > 1e-5:
        ax.annotate("", xy=(x + dx - 0.13, y + 0.06), xytext=(x - dx + 0.13, y + 0.06),
                    arrowprops=dict(arrowstyle="-|>", lw=0.9, color=INK, shrinkA=0, shrinkB=0), zorder=4)
        ax.annotate("", xy=(x - dx + 0.13, y - 0.06), xytext=(x + dx - 0.13, y - 0.06),
                    arrowprops=dict(arrowstyle="-|>", lw=0.9, color=INK, shrinkA=0, shrinkB=0), zorder=4)
        ax.text(x, y + 0.5, f"⟨open⟩ {fmt_t(1/kOC)}  ·  ⟨closed⟩ {fmt_t(1/kCO)}", ha="center", va="bottom",
                fontsize=7.5, color=INK, zorder=5, bbox=dict(fc="white", ec="none", pad=0.5))
    else:
        pass
    for s, col in ((2 * k, CO), (2 * k + 1, CC)):
        r = 0.05 + 0.24 * np.sqrt(pmode[s] / pmode.max())
        ax.add_patch(plt.Circle(xy(s), r, fc=col, ec="white", lw=1.5, zorder=6))
    ax.text(x, y - 0.52, f"{NAMES[k]}  {LONG[k]}  ·  {100*(pmode[2*k]+pmode[2*k+1]):.0f} %" + ("" if Q0[2*k, 2*k+1] > 1e-5 else "  ·  no O⇄C"), ha="center", va="top",
            fontsize=8, color=INK, zorder=5, bbox=dict(fc="white", ec="none", pad=0.5))
for a, b in [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4)]:               # dominant switching times, open channel
    ta, tb = 1 / Q0[2 * a, 2 * b], 1 / Q0[2 * b, 2 * a]
    (xa, ya), (xb, yb) = pos[a], pos[b]; xm, ym = (xa + xb) / 2, (ya + yb) / 2
    off = {(0, 1): (-1.3, -0.3), (1, 2): (-1.05, -0.15), (2, 3): (1.05, -0.15), (3, 0): (1.3, -0.3), (0, 4): (0.95, -0.2)}[(a, b)]
    ax.text(xm + off[0], ym + off[1], f"{NAMES[a]}→{NAMES[b]} {fmt_t(ta)}\n{NAMES[b]}→{NAMES[a]} {fmt_t(tb)}",
            ha="center", va="center", fontsize=7, color=MUTED)
ax.plot([], [], "o", color=CO, ms=8, label="open"); ax.plot([], [], "o", color=CC, ms=8, label="closed")
ax.legend(loc="upper left", frameon=False, fontsize=8, bbox_to_anchor=(-0.02, 1.0))
ax.text(-3.4, -3.75, "node area ∝ population;  line width ∝ log one-way flux π$_i$k$_{ij}$;\n"
        "switching times via the open channel; mode-internal rates at z = 0", fontsize=7, color=MUTED, va="top")
ax.set_xlim(-3.5, 3.5); ax.set_ylim(-4.3, 5.6)
ax.set_title("a   Gating-mode network (ring + frozen branch)", loc="left", fontsize=10, fontweight="bold")

# ------------------------------------------------------------------ (b) disconnectivity graph
ax = fig.add_subplot(gs[1])
G = -np.log(p0); G -= G.min()
edges = []
for (i, j), f in flux.items():
    if i < j or (j, i) not in flux:
        edges.append((G[i] - np.log(Q0[i, j] / K0), i, j))
edges.sort()
parent = {i: i for i in range(10)}; clusters = {i: dict(members=[i], top=None, kids=[], bottom=G[i]) for i in range(10)}
def find(x):
    while parent[x] != x: x = parent[x]
    return x
nid = 10
for Gts, i, j in edges:
    ri, rj = find(i), find(j)
    if ri == rj: continue
    clusters[nid] = dict(members=clusters[ri]["members"] + clusters[rj]["members"], top=None, kids=[ri, rj],
                         merge=Gts, bottom=min(clusters[ri]["bottom"], clusters[rj]["bottom"]))
    clusters[ri]["top"] = Gts; clusters[rj]["top"] = Gts
    parent[ri] = nid; parent[rj] = nid; parent[nid] = nid; nid += 1
root = find(0); xs = {}; counter = [0]
def layout(c):
    if c < 10:
        xs[c] = counter[0]; counter[0] += 1; return xs[c]
    a, b = clusters[c]["kids"]
    if clusters[a]["bottom"] > clusters[b]["bottom"]: a, b = b, a
    xs[c] = 0.5 * (layout(a) + layout(b)); return xs[c]
layout(root)
for c, cl in clusters.items():
    lo_ = G[c] if c < 10 else cl["merge"]; hi_ = cl["top"] if cl["top"] is not None else lo_ + 1.5
    col = (CO if c % 2 == 0 else CC) if c < 10 else INK
    ax.plot([xs[c], xs[c]], [lo_, hi_], color=col if c < 10 else INK, lw=2.0 if c < 10 else 1.2, solid_capstyle="butt")
    if c >= 10:
        a, b = cl["kids"]; ax.plot([xs[a], xs[b]], [cl["merge"], cl["merge"]], color=INK, lw=1.2)
    if c < 10:
        ax.text(xs[c], G[c] - 0.35, f"{'OC'[c % 2]}$_{{{NAMES[c // 2]}}}$", ha="center", va="top", fontsize=8.5,
                color=CO if c % 2 == 0 else CC)
ax.set_ylabel(r"free energy  $G/k_\mathrm{B}T$")
ax.set_xticks([]); ax.spines[["top", "right", "bottom"]].set_visible(False)
ax.set_ylim(-1.2, max(cl.get("merge", 0) for cl in clusters.values()) + 2.0)
ax.set_title("b   Disconnectivity graph", loc="left", fontsize=10, fontweight="bold")
ax.text(0.02, 0.98, f"barriers: $k = k_0 e^{{-\\Delta G^\\ddagger}}$, $k_0$ = 10$^6$ s$^{{-1}}$\n(only differences meaningful)",
        transform=ax.transAxes, ha="left", va="top", fontsize=7, color=MUTED)

# ------------------------------------------------------------------ (c) hidden slow coordinate
ax = fig.add_subplot(gs[2]); z = m.z
pO = pi[:, 0::2].sum(1); pC = pi[:, 1::2].sum(1)
h = z[1] - z[0]
GO = -np.log(pO / h); GC = -np.log(pC / h); g0 = min(GO.min(), GC.min()); GO -= g0; GC -= g0
zz = np.linspace(z.min(), z.max(), 400)
fO = np.polyfit(z, GO, 2); fC = np.polyfit(z, GC, 2)
ax.plot(zz, np.polyval(fO, zz), color=CO, lw=2, label="open (all modes)")
ax.plot(zz, np.polyval(fC, zz), color=CC, lw=2, label="closed (all modes)")
ax.plot(z, GO, "o", color=CO, ms=3.5, mec="white", mew=0.6); ax.plot(z, GC, "o", color=CC, ms=3.5, mec="white", mew=0.6)
zO, zC = -fO[1] / (2 * fO[0]), -fC[1] / (2 * fC[0])
for zm, f, col in ((zO, fO, CO), (zC, fC, CC)):
    ax.plot([zm, zm], [np.polyval(f, zm), np.polyval(f, zm) + 0.6], color=col, lw=0.8, ls=":")
ax.annotate("", xy=(zC, np.polyval(fC, zC) - 0.25), xytext=(zO, np.polyval(fC, zC) - 0.25),
            arrowprops=dict(arrowstyle="->", lw=0.9, color=INK))
ax.text(zC / 2, -0.55, f"Δz = σ = {zC - zO:.2f}", ha="center", va="top", fontsize=8)
ax.set_xlim(-3.2, 3.2); ax.set_ylim(-2.6, 8.0)
ax.set_xlabel("hidden slow coordinate  z  (units of its sd)"); ax.set_ylabel(r"free energy  $G/k_\mathrm{B}T$")
ax.spines[["top", "right"]].set_visible(False); ax.legend(loc="upper left", frameon=False, fontsize=8)
ins = ax.inset_axes([0.72, 0.07, 0.26, 0.18])
ins.plot(zz, np.exp(sig * zz), color=INK, lw=1.4); ins.axhline(1, color=MUTED, lw=0.6, ls=":")
ins.set_xlim(-3, 3); ins.set_ylim(0.4, 2.0); ins.set_yticks([0.5, 1, 1.5, 2]); ins.tick_params(labelsize=7)
ins.set_xlabel("z", fontsize=7, labelpad=1); ins.set_ylabel(r"$k_\mathrm{OC}(z)/k_\mathrm{OC}(0)$", fontsize=7, labelpad=1)
ins.spines[["top", "right"]].set_visible(False)
ax.text(0.27, 0.83, f"σ = {sig:.2f},  τ$_w$ = {tau:.0f} s\nclosing tilts z (detailed balance)", transform=ax.transAxes,
        fontsize=8, va="top", color=INK)
ax.set_title("c   Slow modulation of the opening rate", loc="left", fontsize=10, fontweight="bold")

fig.suptitle(f"Hsp90 (apo) conformational kinetics: best model ringF1 + slow disorder  "
             f"(exact-likelihood ML, log L = {d['ll']:.1f}; 4 × 6 h, 25,571 dwells)", fontsize=10.5, x=0.01, ha="left")
for ext in ("png", "pdf"):
    fig.savefig(f"figures/model_physchem.{ext}", dpi=200, bbox_inches="tight")
print("ok", f"sigma={sig:.3f} tau={tau:.0f}", "zO,zC", zO, zC)
