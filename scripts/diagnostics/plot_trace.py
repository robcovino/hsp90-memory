"""Overview plot of one plasmon-ruler trace (Zenodo 12697153, 10 Hz, 6 h).
Usage: python3 analysis/plot_trace.py 1   -> figures/1Thermal_overview.png
"""
import sys, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

k = sys.argv[1] if len(sys.argv) > 1 else "1"
dt = 0.1
x = np.loadtxt(f"data/{k}Thermal.txt")
I = x / x.mean()                      # relative intensity, as in Ye et al. Fig. 2
t = np.arange(I.size) * dt

fig = plt.figure(figsize=(12, 9))
gs = fig.add_gridspec(3, 4, height_ratios=[1.2, 1, 1], hspace=0.45, wspace=0.35)

# (a) full trace, hours
ax = fig.add_subplot(gs[0, :3])
ax.plot(t/3600, I, lw=0.3, color="C0", rasterized=True)
ax.set_xlabel("time (h)"); ax.set_ylabel("I / <I>")
ax.set_title(f"{k}Thermal.txt — full 6 h at 10 Hz (216,000 points)")
# (a') histogram
axh = fig.add_subplot(gs[0, 3], sharey=ax)
axh.hist(I, bins=150, orientation="horizontal", color="C0", density=True)
axh.set_xlabel("density"); axh.set_title("distribution")
plt.setp(axh.get_yticklabels(), visible=False)

# (b) 10-min window
t0 = 2*3600 + 15*60; w = 600
sl = slice(int(t0/dt), int((t0+w)/dt))
ax = fig.add_subplot(gs[1, :])
ax.plot(t[sl]-t0, I[sl], lw=0.6, color="C0")
ax.set_xlabel(f"time (s) from t = {t0/3600:.2f} h"); ax.set_ylabel("I / <I>")
ax.set_title("10-min window")

# (c) 1-min window with data points
w2 = 60; t1 = t0 + 240
sl2 = slice(int(t1/dt), int((t1+w2)/dt))
ax = fig.add_subplot(gs[2, :])
ax.plot(t[sl2]-t1, I[sl2], lw=0.6, color="0.6")
ax.plot(t[sl2]-t1, I[sl2], ".", ms=3, color="C0")
ax.set_xlabel(f"time (s) from t = {t1/3600:.2f} h"); ax.set_ylabel("I / <I>")
ax.set_title("1-min window, individual 0.1-s frames")

out = f"figures/{k}Thermal_overview.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print("wrote", out)
