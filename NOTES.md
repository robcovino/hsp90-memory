# Hsp90 plasmon-ruler traces — analysis notes

Date: 2026-10-11 (updated same day with the Vollmar et al. reading). Folder: `pulling hsp90/`. Scripts in `analysis/`,
figures in `figures/`, fitted models in `results/`.

## Question

Whether the SBI force-spectroscopy framework (Dingeldein et al., PRL 137, 088402, 2026; SBIsmfs) applies to the Hsp90
single-molecule data, and, after reformulation, whether SBI can discriminate memory models for these data beyond what
hidden/aggregated Markov models (AMMs) allow.

## Data

Zenodo 12697153 (Hugel & Sönnichsen): four plasmon-ruler traces of apo yeast Hsp90 (285C–285C dimer between two 60-nm gold
spheres), 6 h each at 10 Hz (216,000 frames), raw scattering counts, no nucleotide. First analysed in Ye et al., Nano Lett. 18,
6633 (2018); re-analysed for memory in Vollmar et al., J. Phys. A (2024), arXiv 2404.16799. These are not pulling data:
zero external force, intensity read-out.

## What the traces contain

- Two-level telegraph signal, open (low) and closed (high), separated by 9–13 % of the mean intensity; within-state spread
  0.7–1.2 % (SNR ≈ 11). Threshold idealisation is error-free (false crossings ~1e-8 per frame).
- Frames are snapshots: 5–13 frames per trace in the middle third of the gap versus ~2000 expected if the exposure
  covered the full 100 ms. Transitions complete between consecutive frames; no transition-path or barrier information.
- Within-state fluctuations are white at 0.1 s (ACF at lag 1 = 0.00 ± 0.005). No intra-well dynamics resolved.
- Peaks Gaussian over four decades (|skew|, |excess kurtosis| ≤ 0.03). Per-dwell mean intensities vary exactly as
  shot noise predicts (within 10 % for dwells ≥ 3 s and ≥ 10 s, all traces): kinetically distinct substates have identical
  intensity to ~3e-4 of <I> (~0.01 nm). Levels stable to 1e-4 over 6 h.
- One reproducible well-shape signal: the open peak is wider than the closed one beyond any intensity-only noise model
  (relative-width ratio closed/open 0.80–0.83 vs 0.94–0.96 for shot noise, 0.88–0.91 for additive noise); excess ≈ 4–5 %
  of the gap, ≈ 0.2 nm rms if the map is linear; white at 0.1 s; identical in slow and fast gating episodes.

Consequence: the continuous-landscape content of these data is one ambiguous number. All other information is in the
sequence of ~25,000 dwell times. The limitation is bandwidth, not statistics.

## Dwell-sequence diagnostics (model-free)

- Successive dwells are correlated: Spearman r = +0.35 between consecutive open dwells, +0.19 between an open dwell and the
  following closed dwell; open–open correlation decays roughly geometrically with lag and reaches zero after ~50–70 cycles.
  Renewal-type memory is excluded.
- The hidden variable slows both directions: after an open dwell > 10 s the next closed dwell averages 4.3 s (vs 1.0 s);
  after a closed dwell > 10 s the next open dwell averages 16 s (vs 3.5 s).
- Modal gating: local transition rates in 2-min windows are bimodal (peak at ~0.09 s⁻¹, hump at ~0.7 s⁻¹). Slow episodes
  average 16 min (max 60), fast episodes 8 min (max 30). Within the slow mode <open> = 18 s, <closed> = 4.0 s, closed
  fraction 0.18; within the fast mode 1.9 s, 1.0 s, 0.35.
- Not environmental: the four rate time series are mutually uncorrelated (pairwise r inside the circular-shift null).
- Not nanoparticle sticking: levels and within-state widths identical between slow and fast episodes to four decimals.
- Time-reversal symmetry of the joint (open, following closed) distribution holds within the shuffle null: consistent
  with detailed balance.
- The four molecules are statistically alike: survival curves overlay; closed fraction 0.22–0.27.

## Aggregated Markov model fits

Exact likelihood of a continuous-time Markov chain observed at 0.1 s through the two-class aggregation, P = expm(Q dt)
(sub-frame excursions handled exactly), pooled over the four traces (25,571 dwells). Modal(n) = n gating modes in a chain,
each with its own open<->closed rates, mode switching from both conformations, detailed balance imposed.

| model | parameters | log L | BIC | ΔBIC vs modal2 |
|---|---|---|---|---|
| modal1 (Markov two-state) | 2 | −110,887.8 | 221,796 | +24,029 |
| modal2 | 7 | −98,848.1 | 197,767 | 0 |
| modal2 without detailed balance | 8 | −98,848.1 | 197,777 | +10 |
| general 2 open + 2 closed, fully connected | 12 | −98,848.1 | 197,818 | +51 |
| modal3 | 12 | −98,360.5 | 196,843 | −925 |
| modal4 | 17 | −98,133.7 | 196,440 | −1,327 |

- modal2: slow mode kOC = 0.050, kCO = 0.227 s⁻¹ (open 20 s, closed 4.4 s); fast mode 0.713, 1.210 s⁻¹ (1.4 s, 0.83 s);
  slow->fast from open 1/855 s, fast->slow 1/270 s; 71 % of time in the slow mode. Matches the model-free numbers.
- Relaxing detailed balance gains ΔlogL = 0.0. The fully connected 4-state model converges to the modal topology.
- Per-trace modal2 fits: rates agree within ~20 % across molecules; sum of per-trace log L exceeds the pooled fit by only 28
  for 21 extra parameters. No molecule-to-molecule heterogeneity at this level.
- modal3 splits the slow mode into a symmetric slow mode (open 13 s, closed 6.2 s, closed fraction 0.32) and an open-biased
  mode (open 30 s, closed 1.1 s, closed fraction 0.04). The hidden variable modulates the two rates differently: a single
  1D modulator with one coupling direction cannot do this.
- modal4 adds a frozen mode (no open<->closed transitions, exit 1/249 s from open, 1/27 s from closed, 12 % occupancy),
  which is what captures open dwells > 100 s and closed dwells > 30 s. The hierarchy is: frozen (minutes) > slow (10–30 s)
  > fast (1 s), i.e. the paper's hour-scale stuck states are the top of a ladder, not a separate phenomenon.

## Predictive checks (`figures/predictive_checks.png`)

Simulations from the fitted generators (4 × 6 h) vs data: modal3 and modal4 reproduce the dwell-lag correlation curve at
all lags, the bimodal local-rate histogram, both survival curves (modal4 including the extreme tails), and the
Chapman–Kolmogorov divergence profile D_τ(t) (τ = 1 s) up to ~30 s. Residual: the D_τ(t) tail at 30–300 s is
underestimated by about a factor of two (data q0 = closed: 0.0041, 0.0022, 0.0011 nats at 30, 100, 200 s; modal4: 0.0022,
0.0004, 0.0002), at borderline significance given four traces. Memory in the propagator sense persists to 200–300 s in the
data, longer than the ~40–50 s reported by Vollmar et al., whose comparison model lacked any minute-scale modulator.

## Conclusions

1. SBIsmfs (harmonic-spring diffusive model, bead + linker observation) is not applicable: no bead dynamics are observed,
   no intra-well dynamics, no transition frames, and the observable is non-Markovian for minutes. It would return the prior
   for the barrier and D and fail its own predictive checks.
2. The dominant memory is discrete, hierarchical modal gating. A finite AMM with 3–4 modes represents it natively, with a
   tractable likelihood, and reproduces every statistic examined except part of the 30–300 s propagator tail.
3. SBI is the tool for model classes outside finite-state Markov: multi-dimensional or GLE-type hidden dynamics, heavy-tailed
   or aging processes, measurement physics by simulation. A 1D continuous modulator is not outside: discretised it is a large
   HMM. The scientifically open question these data pose is "finite ladder of discrete modes vs continuous distribution of
   hidden rates", and it cannot be settled by fit quality (any dwell statistics over 0.1 s–1 h are a finite exponential
   mixture); it needs discriminating statistics (bimodality of local rates — which favours discrete; shape of lag decay;
   rank of the 2D dwell histogram) and perturbation (ATP traces, temperature).
4. The remaining misfit (CK tail, ladder of ever-slower modes) is exactly where a continuous-disorder model could be more
   parsimonious than modal5, modal6, …; that comparison is the SBI project, if any.

Addendum: the seeded fully connected 6-state fit (30 free rates) reaches log L = −98,290 (two restarts agree to 0.3),
i.e. +70 over the 3-mode chain with 18 more parameters; BIC prefers modal3. The modal chain is the right topology at
6 states; the gain from modal4 comes from the frozen mode, not from extra connectivity.

Relaxation spectra of the fitted generators vs the measured autocorrelation: the data ACF of the two-level signal drops
to 0.55 at 1 s, 0.18 at 10 s, plateaus near 0.11 between 30 and 80 s, and reaches zero near 300 s. modal3 contains slow
components at 106 s (10 % of the ACF amplitude) and 314 s (1 %); modal4 at 68 s (8 %), 189 s (1 %) and 384 s (4 %).
The fitted AMMs therefore reproduce the slow correlation component that Vollmar et al. report as t_corr ≈ 90 s and that
their smFRET-rate HMM (t_corr ≈ 14.5 s) missed.

## Vollmar et al., J. Phys. A 57, 365001 (2024) — critical reading and open questions for SBI

What the paper establishes: a model-free Markovianity test (Kullback–Leibler divergence D_τ,q0(t) between the measured
two-point propagator and its Chapman–Kolmogorov reconstruction with hidden degrees of freedom reset at τ), validated on
simulated 4-state networks (uniform rates: no memory; smFRET-inspired rates, equilibrium or driven at 10 k_BT: memory over
~1 decade, t_mem ≈ 6 s), then applied to these four traces: t_corr ≈ 90 s, t_mem ≈ 37–48 s, memory stronger from the closed
state, a shoulder indicating ≥ 2 time scales. Conclusion drawn: HMMs set the memory profile by model choice and must
reproduce it before being trusted; the field needs model-free methods. The paper detects memory and says nothing about its
mechanism.

Critical points:
- The HMM comparison is a straw man. The reference model has no minute-scale modulator because its rates come from ~100 s
  smFRET traces that cannot contain one. The 3–4-mode AMMs fitted here reproduce the lag correlations, the rate
  bimodality, the CK profile to ~30 s and the slow ACF component. A particular HMM failed, not the class.
- "Memory completely lost after ~50 s" is not supported: in our computation D is ~13 % of its peak at 100 s and ~6 % at
  200 s, about two across-trace standard deviations above zero. The exponential-tail t_mem ignores the shoulder the
  authors themselves note.
- The bootstrap (random deletion of 20 % of frames) underestimates uncertainty for data with a 90 s correlation time; a
  block bootstrap or the across-trace spread is the honest error bar.
- The nanoparticle exclusion ("diffusion time 0.1 µs") rules out hydrodynamic coupling only, not transient surface or
  PEG interactions. Mode-invariance of levels and widths (sticking test above) is a stronger argument, still not a proof.
- Correct and consequential: memory depends on q0, so no convolution-type kernel can describe the observable; any reduced
  description must be state-dependent. Stated and left open.

Open questions an SBI approach could settle (ranked):

1. Mechanism class of the memory — finite ladder of discrete hidden states, continuous slow coordinate (dynamic
   disorder; their refs 49–51), or memory-dependent friction (ref 52). The CK test is blind to mechanism by design; HMM
   likelihoods represent only the first class. The AMM ladder here (+488, +227 log-units for modes 3 and 4) is what a
   finite-exponential approximation of a continuum looks like, while the bimodal local-rate histogram argues for
   discreteness at the top level. SBI: simulators for modal-K, rate-modulation by a diffusing y on U(y), and a
   Markovian-embedded GLE gating the transition; summaries = survival curves, lag-decay shape, rank of the 2D dwell
   histogram, local-rate histogram, D_τ(t); Bayesian model comparison. Go/no-go: simulation study at N ≈ 25,000 dwells.
   Data in hand.
2. Coupling of the hidden coordinate(s) to the two rates — the "kinetic landscape". modal3 (closed fraction 0.32 vs 0.04
   in the two slow modes) and the q0-dependence of D both say opening and closing are modulated differently: ≥ 2 hidden
   dimensions or state-specific couplings. SBI with spline-parameterised k_OC(y), k_CO(y) and U(y) infers the joint
   stationary distribution of the rate pair and its relaxation spectrum (the hidden-space analogue of G0(x) in the PRL).
   Identifiable only up to monotone reparameterisation of y: fix a gauge.
3. ATP: is the memory structure driven, and can dissipation be inferred from a two-level observable? The paper shows its
   test cannot separate equilibrium hidden states from driving and flags the ATPase rate (1/min), which coincides with the
   apo mode-switching time scale (4–14 min episodes). In apo the modes are not hydrolysis-driven; under ATP the nucleotide
   state may couple to them. Inferring affinity/entropy production from coarse-grained traces is an active problem (their
   refs 7, 13, 16); a driven-network simulator with SBI gives a model-based estimator with calibrated uncertainty, to be
   set against model-free lower bounds. Data: the 33,000-transition ATP set from Ye et al., not on Zenodo.
4. Is smFRET "static heterogeneity" the plasmon-ruler modal gating? A 100 s FRET trace cannot see a 10-min mode switch; it
   sees fast molecules and slow molecules. One latent model with two observation models (10-Hz snapshots for hours;
   photon-by-photon FRET for seconds) is awkward for HMM pipelines and natural for simulator-based inference — the
   machinery of the smFRET SBI paper. Also settles whether the slow modes are intrinsic to Hsp90 (same modes without gold
   and tether). Data: Hugel's smFRET traces under matching conditions (Schmid 2016; Vollmar 2024 Nat. Commun.).
5. Ergodicity and aging across the population. Ye et al. claim non-ergodicity in 24 h traces; the four 6 h traces are one
   population. Whether an hour-scale tier exists (hierarchy → effective continuum, or trap model with aging) is open; aging
   processes have no stationary likelihood, SBI handles them, and amortised per-molecule posteriors over the hundreds of
   rulers recorded in parallel give a hierarchical test. Data: 24 h traces and the parallel-ruler ensemble.
6. Relation between CK memory time and GLE kernel time (their stated open question). Mostly theory (compute D_τ(t) for
   embedded GLEs); for a dichotomous observable the GLE is the wrong object; becomes an inference question only with kHz
   data, where the PRL machinery re-enters.

Not settled by SBI: whether the modulator is protein-intrinsic rather than surface-related (needs particle-size,
linker-length and no-gold controls); biological-relevance claims about client interaction (needs cochaperone experiments).

Positioning: the paper frames model-free vs model-biased. SBI dissolves the dichotomy — the model-free statistics become
the summaries the posterior must reproduce, with explicit priors over model classes and calibrated uncertainty. Questions
1–2 are doable now and cheap to kill if unidentifiable; 4 has a structural advantage no HMM pipeline has and builds on
existing work; 3 is the most interesting physics and depends on the ATP data.

## Phase 0 — identifiability study (2026-10-11): GO

Question: within each gating mode, are the open<->closed rates fixed (discrete modes) or spread (continuous disorder)?
Nested simulator (`simulators.simulate_modal_ou`): gating modes and their switching from the fitted modal3 generator;
within each mode ln k_OC = ln k_OC(m) + σ z1, ln k_CO = ln k_CO(m) + σ z2, with z1, z2 unit-variance OU processes of
relaxation time τ_w and correlation ρ. σ = 0 is the discrete model. Grid σ ∈ {0.25, 0.5, 1}, τ_w ∈ {1, 10, 100} s,
ρ ∈ {0, 1}; 50 replicate measurements (4 × 6 h at 10 Hz) per setting; 126 summary statistics
(`summaries.summary_vector`): dwell quantiles, dwell-lag correlations, local-rate histograms (1- and 2-min windows),
within-episode Fano factors, Chapman–Kolmogorov profile, 2D dwell-histogram spectrum, ACF. Separability = 5-fold
cross-validated AUC of logistic regression / random forest. Figure: `figures/phase0_identifiability.png`.

Results:
- Sanity: nested simulator at σ = 0 vs exact AMM: AUC 0.59 / 0.45 (indistinguishable, as required).
- Against the discrete model with the same mode rates: AUC = 1.00 for all σ ≥ 0.5 and for σ = 0.25 with τ_w ≥ 10 s or
  ρ = 1; the only marginal case is σ = 0.25, τ_w = 1 s, ρ = 0 (AUC 0.82 / 0.77) — fast disorder averages out within a
  dwell, as expected.
- The fair test — modal3 refitted by maximum likelihood to each disordered dataset, then compared to it: AUC 1.00 for
  all σ ≥ 0.5; 0.97–1.00 for σ = 0.25 with τ_w ≥ 10 s; 0.85 / 0.68 for σ = 0.25, τ_w = 1 s. Refitting modal4 (8 states)
  to σ = 0.25, τ_w = 10 s still gives 0.97 / 0.90, and 1.00 / 0.99 for σ = 0.5. An extra discrete mode does not absorb
  the continuum at this data volume.
- What carries the separation: the intermediate bins (0.3–0.7 s⁻¹) of the local-rate histogram, the Fano factor of
  transition counts within fast episodes, the low quantiles of the dwell distributions, and, for τ_w = 100 s, the ACF
  at 5–30 s and the dwell-lag correlations.
- Where the data sit (preview, not inference): the data have more within-mode dispersion than either discrete model —
  Fano in fast episodes 5.1 vs 4.2 ± 0.4 (modal4) and 4.4 ± 0.5 (modal3); the 0.3–0.6 s⁻¹ rate bins are over-populated
  by 2–3.5 sd relative to modal4; on these statistics the data are closest to σ ≈ 0.25–0.5 with τ_w ≈ 10–100 s.

Conclusion: within-mode continuous disorder with σ ≥ 0.25 and τ_w ≥ 10 s is identifiable against discrete modal models
of 3 and 4 modes at the experimental data volume, and the data lean in that direction. Phase 1 is justified.
Caveats: the continuous model inherits the modal3 top-level structure, which the data exceed (frozen mode, CK tail);
50 replicates per setting (AUC ± ~0.05 near 0.85); classifier separability is not posterior precision — Phase 1 gives
the posterior widths on (σ, τ_w, ρ), and the coupling question (ρ, or separate σ_OC, σ_CO) is where question 2 starts.

## Next steps

- Phase 1 (repository `hsp90-memory/`, to be continued in Claude Code): NPE over the nested model with priors on the
  mode rates (from modal3/modal4 posteriors), σ_OC, σ_CO, τ_w, ρ; ~1e5 simulations; simulation-based calibration;
  posterior-predictive checks on the real traces; model comparison modal-K vs modal-OU vs (optionally) GLE-gated.
- Ask Hugel, in one email: the ATP traces (33,000 transitions); the smFRET traces under matching conditions; the 24 h
  traces and parallel-ruler ensemble; acquisition details (exposure time; whether the 10-Hz series is decimated from a
  faster acquisition — the absence of straddle frames suggests so). kHz bursts would make intra-well dynamics and the
  diffusive model meaningful.

## Files (repository `hsp90-memory/`)

- `src/hsp90memory/` — `io.py` (loading, idealisation), `amm.py` (AMM likelihood, fitting), `simulators.py` (AMM and
  nested modal-OU simulators), `summaries.py` (summary statistics incl. CK divergence). `tests/` — 7 tests, `pytest`.
- `scripts/` — `run_fits.py`, `run_general_seeded.py`, `predictive_checks.py`, `phase0_simulate.py`, `phase0_refit.py`,
  `phase0_analyse.py`; `scripts/diagnostics/` — the model-free diagnostics (`plot_trace.py`, `within_state.py`,
  `dwell_memory.py`, `sticking_test.py`).
- `results/*.json` — fitted generators Q, stationary distributions, log L, BIC; `results/phase0/` — summary vectors
  per setting (npz), refit log-likelihoods, `phase0_results.json`.
- `figures/` — `1Thermal_overview.png`, `within_state.png`, `predictive_checks.png`, `phase0_identifiability.png`.
- Data: `data/1Thermal.txt` … `4Thermal.txt` (not committed; Zenodo 12697153).

## Running on the VM (de.NBI) — hand-off to Claude Code

Repository: https://github.com/robcovino/hsp90-memory (transfer to covinolab later if wanted; links keep working).

On a fresh Ubuntu VM, one line, then follow what it prints:

    curl -fsSL https://raw.githubusercontent.com/robcovino/hsp90-memory/main/setup_vm.sh | bash

It installs packages and Claude Code, clones this repository, creates the Python environment, downloads the four traces
from Zenodo, and runs the tests. Then `claude` once to sign in (it prints a link to open on the laptop), then
`bash ~/hsp90-memory/start_remote.sh`, which starts Remote Control inside tmux. The session "hsp90-memory" then appears
at https://claude.ai/code and in the Claude phone app; everything from here on is done from there.

First message to give the session (copy verbatim):

> Read CLAUDE.md and NOTES.md before anything else. Then, in this order: (1) Rewrite the frame loop of
> `simulators.simulate_modal_ou` with numba, keeping the numpy version as the reference implementation, and add a test
> that both produce the same summary statistics on a fixed seed within replicate tolerance. Report the speed-up on one
> 6-h trace. (2) Write a checkpointed Phase 1 simulation driver: prior over the modal3 mode rates (log-normal around the
> fitted values, factor 2), sigma_OC and sigma_CO uniform on [0, 1], tau_w log-uniform on [1, 300] s, rho uniform on
> [-1, 1]; batches of 1000 simulated measurements (4 x 6 h each) written to `results/phase1/batch_XXXX.npz` with
> parameters and `summaries.summary_vector`, resumable from the last complete batch. (3) Run a 10^4-simulation smoke
> test, report wall time and the distribution of summaries, and stop. Do not launch the full 10^5 run until I say so.
> Update NOTES.md with what you did and the timings.

After approving the smoke test: "proceed with Phase 1 step 2" (full run, NPE training, simulation-based calibration,
predictive checks on the real traces — see CLAUDE.md, Next steps).
