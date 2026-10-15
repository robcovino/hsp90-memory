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

## Phase 1 — simulation (2026-10-12)

Simulator speed. `simulators.simulate_modal_ou_nb` runs the frame loop of the nested modal-OU model in numba (loop
over replicates and frames, numba's own RNG). `simulate_modal_ou` (numpy, vectorised over replicates) is kept as the
reference. Both now take sigma as a scalar or as a pair (sigma_OC, sigma_CO). One 6-h trace (216,000 frames, modal3
modes, sigma = (0.5, 0.3), tau_w = 10 s): numpy 12.2 s, numba 0.045 s (~270x); one measurement of 4 x 6 h: 12.0 s vs
0.22 s. The numpy version costs the same for 1 or 4 replicates because it is bound by per-frame overhead.
`tests/test_simulators.py` holds the equivalence test: 20 replicate measurements of 4 x 1 h each, sigma =
(0.6, 0.4), tau_w = 20 s, rho = 0.5. All 113 finite summaries agree (max |z| 3.6, 1.8 % beyond 2 sd). The same test
rejects the numba simulator with sigma = 0 (max |z| 9.8), and it also separated sigma = (0.3, 0.2) (6.5).
After the speed-up, `summaries.summary_vector` is the bottleneck: 0.33 s of the 0.45 s per simulated measurement.

Prior (`scripts/phase1_simulate.py`; settings in `results/phase1/settings.json` and in every batch file):
- all 12 free log-rates of modal3 (within-mode kOC, kCO and the 6 switching rates): N(ML value, sd ln 2), i.e.
  log-normal with factor 2. Q is rebuilt with `amm.modal(3)`, so detailed balance holds for every draw;
- sigma_OC, sigma_CO ~ U[0, 1]; tau_w ~ log-U[1, 300] s; rho ~ U[-1, 1].
One simulation is one measurement of 4 x 6 h at 10 Hz reduced to the 126-entry `summary_vector`. Batches of 1000 go to
`results/phase1/batch_XXXX.npz`, written atomically. Parameters and seeds of batch b depend only on (seed, b), so a
rerun skips complete batches and reproduces any missing one bit for bit (checked). A worker uses ~97 MB peak RSS.
Batch files are not committed (`.gitignore`; ~1.2 MB each).

Smoke test, 10^4 simulations (batches 0–9; they are the first 10 batches of the full run):
- 1363 s wall on 4 workers (0.454 s per simulation per worker; 5/50/95 %: 0.35/0.44/0.61 s).
- So 10^5 takes ~2.9 h on this VM.
- Non-finite entries: 12 structural (CK lags <= tau = 1 s, undefined by construction). fano_1 (fast episodes) is
  missing in 1.7 % of draws, those with no fast episodes. sv2d_0 = 1 by normalisation and carries no information.
  The NPE drops the structural and constant entries and imputes fano_1 with an indicator.
- The data sit inside the prior predictive on every statistic, mostly between the 30th and 75th percentiles
  (`results/phase1/smoke_report.json`, `figures/phase1_smoke.png`). The tails are the 99.9 % open-dwell quantile
  (percentile 0.97: the long-open tail that modal3, lacking the frozen mode, under-produces), the high-rate bins of the
  local-rate histograms (0.88–0.92) and the CK tail from q0 = closed at 50–140 s (0.74–0.79). The prior covers
  the data. No summary is degenerate apart from sv2d_0.

Full run, 10^5 simulations (batches 0–99; seed 20261012): 90 batches after the smoke test in 10,950 s wall
(3.04 h; 0.46–0.52 s per simulation per worker), ~97 MB per worker throughout.

## Phase 1 — NPE, calibration and posterior for the data (2026-10-12, tag v1)

Setup (`scripts/phase1_npe.py --tag v1`; `results/phase1/npe_v1.{pt,json}`, `npe_v1_sbc.npz`,
`npe_v1_posterior_data.npz`, logs `npe_v1.log`, `npe_v1_eval.log`; figures `figures/phase1_npe_v1_*.png`).
- Trained on 99,000 simulations, with batch 0 (1000 simulations) held out for calibration.
- Parameters: the 12 modal3 log-rates, sigma_OC, sigma_CO, log10 tau_w and rho.
- Summaries: 115 = 126 − 12 structural CK NaNs − sv2d_0, plus 2 imputation indicators.
- Network: NSF (5 transforms, 128 hidden) on an MLP embedding (3 x 256 -> 64), z-scored inputs.
- Training: converged after 99 epochs in 3021 s on 4 threads; best validation loss −5.94 (held-out mean log q of the
  true parameters +6.02).
- The first evaluation was OOM-killed at 7.5 GB: sbi's `log_prob_batched` with leakage normalisation samples per
  test point. The evaluation now draws posterior samples in chunks, computes ranks and TARP from them, and resumes
  from the saved estimator. It takes 163 s and 1.3 GB.

Calibration (1000 held-out simulations x 1000 posterior samples):
- Joint: TARP ATC −0.011 (KS p = 1.00).
- Marginal coverage of central intervals within 0.01–0.03 of nominal for all 16 parameters (50/80/90/95 %).
- SBC rank KS p > 0.09 for 12 of 16. The exceptions:
  - log10 tau_w (p = 0.002): mildly overconfident, coverage 0.74/0.87/0.93 at 80/90/95 %; posterior medians shrink
    towards the centre for truths near either prior edge (bias +0.21 below 3 s, −0.16 above 100 s).
  - sigma_CO (p = 0.02): cannot resolve sigma_CO < 0.1 (bias +0.10 there).
  - log_a0, log_b0 (p = 0.03, 0.01): the 1/(2–5 days) O_0<->O_1 rates, which the data do not inform.
- Posterior contraction (1 − var_post/var_prior), mean over the held-out set:
  - within-mode rates 0.72–0.98;
  - switching rates 0.30–0.60;
  - sigma_OC 0.89, sigma_CO 0.88, log10 tau_w 0.91, rho 0.66.

Posterior for the data (median, 90 % interval):
- sigma_OC = 0.38 [0.29, 0.45]; P(sigma_OC < 0.1) = 0.003.
- sigma_CO = 0.14 [0.03, 0.26]; P(sigma_CO < 0.1) = 0.34.
- P(both < 0.1) = 2e-4: within the modal3 skeleton the discrete model (sigma = 0) is excluded.
- tau_w = 180 s [37, 291 s]: piled against the prior bound (300 s); P(tau_w > 100 s) = 0.76.
- rho = 0.40 [−0.51, 0.95] (P(rho > 0) = 0.77): weakly informed.
- Within-mode rates (geometric medians) stay close to the modal3 ML values, except kOC in the open-biased mode:
  0.021 [0.016, 0.032] vs 0.034 s^-1 (longer open dwells there). That marginal has a small second mode (4 %) near
  the ML value.
- Leakage outside the prior support for the data: 7.6 % (rejected).

The disorder is in kOC (the rate of leaving the open state) and is slow: minutes, comparable to mode switching
(episodes 8–16 min). kCO is compatible with no disorder. This fits the q0 dependence of the CK memory (stronger from
closed, i.e. after a long open dwell) and the earlier finding that the hidden variable modulates the two rates
differently.

Posterior-predictive checks (`scripts/phase1_ppc.py --tag v1`; `results/phase1/ppc_v1.{json,npz}`,
`figures/phase1_ppc_v1.png`; 253 s). Setup: 200 measurements from 200 posterior draws, against 200 measurements each
from modal3 and modal4 at their ML generators (exact AMM simulator). For each of the 113 finite statistics, the data's
predictive percentile and z score:

| model | data outside 95 % band | \|z\| > 3 | sum z^2 |
|---|---|---|---|
| modal-OU (posterior predictive) | 12 / 113 | 2 | 159 |
| modal4 (ML plug-in) | 18 / 113 | 12 | 479 |
| modal3 (ML plug-in) | 40 / 113 | 29 | 1100 |

- modal-OU reproduces what modal4 misses:
  - the CK tail at 50–200 s, from both q0 (modal4 z = 4–7; modal-OU z ≈ 2.3, data at the upper edge of the band);
  - the ACF at 60–200 s (modal4 z = 4–5; modal-OU 2.0);
  - the fast-episode Fano factor, partly: data 5.1; modal-OU 6.1 [4.8, 7.6] (5–95 %) overshoots it, while modal3
    4.4 [3.7, 5.3] and modal4 4.2 [3.5, 5.1] undershoot it. The data lie between the continuous and discrete
    predictions.
- modal-OU misses what modal4's frozen mode captures:
  - the extreme open tail (99.9 % open-dwell quantile z = +4.2): the data have open dwells up to ~500 s, while
    modal-OU's open survival falls off above ~250 s;
  - too much weight at 50–150 s open (90 % quantile z = −3.0);
  - the 99.9 % closed quantile (z = +2.7).
  The slow kOC disorder is doing part of the frozen mode's job, and does it with the wrong shape.
- The comparison is not symmetric. modal-OU's bands include parameter uncertainty, while the modal3/4 bands are
  plug-in, so the counts favour modal-OU somewhat. A fair model comparison needs simulation-based posteriors (or
  Laplace-approximated likelihood posteriors) for the discrete classes too, plus the classifier over classes.

Conclusions of v1 (these change the plan):
1. Within the modal3 skeleton, continuous within-mode disorder is required, calibrated and specific:
   - opening-rate disorder sigma_OC ≈ 0.4, slow (tau_w of minutes);
   - closing-rate disorder small or absent.
2. The tau_w prior truncates the posterior: extend it (log-U[1, 3000] s) before interpreting tau_w.
3. Neither model reproduces both ends. modal4 (discrete) has the frozen-mode tails but not the 30–300 s CK/ACF
   memory; modal-OU has the latter but not the former. The next class is modal4-OU (the modal4 skeleton with OU
   disorder within modes). It asks whether the slow disorder survives once the frozen mode is present, or whether
   both are needed (ladder plus continuum).
4. A slow 1D OU modulating kOC with tau_w of minutes is close to "more slow modes". The question "continuous or
   discrete" moves to the top-level structure: discrete frozen mode plus continuous slow modulation, or a longer
   discrete ladder. modal5 (discrete) should be in the model set.

## Phase 1 v2 — four model classes (2026-10-12)

Speed (before the runs):
- `simulators.simulate_modal_ou_fast` draws its random numbers with numpy's Generator in the reference's order and
  layout and runs the frame loop in numba. It reproduces `simulate_modal_ou` trajectories exactly (test, 6 settings,
  modal3 and modal4); 0.096 s per 4 x 6 h measurement with disorder (v1 numba 0.19 s, numpy 12 s).
- `summaries.ck_divergence` now builds its propagators from prefix sums and lagged dot products. The output is
  bit-identical to the masked-mean version (test; the stored v1 summary vectors are reproduced exactly), and
  `summary_vector` drops from 0.31 to 0.14 s.
- On 4 shared vCPUs with 4 workers: 0.29 s per simulation per worker; 1.5–2.1 h per set of 10^5.

modal5 ML fit (`scripts/fit_modal5.py`; `results/modal5.json`, `modal5_seeds.json`):
- Five L-BFGS runs, 1151 s: the default initialisation and modal4 with each mode split in two.
- Best log L = −98,104.1 (default and split-3 seeds agree to 0.06), vs modal4 −98,133.7: +29.6 for 5 parameters.
  BIC 196,431.5 vs 196,439.9, so BIC prefers modal5 by 8.4.
- modal5 keeps the frozen, slow-symmetric and open-biased modes of modal4 and splits the fast mode into two nearly
  identical modes:
  - <open> 1.5 s, <closed> 0.8 s, occupancy 0.20;
  - <open> 1.1 s, <closed> 0.8 s, occupancy 0.08.
  The two exchange through the closed state on 27–50 s, and through the open state only on hours.
  This is heterogeneity of kOC within the fast mode: the discrete counterpart of the kOC disorder.

Simulation sets (`scripts/run_phase1_v2.sh`; 10^5 each; batch 0 held out; own prior seed per set; fast simulator):

| set | class | parameters | prior |
|---|---|---|---|
| `results/phase1_m3ou_v2` | modal3-OU | 16 | modal3 rates x/÷2; sigma_OC, sigma_CO U[0,1]; tau_w log-U[1, 3000] s; rho U[-1,1] |
| `results/phase1_m4ou` | modal4-OU | 21 | modal4 rates x/÷2; OU priors as above |
| `results/phase1_m4` | modal4 (discrete) | 17 | modal4 rates x/÷2 |
| `results/phase1_m5` | modal5 (discrete) | 22 | modal5 rates x/÷2 |

NPE per class: same architecture and training as v1 (NSF, 5 transforms, 128 hidden; MLP embedding; 99,000 training
simulations; 50–68 min each). Calibration on the 1000 held-out simulations, 1000 posterior samples each:

| class | TARP ATC | coverage at 80 % / 95 % (range over parameters) | rank KS p < 0.01 | largest median offset / sd of standardised error |
|---|---|---|---|---|
| modal3-OU | +0.001 | 0.75–0.82 / 0.92–0.96 | 8 of 16 | 0.18 sd (kOC0) / 0.96–1.22 (tau_w widest) |
| modal4-OU | +0.007 | 0.76–0.82 / 0.92–0.96 | 5 of 21 | 0.15 sd / 0.97–1.12 |
| modal4 | +0.004 | 0.75–0.81 / 0.93–0.96 | 6 of 17 | 0.50 sd (kCO of the fast mode) / 0.98–1.40 |
| modal5 | −0.002 | 0.76–0.82 / 0.92–0.96 | 2 of 22 | 0.13 sd / 0.97–1.88 (kOC of fast mode A) |

- Joint calibration is good in every class and the coverage of central intervals is close to nominal.
- The rank tests detect small systematic offsets, largest for the sharpest posteriors (fast-mode rates, contraction
  0.99), where the flow over-concentrates by up to a factor 1.4–1.9 in sd.
- None of this affects the disorder parameters at the level reported below: their offsets are <= 0.09 sd, except
  log10 tau_w in modal3-OU (0.17 sd, sd of the standardised error 1.22).
- An NPE ensemble is not warranted for the conclusions here. It is for a final, paper-grade posterior (or sequential
  NPE around the data).

Posteriors for the data (median [90 %]):
- modal3-OU, tau_w prior extended to 3000 s:
  - sigma_OC = 0.38 [0.31, 0.47]; sigma_CO = 0.07 [0.01, 0.18] (P(sigma_CO < 0.1) = 0.69);
  - tau_w = 255 s [42, 1900 s]: no longer at the bound, but only bounded below;
  - rho = 0.83 [−0.24, 0.99].
  The v1 conclusions hold. Leakage outside the prior for the data rose to 19 % (rejected), as the posterior
  approaches sigma_CO = 0.
- modal4-OU:
  - sigma_OC = 0.23 [0.08, 0.32]; sigma_CO = 0.10 [0.01, 0.21]; P(sigma_OC < 0.1) = 0.07; P(both < 0.1) = 0.02;
  - tau_w = 171 s [15, 2160 s]; rho = −0.23 [−0.92, 0.73].
  With the frozen mode present, the kOC disorder shrinks from 0.38 to 0.23 but does not vanish.
- Discrete modal4 and modal5: rates close to the ML fits, with the tension below.

Summary-based posterior vs ML fit: the position of the full-likelihood ML parameters within the summary-based NPE
posterior of the same class (Mahalanobis^2 in the log-rates):
- modal4-OU: 7 over 17 dimensions (sigma = 0 margin);
- modal5: 13 over 22;
- modal3-OU: 16 over 12 (kOC of the open-biased mode at +2.9 sd);
- modal4: 36 over 17 (p ≈ 0.005 if the NPE were exact). kOC of three modes and kCO of the fast mode sit 2.1–3.6 sd
  from the posterior means. To match the summaries (the CK/ACF tail above all), modal4 needs closing rates that the
  dwell-sequence likelihood does not prefer. Read as a moderate tension: the NPE over-concentrates by up to 1.4x in
  sd for some of these rates, which inflates the number.

Symmetric posterior-predictive checks (`scripts/phase1_ppc_classes.py --name v2`;
`results/phase1_classes/ppc_v2.{json,npz}`, `figures/phase1_ppc_v2.png`). 200 measurements per class from its own
posterior; plug-in ML references for modal3 and modal4. Expectation for a correct model: ~6 of 113 statistics
outside the 95 % band.

| model | outside 95 % | \|z\| > 3 | sum z^2 | worst |
|---|---|---|---|---|
| modal3-OU (posterior) | 17 | 8 | 206 | 99.9 % open/closed quantiles (z 3.9, 3.3), CK at 70–140 s (3.0–3.4) |
| modal4-OU (posterior) | 5 | 0 | 96 | CK at 70 s (2.6) |
| modal4 (posterior) | 4 | 0 | 77 | — (all \|z\| < 2.5) |
| modal5 (posterior) | 3 | 0 | 81 | CK at 100 s (2.5) |
| modal3 (ML plug-in) | 40 | 31 | 1023 | open tail, CK, ACF |
| modal4 (ML plug-in) | 18 | 14 | 494 | CK 50–200 s (z up to 8), ACF 60–200 s |

- The large v1 gap between modal4 and modal-OU was mostly the plug-in: with parameter uncertainty and rates fitted to
  the summaries, discrete modal4 and modal5 pass the checks as well as modal4-OU.
- modal3-OU fails both tails: it lacks the frozen mode. Its posterior predictive also overshoots the fast-episode Fano
  factor (6.1 vs 5.1 in the data).
- Common residual: no class reproduces the shape of the CK divergence from q0 = open at 30–300 s. The data have a
  plateau (4.5e-4 nats from 33 to 68 s), then a drop by a factor ~30 between 100 and 300 s; every class predicts a
  smooth decline and leaves the data at the upper edge of its band (z ≈ 2.5 at 70–140 s). This points to a hidden
  process with a fairly well-defined time scale of ~2–3 min that none of the models contain in this form.

Model posterior (`scripts/phase1_classify.py --name v2`; `results/phase1_classes/classify_v2.{json,pt}`,
`figures/phase1_classify_v2.png`; 65 min):
- Setup: 4 x 99,000 training simulations (equal numbers, i.e. equal prior model probabilities); 5 % validation;
  the 4 x 1000 held-out batches as test set. Ensemble of 5 MLPs on robust-scaled, arcsinh-transformed summaries
  (121 inputs); softmax temperature 0.93 fitted on the validation set.
- Held-out: accuracy 0.80, log loss 0.49, ECE 0.009; reliability on the diagonal for every class.
- Confusion (rows true, columns predicted):
  - modal3-OU 0.90 / 0.06 / 0.02 / 0.02;
  - modal4-OU 0.13 / 0.69 / 0.09 / 0.08;
  - modal4 0.01 / 0.01 / 0.90 / 0.07;
  - modal5 0.03 / 0.04 / 0.23 / 0.69.
  modal4 vs modal5 and modal4-OU vs the rest are the hard pairs.
- Data: p(modal5) = 0.64 (members 0.55–0.79), p(modal4) = 0.31 (0.17–0.43), p(modal4-OU) = 0.05 (0.03–0.09),
  p(modal3-OU) = 7e-6.
  - Discrete (modal4 + modal5) vs continuous within-mode disorder (modal4-OU): 19:1 (11–36 over members).
  - modal5 vs modal4: 2:1, not decisive.
- The data are typical for the simulations: their nearest-neighbour distance to the training set is exceeded by
  55 % of the held-out simulations.

Conclusions of v2 (these change the plan; 1-2 superseded by the exact-likelihood programme below):
1. The v1 preference for continuous within-mode disorder came from comparing against modal3 and against ML plug-ins.
   Within the set {modal3-OU, modal4-OU, modal4, modal5}, inferred with the same machinery and checked
   symmetrically, the summaries prefer discrete modal gating:
   - modal5 (a fast mode split in kOC) or modal4, over modal4-OU by ~20:1;
   - modal3-OU is excluded (it lacks the frozen mode).
   The over-dispersion in fast episodes that motivated the continuum (Phase 0: Fano 5.1 vs 4.2–4.4 for the ML
   plug-ins) is reproduced once the rates are inferred: posterior-predictive medians 5.09 (modal5) and 4.94 (modal4).
2. The answer to "fixed or spread rates within modes" is therefore "fixed, with one more fast mode", at moderate
   strength. The continuum is not excluded: modal4-OU passes the predictive checks, keeps sigma_OC = 0.23
   [0.08, 0.32], and loses on Occam (4 extra parameters with broad priors). The odds depend on those priors.
3. modal4 needs closing rates in tension with its own ML fit to match the summaries (Mahalanobis^2 36 over 17
   dimensions); modal5 does not (13 over 22). This favours modal5 independently of the classifier.
4. Residual common to all classes: the shape of the CK divergence from q0 = open at 30–300 s (plateau, then a sharp
   drop). It points to a ~2–3 min process that none of the current classes contains.
Caveats:
- Model posteriors from non-sufficient summaries (Robert et al. 2011): they are posteriors given the summaries.
- All class priors are centred on ML fits to the same data (empirical Bayes, symmetric across classes).
- The NPEs over-concentrate slightly for the sharpest rates.
- The model set is not exhaustive (no GLE-gated class; no modal5-OU; no hour-scale tier).

## Exact-likelihood programme (2026-10-12 to 10-15): continuous kOC modulation on a cyclic mode network

Why. The Phase 1 v2 summary-based model posterior (discrete over continuous, 19:1) turned out to be prior-dominated
(see SBI reconciliation below). Where the likelihood is tractable it must decide. A 1D OU modulator discretised on a
grid is a finite AMM (Conclusions 3), so continuous within-mode disorder and discrete ladders can be compared on the
full dwell sequence.

### Methods

All in `src/hsp90memory/amm.py` and `scripts/`; 16 tests pass.

- **Exact likelihood:** `amm.loglik_Q(Q, classes, dwells)`.
  - Numba forward pass on the eigendecomposition of the within-class blocks of expm(Q dt), O(k^2) per dwell.
  - Equals `amm.loglik` to 1e-10 and is ~10x faster.
  - `amm.loglik_terms` gives per-dwell predictive terms.
- **Grid-OU models:** `amm.modal_ou(..., couple, topology)`.
  - One unit-variance OU coordinate z multiplies kOC_m by exp(sigma_OC z) in every mode.
  - Square-root approximation of the Fokker-Planck operator for U(z) = z^2/2, rescaled to exactly unit variance and
    slowest relaxation time tau_w. **Gauge: unit variance, tau_w = slowest relaxation.**
  - couple='eq': z relaxes in z^2/2 when open and in z^2/2 - sigma z when closed, so detailed balance holds exactly
    (flux asymmetry 1e-16). couple='OC': z is independent of the conformation (the SBI simulator's model; not balanced).
  - Validation, step 1: sigma -> 0 reproduces modal-K exactly.
  - Validation, step 2: an independent cumulative-power implementation and a frame-by-frame HMM agree (the difference
    is the first-frame term).
  - Validation, step 3: simulated exactly from its generator, the model is indistinguishable from
    `simulate_modal_ou_fast` (sigma_CO = 0) on all 113 summaries (max |z| 2.65), and the test has power (half sigma
    rejected, |z| 4.9).
  - Validation, step 4: grid convergence. log L changes by 0.6 from 7 to 15 points and by 0.25 from 15 to 31.
- **Grid policy:** ML on 15 points. Evidences sampled on 9 points and corrected to 21 points
  (log E[exp(ll21 - ll9)], about -0.75).
- **Mode networks:**
  - `amm.modal`: chain.
  - `amm.ring`: cycle, with detailed balance on every square and around the ring; db=False gives a driven cycle.
  - `amm.graph` / `named_graph`: any graph, with one balance constraint per independent cycle through a spanning tree.
    Named graphs: ringL1, ringL2, ringS2 (discrete kOC ladders off the fast mode), ringF1 (ring + a mode on the frozen
    state), full5 (all 10 pairs).
- **Bayesian inference:** `scripts/bayes_exact.py`.
  - Priors, data-independent: every free rate log-uniform on [1e-7, 1e2] s^-1; sigma_OC U[0, 2]; tau_w log-U[1, 1e4] s.
  - MCMC: 4 chains, adaptive Metropolis in logit space. Andrieu-Thoms scale and covariance adaptation from 2026-10-13
    (earlier runs: Haario).
  - Evidence: bridge sampling per chain; the chain spread is the error. Validated against brute-force integration
    for modal1 (both -110901.526).
- **ML search:** `fit_exact.py`, `fit_ladder.py`, `fit_graph.py`, with multi-start and basin-hopping.
  - **Pitfall:** L-BFGS from one family of starts stopped in local optima for modal6, modal7, the chain-5 + OU model
    and the first ringF1 fit.
  - MCMC found the better peaks each time. Reported ML values are the best of optimisation and MCMC
    (`scripts/_fits.best_fit_file`).
- **Checks:** `exact_cv.py` (leave-one-trace-out), `exact_attribution.py`, `exact_coverage.py` (per-data-set MCMC),
  `exact_coverage_laplace.py`, exact posterior-predictive checks via `phase1_ppc_classes.py label:exact:tag`.
- **Table and figure:** `exact_summary.py` writes `results/exact/summary_table.md`. `plot_model_schematic.py` writes
  `figures/model_schematic.{png,svg}`.

### Model comparison

n = 25,571 dwells; k = number of parameters; log Z with OU models corrected to the 21-point grid; held-out =
leave-one-trace-out CV.

| model | description | k | max log L | BIC | log Z | Δlog Z | held-out log L (CV) | note |
|---|---|---|---|---|---|---|---|---|
| modal3 | chain of 3 modes | 12 | -98360.5 | 196842.7 | -98416.3 | -407.8 | — |  |
| modal4 | chain of 4 modes | 17 | -98133.7 | 196439.9 | -98206.3 | -197.8 | -98177.2 (cv_v4) |  |
| modal5 | chain of 5 modes | 22 | -98092.3 | 196407.9 | -98181.8 | -173.4 | -98146.6 (cv_v5) |  |
| modal6 | chain of 6 modes | 27 | -98058.0 | 196390.1 | -98162.5 | -154.1 | -98115.8 (cv_v4) |  |
| modal7 | chain of 7 modes | 32 | -98013.2 | 196351.2 | -98157.2 | -148.7 | -98087.9 (cv_v5) | multimodal posterior: log Z of the dominant sampled mode |
| modal3_ou | chain 3 + OU on kOC | 14 | -98217.0 | 196576.1 | — | — | — |  |
| modal4_ou | chain 4 + OU on kOC (non-eq.) | 19 | -98056.6 | 196306.0 | -98134.0 | -125.6 | -98101.2 (cv_v4) |  |
| modal4_oueq | chain 4 + OU on kOC (eq.) | 19 | -98056.9 | 196306.6 | -98134.3 | -125.8 | -98100.3 (cv_v5) |  |
| modal4_ou2 | chain 4 + OU on kOC and kCO | 20 | -98056.4 | 196315.7 | — | — | — |  |
| modal5_oueqF | chain 5 + OU (eq.) | 24 | -97998.6 | 196240.9 | -98127.8 | -119.4 | -98052.5 (cv_v6) | log Z from an inferior peak (ML -98034.7 at the time) |
| ring4 | ring of 4 modes | 19 | -98050.0 | 196292.8 | -98125.6 | -117.1 | -98097.3 (cv_v6) |  |
| ring4_noDB | ring of 4, cycle not balanced | 20 | -98050.2 | 196303.4 | — | — | — |  |
| ringL1 | ring 4 + 1 discrete kOC ladder mode | 24 | -98007.7 | 196259.0 | -98100.6 | -92.1 | -98064.4 (cv_v5) |  |
| ringL2 | ring 4 + 2-step discrete kOC ladder | 29 | -97977.3 | 196249.0 | -98084.9 | -76.4 | -98037.9 (cv_v5) |  |
| ringS2 | ring 4 + 2 discrete kOC branches | 29 | -97976.3 | 196246.9 | -98084.3 | -75.8 | -98042.0 (cv_v6) |  |
| ringF1 | ring 4 + second frozen state | 24 | -97996.3 | 196236.2 | -98080.3 | -71.8 | -98040.0 (cv_v6) |  |
| ring4_oueq | ring 4 + OU on kOC (eq.) | 21 | -97971.2 | 196155.4 | -98052.1 | -43.7 | -98017.1 (cv_v6) |  |
| ring4_noDB_oueq | ring 4 + OU, cycle not balanced | 22 | -97971.1 | 196165.5 | — | — | -98022.8 (cv_v5) |  |
| ringF1_oueq | ring 4 + second frozen state + OU (eq.) | 26 | -97918.9 | 196101.6 | -98008.4 | 0.0 | -97964.1 (cv_v6) |  |
| full5_oueq | all 10 connections of 5 modes + OU (eq.) | 36 | -97909.2 | 196183.7 | — | — | -97957.0 (cv_v6) | likelihood ceiling for 5 modes |

### Final model: ringF1 + equilibrium OU on kOC

26 parameters; ML -97918.9 (15-point grid); log Z -98008.4 (chains -98007.4 to -98008.0).

Gating network. Five modes:

| mode | mean open | mean closed | share of time |
|---|---|---|---|
| slow | 15 s | 4.6 s | 40 % |
| open-biased | 19 s | 1.0 s | 18 % |
| fast | 1.4 s | 0.8 s | 29 % |
| long frozen state X | no transitions, mostly open | | 11 % (exits after ~4 min) |
| transient frozen hub | no transitions | | 2 % (lifetime <= ~3 s; exit rates bounded only from below) |

- X connects only to the hub. The hub connects to slow and fast.
- Cycle: hub - slow - open-biased - fast - hub.
- Switching events per day: hub-slow ~256, fast-hub ~53, hub-X ~45, slow-open-biased ~42, open-biased-fast ~17.

Within-mode modulation:

- **sigma_OC = 0.22 [0.19, 0.26]** (90 %): kOC varies by about +-25 % at 1 sd.
- **tau_w = 296 s [196, 448]** (R-hat 1.005, ESS > 450).
- The closed->open rate is not modulated: adding sigma_CO gained 0.2 on the chain skeleton.

Figure: `figures/model_schematic.png`.

### Findings

1. **Continuous modulation of kOC within the gating modes is decisive and does not depend on the skeleton.**
   - One OU coordinate (2 parameters) raises log L by 77-79 on every skeleton tested: chain of 4, ring of 4, and
     ring + second frozen state.
   - sigma_OC 0.21-0.22 and tau_w 267-296 s throughout.
   - Within the SBI v2 classes the NPE agrees: modal4-OU sigma_OC 0.23 [0.08, 0.32].
2. **Discrete ladders do not substitute.**
   - Same ring skeleton: two-step discrete kOC ladders (ringL2, ringS2; 29 parameters) reach -97976 to -97977 vs
     -97971 for ring + OU (21 parameters).
   - Evidence: 32 lower. Held-out log L: 21-25 lower.
   - Chain ladders up to 7 modes stay below chain-4 + OU in evidence.
3. **The gating modes form a cycle.**
   - In the chain-5 + OU optimum the two end modes were identical fast modes: the frozen side connects directly to
     the fast mode.
   - The ring of four gains 84 over the chain of four (discrete, +2 parameters) and 85 with OU.
4. **Detailed balance holds; no driven gating cycle in apo Hsp90.**
   - Dropping the balance constraint around the ring gains 0.036 log-units (likelihood-ratio statistic 0.07); out of
     sample it is worse by 5.8.
   - Equilibrium and non-equilibrium couplings of z to the conformation are indistinguishable (log Z -98134.3 vs
     -98134.0 on chain 4).
   - The equilibrium version is the physically appropriate model and is the one reported.
5. **The frozen state is two-staged.**
   - The second frozen state gains 52 log-units for 5 parameters on the ring + OU skeleton; evidence +43.7.
   - It is what reproduces the Chapman-Kolmogorov plateau at 30-140 s, Vollmar's shoulder and the Phase 1 v2 residual:
     ring + OU leaves z ~ 3 there.
   - Posterior-predictive check of the final model: 4/113 statistics outside the 95 % band, none with |z| > 2.2,
     sum z^2 = 40. For comparison: ring + OU 11/113 and 117; chain-5 + OU 7/113 and 84.
6. **The five-mode network is nearly complete.**
   - Connecting all ten pairs (full5 + OU, 36 parameters) gains 9.7 in log L; BIC is 82 worse.
   - It predicts held-out traces better by 7.0, with the same sign on all four traces (+2.5, +3.6, +0.4, +0.5).
   - Weak extra connections are possible; their evidence was not computed. The modulation parameters do not depend on
     it.
7. **Where the gain is (attribution, per-dwell predictive terms at the ML).**
   - Final model vs its discrete version: +77.5, almost all in open dwells (+73.1).
   - Per trace: -3.2, +19.6, +26.6, +34.5. Trace 1 is indifferent; traces 2-4 carry the evidence.
   - The modulation improves the prediction of each open dwell given its history: a slowly drifting closing rate.
8. **Cross-validation (cv_v6, held-out log L relative to the final model).**

   | model | held-out log L |
   |---|---|
   | ring + OU | -53.0 |
   | ringF1 (discrete) | -75.9 |
   | ringS2 | -78.0 |
   | chain-5 + OU | -88.5 |
   | ring4 | -133.3 |
   | full5 + OU | +7.0 |

### SBI reconciliation

- **The Phase 1 v2 classifier odds (discrete 19:1) were a prior-volume effect.**
  - The class priors were centred on each class's own ML fit, and the OU class carried 4 broad extra parameters.
  - Restricting the OU prior to the region holding 85 % of its posterior (sigma_OC <= 0.5, sigma_CO <= 0.2,
    tau_w 30-3000 s) gives 1.3:1 (`results/phase1_classes/prior_sensitivity_v2.json`).
- **The summaries cannot decide what the dwell-sequence likelihood decides by 77 log-units.**
  - Parameter posteriors agree where both exist.
  - Lesson: with non-sufficient summaries and empirical-Bayes class priors, a calibrated classifier can be confidently
    wrong about the class.
  - Use the exact likelihood wherever it exists. SBI remains the tool for classes without one (2D or GLE modulators,
    measurement physics).
- **The Phase 1 v2 conclusions 1-2 are superseded by this section.**

### Calibration

- **MCMC coverage** (`coverage_ringF1_oueq.json`, `coverage_modal4_oueq.json`).
  - Setup: parameter draws from the data posterior, exact simulation on a 15-point grid, re-inference with a
    7-point-grid MCMC.
  - Mean 95 % coverage over parameters: 0.88-0.89. sigma_OC 0.92, tau_w 0.92-0.95.
  - **Inconclusive:** the per-data-set chains had a median ESS of 24-32, so their intervals cover only part of the
    posterior. The U-shaped ranks are those of under-sampled chains.
- **Laplace-interval coverage** on the same 40 data sets (`coverage_laplace_ringF1_oueq.json`): running (2026-10-15). Real data first: the Laplace approximation matches the MCMC posterior for the identified parameters (sd ratio 1.00 for sigma_OC, 0.97 for tau_w, 0.9-1.02 for the well-determined rates; 0.26-0.63 only for the prior-dominated ones), so Laplace intervals are representative of the reported posterior there.
  Result over 40 data sets (10.6 h):
  - Identified parameters (18 of 26: Laplace/MCMC sd ratio 0.9-1.1 on the real data): coverage 0.53 / 0.81 / 0.89 /
    0.94 at 50 / 80 / 90 / 95 %, i.e. calibrated.
  - sigma_OC: 0.57 / 0.80 / 0.90 / 0.92; standardised error mean +0.09, sd 0.97.
  - log tau_w: 0.45 / 0.78 / 0.88 / 0.95; sd 1.00 but mean +0.47. The truth lies on average half a posterior sd
    above the estimate (tau_w ~12 % too short), plausibly from inferring on a 7-point grid data simulated on 15
    points; the real-data analysis used 9 and 15 points.
  - The 8 prior-dominated parameters (hub and frozen-state rates, some switching rates) undercover: 0.75 at 95 %.
    Expected, since Laplace is only 0.26-0.77 as wide as the MCMC posterior there; those parameters are reported as
    unidentified.
  - A fully sampled MCMC coverage (e.g. 12 data sets x 64,000 steps) would take ~8 h on this VM.

### Caveats

- Model posteriors are conditional on the model set. The skeleton search covered chains of 3-7 modes, rings, ring
  ladders, ringF1 and the full five-mode graph; six-mode networks were not explored.
- Evidence of a multimodal posterior (modal7) is a dominant-peak estimate.
- Unidentified rates mix slowly (R-hat up to 1.3):
  - the hub's exit rates and its O<->C rates;
  - the frozen-state O<->C rates;
  - some switching rates.
  Evidences are stable across chains (spread <= 0.6 for the final models).
- The CK plateau and the long-dwell tails rest on few long frozen episodes in four molecules.
- Trace 1 does not prefer the modulation.

## Next steps

Status 2026-10-15. Phase 1 is done (SBI v1/v2 plus the exact-likelihood programme above). Final model: ringF1 + an
equilibrium OU modulator on kOC.

1. **Calibration to publication standard.** The per-data-set MCMC coverage was under-sampled (ESS ~24).
   - Options: longer chains on a subset (e.g. 12 data sets x 64,000 steps, ~6-7 h on this VM), or a sampler with
     better mixing for the 26-dimensional posterior (parallel tempering, or gradients via automatic differentiation of
     the forward pass).
   - The Laplace-interval coverage (above) is the current calibration for the identified parameters.
2. **Robustness of the network.**
   - Evidence of full5 + OU: it is +7 out of sample, so weak extra connections are possible.
   - A six-mode search.
   - Per-trace fits of the final model (molecule-to-molecule variation of sigma_OC and tau_w; trace 1 does not prefer
     the modulation).
3. **Kinetic landscape (Phase 2, partly answered).**
   - What is settled: the slow coordinate couples to kOC only, log-linearly with sigma ~ 0.22, and it is the same
     coordinate in every mode.
   - Next: a spline coupling kOC(z) and a mode-dependent sigma.
   - A second coordinate (2D grid: ~100 x 10 states, still exact), or SBI for 2D/GLE modulators.
4. **ATP traces** (when available).
   - The driven-cycle test (ring without balance) and the equilibrium vs non-equilibrium coupling of z are ready to
     run.
   - The apo result (balanced cycle, eq. coupling) is the reference.
5. **Paper figures:** model schematic (done); the comparison table; posterior-predictive panels of the final model vs
   ring + OU (CK plateau); attribution by trace.
- Ask Hugel, in one email: the ATP traces (33,000 transitions); the smFRET traces under matching conditions; the 24 h
  traces and parallel-ruler ensemble; acquisition details (exposure time; whether the 10-Hz series is decimated from a
  faster acquisition — the absence of straddle frames suggests so). kHz bursts would make intra-well dynamics and the
  diffusive model meaningful.


## Files (repository `hsp90-memory/`)

- `src/hsp90memory/` — `io.py` (loading, idealisation), `amm.py` (AMM likelihood, fitting), `simulators.py` (AMM and
  nested modal-OU simulators, numpy reference and numba), `summaries.py` (summary statistics incl. CK divergence).
  `tests/` — 16 tests, `pytest`.
- `scripts/` — `run_fits.py`, `run_general_seeded.py`, `predictive_checks.py`, `phase0_simulate.py`, `phase0_refit.py`,
  `phase0_analyse.py`, `phase1_simulate.py` (checkpointed simulation driver), `phase1_smoke_report.py`,
  `phase1_npe.py` (NPE, SBC, TARP, posterior for the data), `phase1_ppc.py` (v1 posterior-predictive checks vs
  modal3/4), `phase1_ppc_classes.py` (symmetric checks across classes), `phase1_classify.py` (model posterior),
  `fit_modal5.py`, `run_phase1_v2.sh` (the v2 pipeline); exact-likelihood programme: `fit_exact.py`, `fit_ladder.py`,
  `fit_graph.py`, `bayes_exact.py`, `exact_cv.py`, `exact_attribution.py`, `exact_coverage.py`,
  `exact_coverage_laplace.py`, `exact_summary.py`, `plot_model_schematic.py`, `thin_samples.py`, `_fits.py`,
  `run_exact_pipeline{,2,3,4,5}.sh`;
  `scripts/diagnostics/` — the model-free diagnostics (`plot_trace.py`, `within_state.py`,
  `dwell_memory.py`, `sticking_test.py`).
- `results/*.json` — fitted generators Q, stationary distributions, log L, BIC; `results/phase0/` — summary vectors
  per setting (npz), refit log-likelihoods, `phase0_results.json`.
- `results/phase1/` (v1), `results/phase1_m3ou_v2/`, `results/phase1_m4ou/`, `results/phase1_m4/`,
  `results/phase1_m5/`: `settings.json`, logs, trained NPEs (`npe_<tag>.pt`), calibration (`npe_<tag>.json`,
  `npe_<tag>_sbc.npz`), posterior samples for the data (`npe_<tag>_posterior_data.npz`); simulation batches are not
  committed. `results/phase1_classes/`: symmetric posterior-predictive checks and the classifier.
  `results/modal5.json`, `modal5_seeds.json`. `results/exact/`: ML fits (`fit_*.json`), posteriors and evidences
  (`bayes_*.json`; full samples not committed, thinned copies in `results/exact/thin/`), CV (`cv_v*.json`),
  attribution, coverage, `summary_table.md`, pipeline logs.
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
