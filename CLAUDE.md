# hsp90-memory — working conventions

Project: memory in the conformational dynamics of Hsp90 measured by plasmon rulers (Ye et al. Nano Lett. 2018; Vollmar
et al. J. Phys. A 2024; data Zenodo 12697153). We ask what generates the memory, using aggregated Markov models (AMMs)
with exact likelihoods where they apply and simulation-based inference (SBI) where they do not. Background and the
state of the analysis are in `NOTES.md` — read it first.

## What is established (do not re-derive)

- Traces are snapshot two-level signals at 10 Hz; thresholding is error-free (SNR ~ 11). All information is in the
  dwell sequence. No barrier/transition-path information; within-state fluctuations are white.
- Memory is dominated by hierarchical modal gating: a 3–4-mode AMM (`results/modal3.json`, `modal4.json`) reproduces
  dwell survival, dwell-lag correlations, local-rate bimodality and the Chapman–Kolmogorov profile to ~30 s. Residual:
  CK tail at 30–300 s under by ~2x. Detailed balance holds. The four molecules are one population.
- The open question for SBI is nested: within each gating mode, are the rates fixed (discrete) or spread (continuous
  disorder, parameters sigma, tau_w, rho in `simulators.simulate_modal_ou`)? Phase 0 (`scripts/phase0_*.py`) showed
  the data volume can tell (NOTES.md, "Phase 0"). Phase 1 (NPE over the nested model and model classes) is next.

## Conventions

- Time unit seconds, frame interval `io.DT = 0.1`. Observed class 0 = open (low intensity), 1 = closed.
- AMM states for modal models are ordered O_0, C_0, O_1, C_1, ... with mode 0 the slowest.
- Likelihoods use P = expm(Q dt) (sub-frame excursions exact). Never fit a discrete-time transition matrix directly.
- Rates are optimised in log space. Report log L, BIC and the number of dwells; always compare against modal3/modal4.
- Hidden continuous coordinates are identifiable only up to monotone reparameterisation: fix the gauge (unit-variance OU
  modulators; or U(y) = y^2/2 with free spline couplings) and say so in any result.
- Any inference result must come with (i) simulation-based calibration and (ii) posterior-predictive checks on the
  statistics in `summaries.py` — survival curves, lag correlations, local-rate histograms, CK profile. No exceptions.
- Summary statistics are the ones in `summaries.summary_vector`; add new ones there, with a test.
- Data files are not committed. `data/README.md` says how to get them. Never normalise traces per-trace to [0, 1]
  (that discards calibration); use I/<I>.

## Working rules

- Run `pytest` before and after changes to `amm.py`, `simulators.py`, `summaries.py`.
- Long jobs: write results as JSON/NPZ into `results/` with the settings embedded; log to `results/*.log`.
- Do not delete or overwrite files in `results/` from earlier runs; version by name.
- Keep `NOTES.md` current: every finding that changes the plan goes in, with numbers.
- Do not submit cluster jobs or install packages system-wide without being asked.

## Next steps (as of 2026-10-11)

1. Phase 0 is done (NOTES.md): within-mode disorder with sigma >= 0.25, tau_w >= 10 s is identifiable against
   3- and 4-mode discrete models at 4 x 6 h; the data lean towards sigma ~ 0.25-0.5, tau_w ~ 10-100 s. GO.
2. Phase 1 — NPE per model class (modal-K; modal-OU with sigma, tau_w, rho; optionally a GLE-gated model),
   ~1e5 simulations per class, classifier for model posterior, SBC, predictive checks on the real traces.
3. Phase 2 — kinetic landscape: spline couplings k_OC(y), k_CO(y) with the gauge above; joint distribution of the
   rate pair and its relaxation spectrum.
4. When the ATP traces arrive: time-reversal asymmetry tests, driven-network simulators, mode structure vs nucleotide.
