# hsp90-memory — working conventions

Project: memory in the conformational dynamics of Hsp90 measured by plasmon rulers (Ye et al. Nano Lett. 2018; Vollmar
et al. J. Phys. A 2024; data Zenodo 12697153). We ask what generates the memory, using aggregated Markov models (AMMs)
with exact likelihoods where they apply and simulation-based inference (SBI) where they do not. Background and the
state of the analysis are in `NOTES.md`; read its "Project summary (wrap-up)" and "Open problems" sections first.

## What is established (do not re-derive)

- **Data.** Traces are snapshot two-level signals at 10 Hz; thresholding is error-free (SNR ~ 11). All information is
  in the dwell sequence (25,571 dwells). No barrier/transition-path information; within-state fluctuations are white.
- **Final model (state 2026-10-15).** Discrete gating modes on a cycle plus a slow continuous modulation of the
  closing rate. In code: `amm.modal_ou(0, n_grid=15, couple="eq", topology="ringF1")`, with parameters in
  `results/exact/fit_ringF1_oueq_global.json`.
  - Network: modes slow, open-biased and fast, a transient frozen hub and a long frozen state. Cycle
    hub-slow-open-biased-fast-hub; the frozen state is attached to the hub.
  - Modulation: kOC_m exp(sigma z), with z a unit-variance equilibrium OU coordinate (gauge: unit variance,
    tau_w = slowest relaxation). sigma = 0.22 [0.19, 0.26], tau_w = 296 s [196, 448].
  - Fit quality: ML log L -97918.9 (15-point grid), log Z -98008.4. It reproduces all 113 summaries, including the
    2-3 min CK plateau.
- **What is settled.**
  - Fixed rates within discrete modes are ruled out: +77 log L for the 2 modulation parameters on every skeleton,
    evidence +72, held-out +76.
  - Detailed balance holds, including around the mode cycle.
  - Summary-statistic SBI could not decide this; its model posterior was prior-dominated. Use the exact likelihood
    where it exists.
- **Limits.**
  - The modulation is uneven between molecules (sigma 0.12-0.31; molecule 1 ~ none) and concentrated in time.
  - The local closing rate decorrelates faster at 5-10 min than any model predicts.
  - sigma and tau_w are effective parameters.

## Conventions

- Time unit seconds, frame interval `io.DT = 0.1`. Observed class 0 = open (low intensity), 1 = closed.
  kOC = open -> closed (closing rate).
- AMM states are ordered O_0, C_0, O_1, C_1, ... (pairs per mode). Grid-OU models: (z_0: all modal states),
  (z_1: ...), ...
- Likelihoods use P = expm(Q dt) (sub-frame excursions exact); use `amm.loglik_Q`. Never fit a discrete-time
  transition matrix directly.
- Rates are optimised in log space. Report log L, BIC and the number of dwells; always compare against the final model
  (ringF1_oueq) and its discrete version (ringF1).
- **Maximum-likelihood fits:** multi-start plus basin-hopping. Local optima are common for networks of five or more
  modes; check against MCMC samples. Use the best stored fit (`scripts/_fits.best_fit_file`).
- **Evidence:** `scripts/bayes_exact.py` (data-independent priors, adaptive MCMC, bridge sampling). OU models are
  sampled on 9 grid points and corrected to 21.
- Hidden continuous coordinates are identifiable only up to monotone reparameterisation: fix the gauge (unit-variance OU
  modulators; or U(y) = y^2/2 with free spline couplings) and say so in any result.
- **Any inference result must come with:**
  - (i) calibration: simulation-based calibration, or coverage from posterior draws (`exact_coverage*.py`);
  - (ii) posterior-predictive checks on the statistics in `summaries.py`: survival curves, lag correlations,
    local-rate histograms, CK profile (`phase1_ppc_classes.py`, with exact posteriors via `label:exact:tag`).
  No exceptions.
- Summary statistics are the ones in `summaries.summary_vector`; add new ones there, with a test.
- Data files are not committed. `data/README.md` says how to get them. Never normalise traces per-trace to [0, 1]
  (that discards calibration); use I/<I>.

## Working rules

- Run `pytest` before and after changes to `amm.py`, `simulators.py`, `summaries.py`.
- Long jobs: write results as JSON/NPZ into `results/` with the settings embedded; log to `results/*.log`.
- Do not delete or overwrite files in `results/` from earlier runs; version by name.
- Keep `NOTES.md` current: every finding that changes the plan goes in, with numbers.
- Do not submit cluster jobs or install packages system-wide without being asked.

## Next steps (as of 2026-10-15, project wrapped up)

See `NOTES.md`, "Open problems and next steps". In order:
1. the 5-10 min misfit of the local closing rate;
2. more molecules (hierarchical sigma);
3. fully sampled calibration;
4. network completeness (full five-mode evidence, six modes);
5. ATP traces (driven-cycle and non-equilibrium-coupling tests are ready);
6. questions to Hugel (simultaneous recording, acquisition details).
