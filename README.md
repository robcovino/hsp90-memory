# hsp90-memory

Memory models for two-level single-molecule traces: the Hsp90 plasmon-ruler data (Zenodo 12697153), aggregated Markov
models with exact likelihoods, and the simulation-based-inference study of what generates the memory.

**Result.** The memory in the open-closed dynamics of apo Hsp90 comes from gating modes on a cycle, in detailed
balance. In addition, within every mode the closing rate drifts slowly: by about ±25 %, over ~5 min. Fixed rates
within discrete modes are ruled out by the exact dwell-sequence likelihood. Details, numbers and limits: `NOTES.md`,
section "Project summary (wrap-up, 2026-10-15)". Schematic: `figures/model_schematic.png`.

On a fresh Ubuntu VM (one command, then follow the printed instructions):

```
curl -fsSL https://raw.githubusercontent.com/robcovino/hsp90-memory/main/setup_vm.sh | bash
```

Manually: `pip install -e ".[dev,sbi]"` and `python -m pytest`.

Data: four 6-h traces at 10 Hz, see `data/README.md`. Analysis notes and current state: `NOTES.md`.
Conventions for working in this repository (also read by Claude Code): `CLAUDE.md`.

Layout
- `src/hsp90memory/io.py`         loading, idealisation, dwell extraction
- `src/hsp90memory/amm.py`        aggregated Markov models:
  - networks: `modal` (chain), `ring`, `graph` / `named_graph`, `split_mode`;
  - grid-OU modulation: `modal_ou`, equilibrium coupling `couple="eq"`;
  - exact likelihood: `loglik_Q`, per-dwell terms `loglik_terms`; fitting.
- `src/hsp90memory/simulators.py` AMM simulator (expm) and the nested modal-OU simulator (numpy reference,
  numba versions)
- `src/hsp90memory/summaries.py`  summary statistics incl. Chapman-Kolmogorov divergence (Vollmar et al. 2024)
- `scripts/`                      AMM fits and predictive checks; Phase 0 identifiability; Phase 1 SBI (simulation,
  NPE, classifier); exact-likelihood programme (`fit_exact.py`, `bayes_exact.py`, `exact_cv.py`,
  `exact_coverage*.py`, `exact_attribution.py`, `exact_summary.py`); audit (`audit_*.py`, `step0_*.py`,
  `fit_step0.py`, `fit_pertrace.py`); figures (`plot_model_schematic.py`, `plot_model_physchem.py`,
  `plot_paper_figure.py`)
- `results/`                      fitted generators (JSON); `phase0/`, `phase1*/` (SBI); `exact/` (ML fits, evidences,
  CV, coverage, summary table; thinned posterior samples in `exact/thin/`); `audit/`; `step0/`
- `figures/`                      model schematic, paper figure, physical-chemistry figure

Reproducing the main result: `python scripts/fit_exact.py`, `scripts/bayes_exact.py` and `scripts/exact_cv.py` as
recorded in `scripts/run_exact_pipeline*.sh`; `python scripts/exact_summary.py` rebuilds the comparison table from the
stored results.
