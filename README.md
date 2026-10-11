# hsp90-memory

Memory models for two-level single-molecule traces: the Hsp90 plasmon-ruler data (Zenodo 12697153), aggregated Markov
models with exact likelihoods, and the simulation-based-inference study of what generates the memory.

On a fresh Ubuntu VM (one command, then follow the printed instructions):

```
curl -fsSL https://raw.githubusercontent.com/robcovino/hsp90-memory/main/setup_vm.sh | bash
```

Manually: `pip install -e ".[dev,sbi]"` and `python -m pytest`.

Data: four 6-h traces at 10 Hz, see `data/README.md`. Analysis notes and current state: `NOTES.md`.
Conventions for working in this repository (also read by Claude Code): `CLAUDE.md`.

Layout
- `src/hsp90memory/io.py`         loading, idealisation, dwell extraction
- `src/hsp90memory/amm.py`        aggregated Markov models: modal(n), general(kO,kC); exact likelihood; fitting
- `src/hsp90memory/simulators.py` AMM simulator (expm) and the nested modal-OU simulator (within-mode disorder)
- `src/hsp90memory/summaries.py`  summary statistics incl. Chapman-Kolmogorov divergence (Vollmar et al. 2024)
- `scripts/`                      fits, predictive checks, Phase 0 identifiability study
- `results/`                      fitted generators (JSON), Phase 0 outputs
