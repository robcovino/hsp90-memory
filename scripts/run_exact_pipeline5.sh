#!/bin/bash
# Final skeleton (ring + second frozen sub-state X, 'ringF1') with equilibrium OU on kOC: evidence, PPC, CV, attribution,
# coverage of the final model, thinned samples. Steps skip existing outputs.
#   setsid nohup bash scripts/run_exact_pipeline5.sh > results/exact/pipeline5.log 2>&1 < /dev/null &
cd "$(dirname "$0")/.." && source .venv/bin/activate
step() { echo "=== $(date -u +%FT%TZ) START $1"; shift; "$@" && echo "=== $(date -u +%FT%TZ) OK" || echo "=== $(date -u +%FT%TZ) FAILED"; }
E=results/exact
for f in fit_ringF1_oueq_global fit_full5_oueq; do while [ ! -f $E/$f.json ]; do sleep 60; done; done
BESTF=$(python -c "import sys; sys.path.insert(0,'scripts'); from _fits import best_fit_file; print(best_fit_file('ringF1_oueq'))")
echo "ringF1_oueq start: $BESTF"
[ -f $E/bayes_ringF1.json ] || step "bayes ringF1 (discrete)" python scripts/bayes_exact.py ringF1 --start $E/fit_ringF1.json --steps 80000
[ -f $E/bayes_ringF1_oueq.json ] || step "bayes ringF1_oueq" python scripts/bayes_exact.py ringF1_ou --couple eq --start $BESTF \
  --grid 9 --grid-correct 15 21 --steps 60000
[ -f results/phase1_classes/ppc_ringF1.json ] || step "ppc ringF1" python scripts/phase1_ppc_classes.py --name ringF1 --n 200 \
  --classes ringF1+OUeq:exact:ringF1_oueq ring4+OUeq:exact:ring4_oueq ringF1:exact:ringF1 chain5+OUeq:exact:modal5_oueqF
[ -f $E/cv_v6.json ] || step "cv v6" python scripts/exact_cv.py --out $E/cv_v6.json --models ring4 ringS2 ringF1 ring4_oueq ringF1_oueq \
  full5_oueq modal5_oueqF
[ -f $E/attribution_v6.json ] || step "attribution v6" python scripts/exact_attribution.py --out $E/attribution_v6.json --ref ringF1 \
  --models ringF1 ring4_oueq ringF1_oueq
FINAL=$(python -c "import json; a=json.load(open('$E/bayes_ring4_oueq.json')); b=json.load(open('$E/bayes_ringF1_oueq.json')); g=lambda d: d['logZ_corrected']['21']; print('ringF1' if g(b) > g(a) else 'ring4')")
echo "final model skeleton: $FINAL"
[ -f $E/coverage_${FINAL}_oueq.json ] || step "coverage ${FINAL}_oueq" python scripts/exact_coverage.py ${FINAL}_ou --couple eq \
  --n-sets 40 --steps 16000 --fit-grid 7 --sim-grid 15
step "thin samples" python scripts/thin_samples.py
echo "=== $(date -u +%FT%TZ) PIPELINE DONE"
