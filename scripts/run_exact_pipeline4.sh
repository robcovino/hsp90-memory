#!/bin/bash
# Ring-skeleton programme: evidences of the discrete competitors, predictive checks, CV, attribution, coverage of the
# final model (ring4 + equilibrium OU on kOC), thinned samples. Steps skip existing outputs.
#   setsid nohup bash scripts/run_exact_pipeline4.sh > results/exact/pipeline4.log 2>&1 < /dev/null &
cd "$(dirname "$0")/.." && source .venv/bin/activate
step() { echo "=== $(date -u +%FT%TZ) START $1"; shift; "$@" && echo "=== $(date -u +%FT%TZ) OK" || echo "=== $(date -u +%FT%TZ) FAILED"; }
E=results/exact
for m in ring4 ringS2 ringL2 ringL1; do
  [ -f $E/bayes_$m.json ] || step "bayes $m" python scripts/bayes_exact.py $m --start $E/fit_$m.json --steps 80000
done
while [ ! -f $E/bayes_ring4_oueq.json ]; do sleep 60; done
[ -f results/phase1_classes/ppc_ring.json ] || step "ppc ring" python scripts/phase1_ppc_classes.py --name ring --n 200 \
  --classes ring4+OUeq:exact:ring4_oueq ring4:exact:ring4 ringS2:exact:ringS2 ringL2:exact:ringL2 \
            chain5+OUeq:exact:modal5_oueqF modal7:exact:modal7_v3
[ -f $E/cv_v5.json ] || step "cv v5" python scripts/exact_cv.py --out $E/cv_v5.json --models modal5 modal7 modal4_oueq modal5_oueqF \
  ring4 ringL1 ringL2 ringS2 ring4_oueq ring4_noDB_oueq
[ -f $E/attribution_v5.json ] || step "attribution v5" python scripts/exact_attribution.py --out $E/attribution_v5.json --ref ring4 \
  --models ring4 ringS2 ring4_oueq modal5_oueqF
[ -f $E/coverage_ring4_oueq.json ] || step "coverage ring4_oueq" python scripts/exact_coverage.py ring4_ou --couple eq --n-sets 48 \
  --steps 16000 --fit-grid 7 --sim-grid 15
step "thin samples" python scripts/thin_samples.py
echo "=== $(date -u +%FT%TZ) PIPELINE DONE"
