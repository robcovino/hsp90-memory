#!/bin/bash
# Final consistent runs with the Andrieu-Thoms adaptive sampler (suffix _v3), the combined model (continuous kOC
# modulation + two-state frozen mode, 'F'), predictive checks, CV, attribution, coverage of the final model.
#   setsid nohup bash scripts/run_exact_pipeline3.sh > results/exact/pipeline3.log 2>&1 < /dev/null &
cd "$(dirname "$0")/.." && source .venv/bin/activate
step() { echo "=== $(date -u +%FT%TZ) START $1"; shift; "$@" && echo "=== $(date -u +%FT%TZ) OK" || echo "=== $(date -u +%FT%TZ) FAILED"; }
E=results/exact
for f in fit_modal5_ouF fit_modal5_oueqF; do while [ ! -f $E/$f.json ]; do sleep 60; done; done
for K in 4 5 6 7; do
  [ -f $E/bayes_modal${K}_v3.json ] || step "bayes modal$K v3" python scripts/bayes_exact.py modal$K --suffix _v3 --steps 80000
done
[ -f $E/bayes_modal5_oueqF.json ] || step "bayes modal5_oueqF" python scripts/bayes_exact.py modal5_ou --couple eq --suffix F \
  --start $E/fit_modal5_oueqF.json --grid 9 --grid-correct 15 21 --steps 60000
while [ ! -f $E/coverage_modal4_oueq.json ]; do sleep 60; done
[ -f $E/bayes_modal4_oueq_v3.json ] || step "bayes modal4_oueq v3" python scripts/bayes_exact.py modal4_ou --couple eq --suffix _v3 \
  --grid 9 --grid-correct 15 21 --steps 60000
[ -f results/phase1_classes/ppc_exact_v3.json ] || step "ppc exact v3" python scripts/phase1_ppc_classes.py --name exact_v3 --n 200 \
  --classes modal4:exact:modal4_v3 modal5:exact:modal5_v3 modal6:exact:modal6_v3 modal7:exact:modal7_v3 \
            modal4-OUeq:exact:modal4_oueq_v3 modal5-OUeq-F:exact:modal5_oueqF
[ -f $E/cv_v3.json ] || step "cross-validation" bash -c "python scripts/exact_cv.py --models modal4 modal5 modal6 modal7 modal4_ou modal4_oueq modal5_ouF modal5_oueqF && mv $E/cv.json $E/cv_v3.json"
step "attribution" python scripts/exact_attribution.py --models modal4 modal5 modal6 modal7 modal4_oueq modal5_oueqF --ref modal5
BEST=$(python -c "import json; a=json.load(open('$E/bayes_modal4_oueq_v3.json')); b=json.load(open('$E/bayes_modal5_oueqF.json')); g=lambda d: d['logZ_corrected'].get('21', d['logZ']); print('F' if g(b) > g(a) else 'none')")
[ $BEST = F ] && { [ -f $E/coverage_modal5_oueqF.json ] || step "coverage modal5_oueqF" python scripts/exact_coverage.py modal5_ou --couple eq --suffix F --n-sets 48 --steps 8000 --fit-grid 7 --sim-grid 15; }
echo "=== $(date -u +%FT%TZ) PIPELINE DONE"
