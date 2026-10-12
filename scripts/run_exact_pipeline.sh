#!/bin/bash
# Exact-likelihood programme: posteriors + evidences, posterior-predictive checks, coverage, attribution.
# Waits for the ML fits and the discrete Bayes chain started separately. Resumable: finished steps are skipped
# (each script refuses to overwrite its output). Run detached:
#   setsid nohup bash scripts/run_exact_pipeline.sh > results/exact/pipeline.log 2>&1 < /dev/null &
cd "$(dirname "$0")/.." && source .venv/bin/activate
step() { echo "=== $(date -u +%FT%TZ) START $1"; shift; "$@" && echo "=== $(date -u +%FT%TZ) OK" || echo "=== $(date -u +%FT%TZ) FAILED ($?)"; }
E=results/exact
for f in fit_modal6 fit_modal4_oueq fit_modal4_ou2 fit_modal5_ou bayes_modal5; do
  while [ ! -f $E/$f.json ]; do sleep 60; done
done
[ -f $E/bayes_modal6.json ]      || step "bayes modal6" python scripts/bayes_exact.py modal6 --steps 40000
[ -f $E/bayes_modal4_ou.json ]   || step "bayes modal4_ou" python scripts/bayes_exact.py modal4_ou --grid 9 --grid-correct 15 21 --steps 40000
[ -f $E/bayes_modal4_oueq.json ] || step "bayes modal4_oueq" python scripts/bayes_exact.py modal4_ou --couple eq --grid 9 --grid-correct 15 21 --steps 40000
[ -f results/phase1_classes/ppc_exact.json ] || step "ppc exact" python scripts/phase1_ppc_classes.py --name exact --n 200 \
  --classes modal4:exact:modal4 modal5:exact:modal5 modal6:exact:modal6 modal4-OU:exact:modal4_ou modal4-OUeq:exact:modal4_oueq
step "attribution" python scripts/exact_attribution.py --models modal4 modal5 modal6 modal4_ou modal4_oueq --ref modal5
[ -f $E/coverage_modal5.json ] || step "coverage modal5" python scripts/exact_coverage.py modal5 --n-sets 100 --steps 8000
BEST=$(python -c "import json; a=json.load(open('$E/bayes_modal4_ou.json')); b=json.load(open('$E/bayes_modal4_oueq.json')); g=lambda d: d['logZ_corrected'].get('21', d['logZ']); print('eq' if g(b) >= g(a) - 1 else 'OC')")
[ -f $E/coverage_modal4_ou$([ $BEST = eq ] && echo eq).json ] || step "coverage modal4_ou ($BEST)" python scripts/exact_coverage.py modal4_ou --couple $BEST --n-sets 60 --steps 8000 --fit-grid 7 --sim-grid 15
[ -f $E/bayes_modal5_ou.json ] || step "bayes modal5_ou" python scripts/bayes_exact.py modal5_ou --grid 9 --grid-correct 15 --steps 40000
echo "=== $(date -u +%FT%TZ) PIPELINE DONE"
