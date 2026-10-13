#!/bin/bash
# Exact-likelihood programme, revised order (run_exact_pipeline.sh stopped after 'bayes modal4_oueq': the modal6 ML
# from modal5 splits was a local optimum). Steps skip existing outputs. Run detached:
#   setsid nohup bash scripts/run_exact_pipeline2.sh > results/exact/pipeline2.log 2>&1 < /dev/null &
cd "$(dirname "$0")/.." && source .venv/bin/activate
step() { echo "=== $(date -u +%FT%TZ) START $1"; shift; "$@" && echo "=== $(date -u +%FT%TZ) OK" || echo "=== $(date -u +%FT%TZ) FAILED"; }
E=results/exact
while [ ! -f $E/bayes_modal4_oueq.json ]; do sleep 60; done
step "global ladder fits" python scripts/fit_ladder.py 5 6 7 --hops 3
for K in 5 6 7; do
  [ -f $E/bayes_modal${K}_g.json ] || step "bayes modal$K (global start)" python scripts/bayes_exact.py modal$K --suffix _g --steps 60000
done
step "attribution" python scripts/exact_attribution.py --models modal4 modal5 modal6 modal7 modal4_ou modal4_oueq --ref modal5
[ -f $E/cv.json ] || step "cross-validation" python scripts/exact_cv.py --models modal4 modal5 modal6 modal7 modal4_ou modal4_oueq
[ -f results/phase1_classes/ppc_exact.json ] || step "ppc exact" python scripts/phase1_ppc_classes.py --name exact --n 200 \
  --classes modal4:exact:modal4 modal5:exact:modal5_g modal6:exact:modal6_g modal7:exact:modal7_g \
            modal4-OU:exact:modal4_ou modal4-OUeq:exact:modal4_oueq
BEST=$(python -c "import json; a=json.load(open('$E/bayes_modal4_ou.json')); b=json.load(open('$E/bayes_modal4_oueq.json')); g=lambda d: d['logZ_corrected'].get('21', d['logZ']); print('eq' if g(b) >= g(a) - 1 else 'OC')")
SUF=$([ $BEST = eq ] && echo eq)
[ -f $E/coverage_modal4_ou$SUF.json ] || step "coverage modal4_ou ($BEST)" python scripts/exact_coverage.py modal4_ou --couple $BEST --n-sets 60 --steps 8000 --fit-grid 7 --sim-grid 15
[ -f $E/coverage_modal5_g.json ] || step "coverage modal5" python scripts/exact_coverage.py modal5 --suffix _g --n-sets 100 --steps 8000
[ -f $E/bayes_modal4_ou${SUF}_long.json ] || step "long posterior modal4_ou ($BEST)" python scripts/bayes_exact.py modal4_ou --couple $BEST --suffix _long --grid 9 --grid-correct 15 21 --steps 100000
echo "=== $(date -u +%FT%TZ) PIPELINE DONE"
