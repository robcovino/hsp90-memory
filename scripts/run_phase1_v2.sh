#!/bin/bash
# Phase 1 v2: four model classes, 10^5 simulations each, NPE + SBC per class, symmetric posterior-predictive checks,
# model-posterior classifier. Every step is resumable (complete batches and trained estimators are reused); a failed
# step is logged and the pipeline continues. Run detached:
#   setsid nohup bash scripts/run_phase1_v2.sh > results/phase1_pipeline_v2.log 2>&1 < /dev/null &
cd "$(dirname "$0")/.." && source .venv/bin/activate
step() { echo "=== $(date -u +%FT%TZ) START $1"; shift; "$@" && echo "=== $(date -u +%FT%TZ) OK" || echo "=== $(date -u +%FT%TZ) FAILED ($?)"; }
sim() { python scripts/phase1_simulate.py --n-batches 100 --workers 4 --simulator simulate_modal_ou_fast "$@"; }

step "sim modal3-OU v2" sim --out results/phase1_m3ou_v2 --amm results/modal3.json --tau-max 3000 --seed 20261013
step "sim modal4-OU"    sim --out results/phase1_m4ou    --amm results/modal4.json --tau-max 3000 --seed 20261014
step "sim modal4"       sim --out results/phase1_m4      --amm results/modal4.json --discrete     --seed 20261015
while [ ! -f results/modal5.json ]; do sleep 60; done
step "sim modal5"       sim --out results/phase1_m5      --amm results/modal5.json --discrete     --seed 20261016

for c in m3ou_v2 m4ou m4 m5; do
  step "npe $c" python scripts/phase1_npe.py --tag $c --sims results/phase1_$c --outdir results/phase1_$c --test-batches 0
done
step "ppc classes" python scripts/phase1_ppc_classes.py --name v2 --n 200 --plugin modal3 modal4 --classes \
  modal3-OU:results/phase1_m3ou_v2:m3ou_v2 modal4-OU:results/phase1_m4ou:m4ou modal4:results/phase1_m4:m4 \
  modal5:results/phase1_m5:m5
step "classifier" python scripts/phase1_classify.py --name v2 --classes \
  modal3-OU:results/phase1_m3ou_v2 modal4-OU:results/phase1_m4ou modal4:results/phase1_m4 modal5:results/phase1_m5
echo "=== $(date -u +%FT%TZ) PIPELINE DONE"
