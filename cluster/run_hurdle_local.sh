#!/bin/bash
# Hurdle fit for all KOs (RESULTS.md R65), 4 shards in parallel, resumable (rerun the same command after an interruption).
# Run from the repo root in your own terminal:  bash cluster/run_hurdle_local.sh
# caffeinate keeps the Mac awake while the shards run. Then:
#   OMP_NUM_THREADS=4 ~/miniforge3/bin/python3 scripts/hurdle_kos.py --parts summary figure
cd "$(dirname "$0")/.." || exit 1
for i in 0 1 2 3; do
  OMP_NUM_THREADS=3 nohup caffeinate -i ~/miniforge3/bin/python3 scripts/hurdle_kos.py --parts fit --shard $i 4 \
    >> data/results/hurdle/fit_shard$i.log 2>&1 &
done
wait
echo "hurdle shards done"
