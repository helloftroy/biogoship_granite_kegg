#!/bin/bash
# R69 base-model fit for all hurdle KOs, 10 shards in parallel, resumable (rerun the same command after an interruption).
# Needs data/results/hurdle/c_nested.npz (scripts/hurdle_base.py --parts cpred). About 50 min on a 12-core Mac. Then:
#   OMP_NUM_THREADS=4 ~/miniforge3/bin/python3 scripts/hurdle_base.py --parts summary
cd "$(dirname "$0")/.." || exit 1
N=10
for ((i = 0; i < N; i++)); do
  OMP_NUM_THREADS=1 nohup caffeinate -i ~/miniforge3/bin/python3 scripts/hurdle_base.py --parts fit --shard $i $N \
    >> data/results/hurdle/base_shard$i.log 2>&1 &
done
wait
echo "hurdle base shards done"
