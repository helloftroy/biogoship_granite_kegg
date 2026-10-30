# Ptolemy (HPC2) first checks

For a planned satellite foundation model trained on 4 A100s at once, with 10+ TB
of NASA and Copernicus imagery. These checks answer two questions before any
setup: how much scratch can this account use, and can one job hold 4 GPUs?

```bash
# on a Ptolemy login node
mkdir -p /scratch/ptolemy/users/$USER && cd /scratch/ptolemy/users/$USER
git clone https://github.com/helloftroy/biogoship_granite_kegg.git
cd biogoship_granite_kegg && mkdir -p logs

bash cluster/ptolemy/check_cluster.sh 2>&1 | tee logs/ptolemy_check.txt     # storage, quotas, partitions, GPUs, accounts
sbatch --account=<account> --partition=<gpu partition> cluster/ptolemy/gpu4_test.sbatch
cat logs/gpu4_test_*.out                                                    # ends with "RESULT: 4 GPUs OK"
```

Take the account and GPU partition from the `sacctmgr` and `sinfo` sections of
the check output. `ddp_check.py` (run by the GPU test when `TORCH_ENV` is set)
confirms the 4 GPUs can work together and how fast they exchange data. It is for
later, once a PyTorch environment exists.

Still to ask HPC2 support if the check does not show it: the scratch quota, the
purge policy (files unused for N days deleted?) and whether 10+ TB is allowed
there or needs a project allocation.
