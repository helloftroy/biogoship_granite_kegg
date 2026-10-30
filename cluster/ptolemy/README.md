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
sbatch --account=191001-366746 cluster/ptolemy/gpu4_test.sbatch
cat logs/gpu4_test_*.out                                                    # ends with "RESULT: 4 GPUs OK"
```

What the check found (2026-10-30):

- **Account** `191001-366746` (Morrill's is `191001-364393`), QOS `normal`: at most
  4 whole A100s, 128 CPUs and 1 TB memory per user, 2 days per job, 128 jobs
  running and 256 queued.
- **Partitions:** `gpu-a100` (4 nodes, 8 whole A100s each, 1 TB memory) is the one
  for training. `gpu-a100-mig2` and `gpu-a100-mig7` (default) split each GPU into
  40 GB or 10 GB slices; `development` has small slices for tests; `service`
  has no GPUs.
- **Scratch** `/scratch/ptolemy/users/hmp278` is GPFS, 300 TB, 248 TB free; no
  per-user quota was reported. `/home` is a separate 5 TB GPFS (4.2 TB free).
- **Modules:** python 3.9/3.10, miniconda3 4.10/24.3, cuda 11.7/12.1, cuda-toolkit
  12.9, gcc 12.2.
- **Internet** from the login node: GitHub and NASA Earthdata reachable. `ddp_check.py` (run by the GPU test when `TORCH_ENV` is set)
confirms the 4 GPUs can work together and how fast they exchange data. It is for
later, once a PyTorch environment exists.

Still to ask HPC2 support if the check does not show it: the scratch quota, the
purge policy (files unused for N days deleted?) and whether 10+ TB is allowed
there or needs a project allocation.
