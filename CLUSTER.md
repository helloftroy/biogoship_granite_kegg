# Running on the morrill cluster (HPC2, Mississippi State)

Operating notes for this project on the cluster, for the project owner and for any
assistant or agent picking the work up. General step-by-step instructions are in
[README.md](README.md); this file is what is specific to this cluster and this
account.

## Access

- Host `orion-login.hpc.msstate.edu` (SSH alias `hpc2`), user `hmp278`.
- Login is password plus Duo; there is no key-based or non-interactive path. **An
  agent cannot log in.** The working loop is: the agent writes and pushes code, the
  owner pulls on the cluster, submits jobs, and pastes job logs back.
- SLURM allocation: `191001-364393`. Pass it on the command line
  (`sbatch --account=191001-364393 ...`); it is not written into `#SBATCH` headers.
- Repo on the cluster: `/scratch/morrill/users/hmp278/biogoship_granite_kegg`.
  There is no sync with the local copy: `git pull` there after every push.

## Rules that have bitten before

- **`$HOME` is at quota. Nothing may be written there.** Everything goes under
  `/scratch/morrill/users/hmp278/`: data, model weights, Python environments,
  caches, temporary files. A write of even a few kilobytes to `$HOME` fails with
  `[Errno 122] Disk quota exceeded`. `cluster/env.sh` redirects `XDG_CACHE_HOME`,
  `MPLCONFIGDIR`, `HF_HOME`, `PIP_CACHE_DIR`, `CONDA_PKGS_DIRS` and `TMPDIR` to
  scratch, and sets `PYTHONNOUSERSITE=1`.
- **Do not create environments with `conda create`.** Even with `-p` pointing at
  scratch and the package cache on scratch, conda writes its own bookkeeping
  under `$HOME` (`~/.conda/environments.txt`, among others) and dies on the
  quota; this is what failed on 2026-10-13. Build environments as plain virtual
  environments instead: `<python on scratch> -m venv /scratch/.../conda_envs/<name>`
  followed by `pip install` with `PIP_CACHE_DIR` and `TMPDIR` on scratch.
  Activating an *existing* conda environment is fine.
- **Always `cd` into the repo before `sbatch`.** SLURM runs a spooled copy of the
  script, so `$0` does not point at the repo; the jobs use `$SLURM_SUBMIT_DIR`.
- **`conda` is not on `PATH`.** To use a conda environment, jobs source
  `/scratch/morrill/users/hmp278/miniforge3/etc/profile.d/conda.sh` and activate
  it by full path, then check that `sys.prefix` matches. Never activate by name:
  same-named environments exist under `$HOME`.
- **Every job requests `--time=48:00:00`.** There is no penalty, and it avoids jobs
  being killed during slow downloads.
- **Only the `service` partition has internet.** GPU and CPU nodes do not, so
  anything that downloads (data, pip packages, model weights) must run there first.
- Check free space before big downloads: `df -h /scratch/morrill/users/hmp278/`.

## Partitions

| Partition | Internet | Use |
|---|---|---|
| `service` | yes | downloads, environment setup, light processing (lands on a data-transfer node) |
| `gpu-a100` | no | GPU jobs (`--gres=gpu:1`) |
| `morrill` | no | CPU-only work |

Confirm with `sinfo` if a submission is rejected.

## Environments

| Environment | Path | Used by |
|---|---|---|
| `faire-agent` (pre-existing, shared) | `/scratch/morrill/users/hmp278/conda_envs/faire-agent` | product download and tiling (needs only numpy, h5py) |
| `biogoship` (this project) | `/scratch/morrill/users/hmp278/conda_envs/biogoship` | back-tracking (OceanParcels, Copernicus Marine toolbox) and Granite embeddings (torch, TerraTorch) |

`biogoship` is a plain virtual environment (pip, not conda), built by
`cluster/setup_env.sbatch` on top of the Python 3.11 in `faire-agent`;
`cluster/cache_granite.sbatch` adds TerraTorch to it. `cluster/env.sh` activates
either kind: a directory with `pyvenv.cfg` is sourced as a virtual environment,
anything else goes through `conda activate`. Any job can be pointed at another environment with
`--export=ALL,BGS_CONDA_ENV=/full/path`.

## Credentials

`cdse_credentials.sh` in the repo root (gitignored, `chmod 600`), copied from
`cdse_credentials.sh.example`. It holds two logins: Copernicus Data Space
(Sentinel-3) and Copernicus Marine (currents).

## Jobs, in order

All submitted from the repo root as `sbatch --account=191001-364393 cluster/<job>`.
Logs go to `logs/<job name>_<job id>.out` and `.err`.

| # | Job | Partition | What it does | Status (2026-10-13) |
|---|---|---|---|---|
| 1 | `download_products.sbatch` | service | all candidate OLCI scenes and SST files, about 300 GB, to `data/products/` | completed |
| 2 | `setup_env.sbatch` | service | creates the `biogoship` environment | first attempt (conda) failed on the `$HOME` quota; rewritten to use venv + pip |
| 3 | `backtrack.sbatch` | service | downloads currents, back-tracks every site, writes `data/manifests/advected_positions.csv` | not run |
| 4 | `download_tiles.sbatch` | service | cuts tiles at fixed and advected positions; checks cloud cover | not run |
| 5 | `cache_granite.sbatch` | service | installs TerraTorch, caches the Granite checkpoint, test-loads the encoder | not run |
| 6 | `embed_tiles.sbatch` | gpu-a100 | embeddings for every clear tile | not run |

Dependencies: 3 needs 2. 4 needs 1, and 3 for the advected tiles (it can run with
`POSITIONS=fixed` before 3, then again later to add the advected ones). 5 needs 2.
6 needs 4 and 5. `search_products.sbatch` is optional; its outputs are committed.

Every job resumes when resubmitted: finished downloads, finished (site, scene,
position) checks and existing currents files are skipped.

## Checking a run

```bash
squeue -u $USER                                   # is it running
tail -n 30 logs/bgs_tiles_<jobid>.out             # progress and the final summary
tail -n 30 logs/bgs_tiles_<jobid>.err             # errors
python scripts/download_products.py --dry-run     # how many archives are on disk (no login needed)
```

When reporting back, paste the tail of both `.out` and `.err`. The summaries at
the end of the tile and manifest steps are the results that go into
[METHODS_RESULTS.md](METHODS_RESULTS.md).

To bring result tables back to the local copy, commit the small CSVs from the
cluster (`git add data/manifests && git commit && git push`) or copy them with
`scp`. Tiles, products, currents and embeddings are gitignored and stay on scratch.

## Not yet verified on this cluster

- Authenticated Sentinel-3 download: verified (job 1 completed).
- Environment creation, the Copernicus Marine download, OceanParcels on
  real currents, reading real Sentinel-3 archives in the tile step, and loading
  Granite through TerraTorch in this environment have only been tested offline
  with synthetic data or carried over from the Florida project.
