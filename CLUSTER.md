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

Version constraint to keep: OceanParcels 3 requires zarr 2, and xarray from 2026
onward only works with zarr 3. Unpinned, pip installs both and the Copernicus
Marine download fails with `open_group() got an unexpected keyword argument
'zarr_format'`. `cluster/env_packages.txt` therefore holds `xarray<2026` and
`zarr<3`. Re-running `setup_env.sbatch` applies changes to that file.

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
| 2 | `setup_env.sbatch` | service | creates the `biogoship` environment | built 2026-10-13 with venv + pip (the first, conda, attempt failed on the `$HOME` quota). Needs one re-run to pick up the xarray pin below |
| 3 | `backtrack.sbatch` | service | downloads currents, back-tracks every site, writes `data/manifests/advected_positions.csv` | done 2026-10-13 on real currents: all 3,727 pairs resolved, 7-day tracks complete |
| 4 | `download_tiles.sbatch` | service | cuts tiles at fixed and advected positions; checks cloud cover | fixed position run 2026-10-13 done 2026-10-13 for both positions within the original window: fixed 614 samples, advected 611, both 589 |
| 4b | `extend_search.sbatch` | service | for sites without a clear advected tile, finds scenes further out along the water's track; then rerun 1 and 4 | not run |
| 5 | `cache_granite.sbatch` | service | installs TerraTorch, caches the Granite checkpoint, test-loads the encoder | done 2026-10-13; encoder loads cleanly |
| 6 | `embed_tiles.sbatch` | gpu-a100 | embeddings for every clear tile | done 2026-10-20: 3,503 tiles in 101 s |
| 6b | `check_granite_inputs.sbatch` | morrill (CPU) | due diligence: tiles vs the model's config, and MAE reconstruction under different gap-filling conventions | done 2026-10-20 (job 80958): inputs match; keep current gap filling |
| 5c | `fill_missing_sst.sbatch` | service | fills tiles without SST from the nearest-in-time SLSTR overpass of the same spot; refreshes manifest and raw features; then re-run 6 | not run |
| N1 | `new_cruises.sbatch` | service | newer cruises: CCHDO + SeaBASS candidate events, Sentinel-3 scene search, download size (no download) | done 2026-10-22 for CCHDO (917 sites; 723 OLCI + 489 SST, ~610 GB); SeaBASS failed at the Earthdata login, fixed in code, to re-run |
| N2 | `download_products.sbatch` with `MANIFEST_DIR=data/new_cruises/manifests` | service | Sentinel-3 archives for the newer cruises | not run |
| 7 | (no job) `python scripts/pool_tile_features.py` | login node is fine (a minute of file reading) | raw-channel summary per tile, for the baseline; copy `data/embeddings/tile_raw_features.csv` back | done 2026-10-20 |

Where things stand (2026-10-20): jobs 1 to 6 done (advected tiles for 663 of
721 sites; embeddings for all 3,503 tiles). The embedding bundle was copied back
and the first CatBoost comparison run locally. Next: step 7 on the cluster for
the raw-channel baseline.

Note: 4b rewrites `data/manifests/site_products.csv`, `products.csv` and
`advected_positions.csv` on the cluster, so the cluster copies of those tracked
files now differ from GitHub. `git pull` still works as long as nobody changes
them upstream; do not `git checkout` or `git reset` them on the cluster. To
publish the final versions, commit them from the cluster or copy them back.

Harmless messages in `.err` files: "unauthenticated requests to the HF Hub"
(set `HF_TOKEN` to silence it), "`torch.jit.script` is deprecated", and "Can't
initialize NVML" on nodes without a GPU.

Original plan (2026-10-13): finish the fixed case (4 with `RECHECK="cloudy edge"`,
to apply the 20 km rule to combinations tested before it existed), then 3, then 4
for the advected position, then 4b -> 1 -> 4 to widen the window for sites still
missing an advected tile, then 5 and 6.

Dependencies: 3 needs 2. 4 needs 1, and 3 for the advected tiles (it can run with
`POSITIONS=fixed` before 3, then again later to add the advected ones). 5 needs 2.
4b needs 3 and an advected run of 4. 6 needs 4 and 5. `search_products.sbatch` is optional; its outputs are committed.

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

## Lessons from this project

- A log file that is appended to across code versions must tolerate its own
  history. `data/tile_log.csv` gained columns twice; rows appended under an older
  header were misread. It is now read by field count and rewritten in the current
  layout automatically (`biogoship/tables.py`), with the original kept once as
  `tile_log.before_layout_update.csv`.
- A job that skips work because "the file exists" must check the file is still
  what is needed. The currents files were kept from a run with a shorter time
  span; `download_currents.py` now checks coverage.

## Not yet verified on this cluster

- Authenticated Sentinel-3 download: verified (job 1 completed).
- Verified since: environment creation (venv), the Copernicus Marine download,
  OceanParcels on real currents, and reading real Sentinel-3 archives in the
  tile step.
- Still only tested offline or carried over from the Florida project: loading
  Granite through TerraTorch in this environment, and the wider scene search
  along tracks at full scale.
