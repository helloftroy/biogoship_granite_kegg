#!/usr/bin/env python3
"""Export the Granite benchmark bundle for a new foundation model (RESULTS.md R64).

Writes data/results/benchmark_export/ (copy the whole folder to the other HPC):
  samples.csv          every metagenome sample with KO targets (1,990; Protocol A excluded): position, time, depth,
                       cruise, sequencing platform, KO reads, and, for the 1,416 samples in the Granite benchmark,
                       the Sentinel-3 tile they were matched to (granule, tile date and centre) and their tile_key
  tiles.csv            the 1,215 benchmark tiles (rows of every target matrix, in this order): cruise, platform
                       group, samples averaged (k), mean sample position / time, Granite tile date and centre,
                       log10 KO reads (technical covariate), in-situ temperature, in_B_De (1,169 tiles have B / De)
  targets/             tile-level targets (samples sharing a tile averaged), raw (before technical removal), index
                       tile_key: clr_all.csv.gz (3,933 KOs), A, B, C, D, Dp, De .csv.gz
  targets_sample/      sample-level targets for every sample with KO counts (incl. those outside the benchmark),
                       index sample_key; used for within-tile variability and for extending to new tiles
  within_tile_pairs.csv  the 305 same-tile sample pairs (sample_a, sample_b, tile_key)
  granite/             tile_embeddings.npy (Granite, rows = tiles.csv), raw_satellite.npy (16 OLCI band means + SST),
                       predictions_<target>.npz (held-out main-model predictions on the adjusted targets),
                       scores.csv, partition.csv, ko_r2.csv, ko_ceiling.csv, module_predictability.csv,
                       category_predictability.csv, variance_budget.csv, cruise_decomposition.csv,
                       added_fields.csv
  kegg/                KO -> module links, module hierarchy, BRITE ko00001 table, KO names
  benchmark_eval.py    standalone evaluation (copy of scripts/benchmark_eval.py)
Usage: OMP_NUM_THREADS=4 ~/miniforge3/bin/python3 scripts/export_benchmark.py
"""
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
import main_model as mm  # noqa: E402
from water_mass_test import features  # noqa: E402

EX = REPO / "data" / "results" / "benchmark_export"
PLATFORM = {"OBI_2131": "HiSeq 4000", "OBI_2631": "NovaSeq 6000", "OBI_3310": "NovaSeq X"}


def main():
    for sub in ("targets", "targets_sample", "granite", "kegg"):
        (EX / sub).mkdir(parents=True, exist_ok=True)
    ks = pd.read_csv(REPO / "data" / "targets" / "ko_samples.csv").set_index("sample_key")
    fs = features()
    d = mm.setup("clr_all", fs, ks)
    unit, rows = d["unit"], d["rows"]
    tk = unit.index.to_numpy()

    # samples
    meta = pd.concat([pd.read_csv(REPO / "data/manifests/biogoship_samples.csv"),
                      pd.read_csv(REPO / "data/faire/manifests/biogoship_samples.csv")]).drop_duplicates("sample_id").set_index("sample_id")  # fmt: skip
    s = ks.rename(columns={"total": "ko_reads", "instrument": "instrument_obi"}).copy()
    s["platform"] = s["instrument_obi"].map(PLATFORM)
    for c in ("biosample", "collection_time_utc", "lat", "lon", "section", "expocode", "sample_method", "metagenome_runs"):
        s[c] = meta[c].reindex(s.index) if c in meta else np.nan
    r = rows.set_index("sample_id")
    s["in_granite_benchmark"] = s.index.isin(r.index)
    s["tile_key"] = r["tile_key"].reindex(s.index)
    for c in ("granule_id", "tile_date", "tile_lat", "tile_lon", "match_day_gap", "match_distance_km", "valid_fraction",
              "sst_valid_fraction", "temperature_c"):  # fmt: skip
        s[f"granite_{c}" if c.startswith(("tile", "match", "valid", "sst")) else c] = r[c].reindex(s.index)
    s.index.name = "sample_key"
    s.to_csv(EX / "samples.csv")

    # tiles
    g = rows.groupby("tile_key")
    t = pd.DataFrame(index=pd.Index(tk, name="tile_key"))
    t["cruise"] = d["cr"]
    t["platform_group"] = np.where(d["plat"], "NovaSeq 6000 group", "NovaSeq X group")
    t["k_samples"] = d["k"]
    t["sample_keys"] = g["sample_id"].agg(lambda x: ";".join(x)).reindex(tk).to_numpy()
    t["sample_time_mean"] = pd.to_datetime(g["sample_time"].agg(lambda x: pd.to_datetime(x, utc=True).mean())).reindex(tk).astype(str).to_numpy()
    t["sample_lat"] = g["sample_lat"].mean().reindex(tk).to_numpy()
    t["sample_lon"] = g["sample_lon"].mean().reindex(tk).to_numpy()
    for c in ("granule_id", "tile_date", "tile_lat", "tile_lon"):
        t[f"granite_{c}"] = unit[c].to_numpy()
    t["log10_ko_reads"] = d["logtot"]
    t["temperature_c"] = d["temp"][:, 0]
    dB = mm.setup("B", fs, ks)
    # tile_key is renumbered per target set by the loader (B / De have 1,169 tiles); map every set to the
    # all-KO keys through the embedding row
    to_main = pd.Series(tk, index=unit["embedding_row"].to_numpy())
    remap = lambda dd: to_main.loc[dd["unit"]["embedding_row"].to_numpy()].to_numpy()  # noqa: E731
    t["in_B_De"] = t.index.isin(remap(dB))
    t.to_csv(EX / "tiles.csv")

    # targets
    for name in mm.ORDER:
        dd = d if name == "clr_all" else (dB if name == "B" else mm.setup(name, fs, ks))
        keys = remap(dd)
        pd.DataFrame(dd["y"], index=pd.Index(keys, name="tile_key"), columns=dd["kos"]).to_csv(
            EX / "targets" / f"{name}.csv.gz", float_format="%.6g")  # fmt: skip
        src = mm.TARGETS[name]
        st = pd.read_csv(src, index_col="sample_key")
        st.to_csv(EX / "targets_sample" / f"{name}.csv.gz", float_format="%.6g")  # all samples (1,990 for clr_all)
        z = np.load(mm.OUT / f"predictions_{name}_granite_cv.npz", allow_pickle=True)
        zz = {k: z[k] for k in z.files}
        assert len(zz["tile"]) == len(keys)
        zz["tile"] = keys  # all-KO tile keys (tiles.csv)
        np.savez_compressed(EX / "granite" / f"predictions_{name}.npz", **zz)
    pairs = []
    for key, samples in rows.groupby("tile_key")["sample_id"].agg(list).items():
        pairs += [(samples[i], samples[j], key) for i in range(len(samples)) for j in range(i + 1, len(samples))]
    pd.DataFrame(pairs, columns=["sample_a", "sample_b", "tile_key"]).to_csv(EX / "within_tile_pairs.csv", index=False)

    # Granite reference
    np.save(EX / "granite" / "tile_embeddings.npy", d["emb"].astype(np.float32))
    np.save(EX / "granite" / "raw_satellite.npy", d["sat"].astype(np.float32))
    for f in ("scores.csv", "partition.csv", "ko_ceiling.csv", "module_predictability.csv", "category_predictability.csv",
              "variance_budget.csv", "cruise_decomposition.csv", "added_fields.csv", "c_per_cruise.csv", "c_per_cruise_r2.csv"):  # fmt: skip
        shutil.copy(mm.OUT / f, EX / "granite" / f)
    shutil.copy(mm.OUT / "ko_temperature.csv", EX / "granite" / "ko_r2.csv")
    for f in ("link_module_ko.tsv", "module_list.tsv", "ko00002.json"):
        shutil.copy(REPO / "data/reference/kegg_modules" / f, EX / "kegg" / f)
    for f in ("kegg_brite_ko00001.csv", "kegg_ko_list.tsv"):
        shutil.copy(REPO / "data/reference" / f, EX / "kegg" / f)
    shutil.copy(REPO / "scripts" / "benchmark_eval.py", EX / "benchmark_eval.py")
    print(f"samples {len(s):,} (benchmark {int(s['in_granite_benchmark'].sum()):,}); tiles {len(t):,}; pairs {len(pairs)}; "
          f"embedding {d['emb'].shape}; written to {EX}")  # fmt: skip


if __name__ == "__main__":
    main()
