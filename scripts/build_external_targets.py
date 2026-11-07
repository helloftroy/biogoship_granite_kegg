#!/usr/bin/env python3
"""KO targets for an external test set on the KO set it shares with Bio-GO-SHIP.

CLR depends on the set of KOs it is taken over (the log values are centred on their
mean), so a model trained on Bio-GO-SHIP CLR can only be scored on another dataset
if both are on the same KO set and pseudocount rule. Here:

  shared KOs  = Bio-GO-SHIP target KOs (data/targets/ko_clr.csv) present in the
                external table, minus any with --min-prevalence below the cut in it.
  Bio-GO-SHIP = the existing CLR restricted to the shared KOs and re-centred. Equal to
                CLR taken over the shared KOs from the start (same pseudocount, half a
                count of the library total).
  external    = proportions of the library total (all KOs in that table), half a
                count as pseudocount, CLR over the shared KOs.

Single-copy-marker normalisation before CLR would change nothing but the pseudocount
(it divides each sample by one number, which CLR removes), so it is not applied; the
marker KOs are reported as a check of scale.

Only JGI (Flinkstrom et al. 2024, IMG annotation) has a ready KO table. Tara Pacific
and EMO BON need annotating first (see RESULTS.md R14).

Outputs in data/external/<set>/targets/:
  ko_clr_shared.csv           the external samples (sample_id as in the manifests)
  ko_clr_shared_biogoship.csv Bio-GO-SHIP on the same KOs (sample_key as data/targets)
  ko_compare.csv              per KO: prevalence and mean CLR in each, difference
  samples.csv                 per sample: total count, marker KO copies, kept or not

Usage:
    python scripts/build_external_targets.py --set jgi_flinkstrom
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship.tables import REPO_ROOT  # noqa: E402

EXTERNAL = REPO_ROOT / "data" / "external"
MARKERS = ["K06942", "K01889", "K01887", "K01875", "K01883", "K01869", "K01873", "K01409", "K03106", "K03110"]


def jgi_counts(sample_ids: list[str]) -> pd.DataFrame:
    """IMG KO gene counts of the given samples (JGI_<IMG genome id>) from the filtered h5ad."""
    import h5py

    path = EXTERNAL / "jgi_flinkstrom" / "adata_raw_filtered_v2.h5ad"
    with h5py.File(path, "r") as f:
        ids = ["JGI_" + x.decode() for x in f["obs"][f["obs"].attrs["_index"]][:]]
        kos = [x.decode().replace("KO:", "") for x in f["var"][f["var"].attrs["_index"]][:]]
        rows = np.flatnonzero(np.isin(ids, sample_ids))
        x = f["X"][np.sort(rows), :]
    return pd.DataFrame(x, index=np.array(ids)[np.sort(rows)], columns=kos)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--set", default="jgi_flinkstrom", choices=["jgi_flinkstrom"])
    parser.add_argument("--targets", type=Path, default=REPO_ROOT / "data" / "targets" / "ko_clr.csv")
    parser.add_argument("--target-samples", type=Path, default=REPO_ROOT / "data" / "targets" / "ko_samples.csv")
    parser.add_argument("--min-prevalence", type=float, default=0.9, help="Share of external samples a KO must be in.")
    parser.add_argument("--min-total", type=float, default=1e5, help="Smallest library total kept.")
    args = parser.parse_args()

    root = EXTERNAL / args.set
    manifest = pd.read_csv(root / "manifests" / "biogoship_samples.csv")
    counts = jgi_counts(manifest["sample_id"].tolist())
    ours = pd.read_csv(args.targets, index_col="sample_key")
    our_totals = pd.read_csv(args.target_samples, index_col="sample_key")["total"].reindex(ours.index)

    total = counts.sum(axis=1)
    info = pd.DataFrame({"total": total})
    for ko in MARKERS:
        info[ko] = counts[ko] if ko in counts else np.nan
    info["marker_median"] = info[MARKERS].median(axis=1)
    info["kept"] = total >= args.min_total
    counts = counts[info["kept"]]
    print(f"{args.set}: {len(info)} samples in the manifest with KO counts, {int(info['kept'].sum())} with >= {args.min_total:,.0f} counts")

    present = [ko for ko in ours.columns if ko in counts.columns]
    prevalence = (counts[present] > 0).mean()
    shared = prevalence[prevalence >= args.min_prevalence].index.tolist()
    print(f"Bio-GO-SHIP target KOs: {ours.shape[1]:,}; in the {args.set} table: {len(present):,}; "
          f"in >= {args.min_prevalence:.0%} of its samples: {len(shared):,}")
    missing_markers = [ko for ko in MARKERS if ko not in shared]
    if missing_markers:
        print(f"Marker KOs not in the shared set: {missing_markers}")

    ext_logs = np.log(counts[shared].div(counts.sum(axis=1), axis=0) + (0.5 / counts.sum(axis=1)).to_numpy()[:, None])
    ext_clr = ext_logs.sub(ext_logs.mean(axis=1), axis=0)
    our_clr = ours[shared].sub(ours[shared].mean(axis=1), axis=0)

    out = root / "targets"
    out.mkdir(parents=True, exist_ok=True)
    ext_clr.round(5).rename_axis("sample_id").to_csv(out / "ko_clr_shared.csv")
    our_clr.round(5).rename_axis("sample_key").to_csv(out / "ko_clr_shared_biogoship.csv")
    info.rename_axis("sample_id").join(manifest.set_index("sample_id")[["section", "site_id"]]).to_csv(out / "samples.csv")

    compare = pd.DataFrame({
        "prevalence_external": prevalence.reindex(shared), "mean_clr_external": ext_clr.mean(),
        "mean_clr_biogoship": our_clr.mean(), "sd_clr_external": ext_clr.std(), "sd_clr_biogoship": our_clr.std(),
    })  # fmt: skip
    compare["difference"] = compare["mean_clr_external"] - compare["mean_clr_biogoship"]
    compare = compare.sort_values("difference")
    compare.rename_axis("ko").round(4).to_csv(out / "ko_compare.csv")

    r = np.corrcoef(compare["mean_clr_external"], compare["mean_clr_biogoship"])[0, 1]
    print(f"\nMean CLR per KO, {args.set} vs Bio-GO-SHIP: r = {r:.3f}; difference median {compare['difference'].median():+.2f}, "
          f"|difference| > 2 (about 7-fold): {int((compare['difference'].abs() > 2).sum())} KOs")
    print("SD of CLR within dataset (median over KOs): "
          f"{args.set} {compare['sd_clr_external'].median():.2f}, Bio-GO-SHIP {compare['sd_clr_biogoship'].median():.2f}")
    print("Marker KOs, mean CLR (external / Bio-GO-SHIP):")
    print(compare.reindex([k for k in MARKERS if k in shared])[["mean_clr_external", "mean_clr_biogoship", "difference"]].round(2).to_string())
    print("\nLargest differences:")
    print(pd.concat([compare.head(6), compare.tail(6)])[["mean_clr_external", "mean_clr_biogoship", "difference"]].round(2).to_string())
    print(f"\nPer-sample totals: median {info.loc[info['kept'], 'total'].median():,.0f} genes with a KO "
          f"(Bio-GO-SHIP {our_totals.median():,.0f}); median marker KO count {info.loc[info['kept'], 'marker_median'].median():,.0f}")
    print(f"Wrote {out}/ko_clr_shared.csv, ko_clr_shared_biogoship.csv, ko_compare.csv, samples.csv")


if __name__ == "__main__":
    main()
