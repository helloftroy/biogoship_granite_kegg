#!/usr/bin/env python3
"""Run every tile through the Granite encoder: one 512-number embedding per tile.

Reads data/manifests/biogoship_tile_manifest.csv and embeds each distinct tile
once, for both positions (`fixed` and `advected`) and for every clear overpass,
not just the closest. Only tiles that passed the clear-ocean check exist, so
cloudy and out-of-scene cases are already excluded; they are listed per site in
data/manifests/site_tile_status.csv.

Outputs in data/embeddings/:
  tile_embeddings.npy      float32 (n_tiles, 512)
  tile_embeddings.csv      one row per tile, in the same order: tile_id, site_id,
                           position, scene, time gap, valid fractions
  sample_embeddings.csv    one row per (sample, tile): the tile manifest plus
                           `embedding_row`, the row of that tile in the .npy
  granite_embeddings_bundle.zip   the three files above plus summary.json

Usage:
    python scripts/embed_tiles.py
    python scripts/embed_tiles.py --local-files-only     # no internet (GPU nodes)
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import zipfile
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship import granite  # noqa: E402
from biogoship.tables import MANIFEST_DIR, REPO_ROOT, from_repo, read_csv, write_csv  # noqa: E402

TILE_COLUMNS = [
    "site_id",
    "position",
    "granule_id",
    "tile_date",
    "tile_lat",
    "tile_lon",
    "match_hour_gap",
    "offset_km",
    "valid_fraction",
    "sst_valid_fraction",
    "is_best",
    "same_scene_as_best_fixed",
    "displacement_km",
    "tile_path",
]


def main() -> None:
    args = parse_args()
    manifest = read_csv(args.manifest_csv)
    if args.positions:
        manifest = [row for row in manifest if row["position"] in args.positions]
    tiles: dict[str, dict[str, str]] = {}
    for row in manifest:  # replicate samples share a tile; embed it once
        tiles.setdefault(row["tile_path"], row)
    present = {path: row for path, row in tiles.items() if from_repo(path).exists()}
    if len(present) < len(tiles):
        print(f"WARNING: {len(tiles) - len(present):,} tiles in the manifest are not on disk and are skipped.")
    if not present:
        raise SystemExit(f"No tiles to embed from {args.manifest_csv}.")
    paths = sorted(present)
    counts = {p: sum(1 for path in paths if present[path]["position"] == p) for p in sorted({r["position"] for r in present.values()})}
    print(f"Tiles to embed: {len(paths):,} (" + ", ".join(f"{n:,} {p}" for p, n in counts.items()) + ")")

    device = torch.device(args.device or ("cuda" if torch.cuda.is_available() else "cpu"))
    encoder, report = granite.build_encoder(args.model_name, args.local_files_only)
    print(report)
    print(f"Device: {device}", flush=True)
    encoder.eval().to(device)

    started = time.time()
    batches = []
    with torch.no_grad():
        for start in range(0, len(paths), args.batch_size):
            batch = granite.load_tiles([from_repo(path) for path in paths[start : start + args.batch_size]]).to(device)
            batches.append(granite.encode(encoder, batch).cpu().numpy().astype(np.float32))
            if (start // args.batch_size) % 20 == 0:
                print(f"  embedded {min(start + args.batch_size, len(paths)):,}/{len(paths):,}", flush=True)
    embeddings = np.concatenate(batches, axis=0)
    if not np.isfinite(embeddings).all():
        raise SystemExit("Non-finite values in the embeddings; inspect the input tiles.")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    row_of = {path: index for index, path in enumerate(paths)}
    tile_rows = [
        {"embedding_row": index, "tile_id": Path(path).stem, **{c: present[path].get(c, "") for c in TILE_COLUMNS}}
        for index, path in enumerate(paths)
    ]
    sample_rows = [{"embedding_row": row_of[row["tile_path"]], **row} for row in manifest if row["tile_path"] in row_of]
    summary = {
        "model_name": args.model_name,
        "n_tiles": int(embeddings.shape[0]),
        "embedding_dim": int(embeddings.shape[1]),
        "tiles_by_position": counts,
        "n_sample_rows": len(sample_rows),
        "pooling": "mean over encoder output tokens",
        "nan_fill_value": 0.0,
        "device": str(device),
        "seconds": round(time.time() - started, 1),
    }
    np.save(args.out_dir / "tile_embeddings.npy", embeddings)
    write_csv(args.out_dir / "tile_embeddings.csv", tile_rows)
    write_csv(args.out_dir / "sample_embeddings.csv", sample_rows)
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    bundle = args.out_dir / "granite_embeddings_bundle.zip"
    with zipfile.ZipFile(bundle, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name in ("tile_embeddings.npy", "tile_embeddings.csv", "sample_embeddings.csv", "summary.json"):
            archive.write(args.out_dir / name, arcname=name)
    print(f"Embeddings: {embeddings.shape[0]:,} tiles x {embeddings.shape[1]} dimensions in {summary['seconds']}s")
    print(f"Wrote {bundle}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest-csv", type=Path, default=MANIFEST_DIR / "biogoship_tile_manifest.csv")
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "data" / "embeddings")
    parser.add_argument("--positions", nargs="+", choices=["fixed", "advected"], help="Default: every position present.")
    parser.add_argument("--model-name", default=granite.MODEL_NAME)
    parser.add_argument("--local-files-only", action="store_true", help="Use the cached model; do not contact Hugging Face.")
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--device", help="cuda or cpu. Default: cuda if available.")
    return parser.parse_args()


if __name__ == "__main__":
    main()
