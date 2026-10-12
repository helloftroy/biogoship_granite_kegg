#!/usr/bin/env python3
"""Download the Sentinel-3 products listed in data/manifests/products.csv.

Fetches every candidate OLCI scene (all overpasses within the search window of a
site, not just the closest) and each scene's paired SST orbit file, as complete
product archives:

    data/products/olci/<product name>.zip
    data/products/sst/<product name>.zip

Nothing is tiled or deleted here; whole scenes are kept so that later steps can
read any location in them. Re-running resumes: archives already on disk are
skipped, failed ones are retried. Every attempt is recorded in
data/download_log.csv.

Credentials: see cdse_credentials.sh.example.

Usage:
    python scripts/download_products.py --dry-run          # sizes only
    python scripts/download_products.py --max-products 2   # smoke test
    python scripts/download_products.py
"""

from __future__ import annotations

import argparse
import sys
import threading
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship import cdse  # noqa: E402
from biogoship.tables import MANIFEST_DIR, REPO_ROOT, append_csv, read_csv  # noqa: E402

LOG_FIELDS = ["product_id", "product_name", "kind", "status", "megabytes", "seconds", "message"]


def main() -> None:
    args = parse_args()
    wanted = wanted_products(read_csv(args.manifest_dir / "products.csv"), args)
    todo = [item for item in wanted if not item["dest"].exists()]
    total_gb = sum(item["size_bytes"] for item in wanted) / 1e9
    todo_gb = sum(item["size_bytes"] for item in todo) / 1e9
    for kind in ("olci", "sst"):
        of_kind = [item for item in wanted if item["kind"] == kind]
        on_disk = sum(1 for item in of_kind if item["dest"].exists())
        print(f"{kind.upper():4s}: {len(of_kind):,} products, {on_disk:,} already on disk")
    print(f"Total {total_gb:,.0f} GB; still to download {len(todo):,} products, about {todo_gb:,.0f} GB")
    print(f"Destination: {args.products_dir}")
    if args.dry_run or not todo:
        return

    tokens = cdse.TokenProvider()
    log_path = args.products_dir.parent / "download_log.csv"
    lock = threading.Lock()
    progress = {"done": 0, "failed": 0, "bytes": 0}
    started = time.time()

    def fetch(item: dict) -> None:
        began = time.time()
        status, message = "ok", ""
        try:
            cdse.download(
                cdse.product_zip_url(item["product_id"]),
                item["dest"],
                tokens,
                timeout=args.timeout,
                retries=args.retries,
                retry_sleep=args.retry_sleep,
            )
            if not zipfile.is_zipfile(item["dest"]):
                item["dest"].unlink(missing_ok=True)
                raise RuntimeError("downloaded file is not a zip archive")
        except Exception as exc:  # noqa: BLE001 - log it and move on; a re-run retries
            status, message = "error", str(exc)[:300]
        megabytes = item["dest"].stat().st_size / 1e6 if item["dest"].exists() else 0.0
        with lock:
            progress["done" if status == "ok" else "failed"] += 1
            progress["bytes"] += int(megabytes * 1e6)
            append_csv(
                log_path,
                [{**item, "status": status, "megabytes": f"{megabytes:.0f}", "seconds": f"{time.time() - began:.0f}", "message": message}],
                LOG_FIELDS,
            )
            finished = progress["done"] + progress["failed"]
            rate = progress["bytes"] / 1e6 / max(time.time() - started, 1)
            print(
                f"[{finished}/{len(todo)}] {status:5s} {item['kind']} {item['product_name']} "
                f"({megabytes:,.0f} MB; {progress['bytes'] / 1e9:,.1f} GB so far at {rate:.0f} MB/s) {message}",
                flush=True,
            )

    with ThreadPoolExecutor(max_workers=args.connections) as pool:
        list(pool.map(fetch, todo))

    print(f"\nDownloaded {progress['done']:,} products; {progress['failed']:,} failed (re-run to retry). Log: {log_path}")
    if progress["failed"] == len(todo):
        raise SystemExit("Every download failed; check credentials and network before re-running.")


def wanted_products(products: list[dict[str, str]], args: argparse.Namespace) -> list[dict]:
    """OLCI scenes in time order, each followed by its SST orbit file the first time it is needed."""
    if args.max_products:
        products = products[: args.max_products]
    wanted, seen = [], set()
    for row in products:
        pairs = [("olci", row["product_id"], row["product_name"], row["size_bytes"])]
        if not args.no_sst and row.get("sst_product_id"):
            pairs.append(("sst", row["sst_product_id"], row["sst_product_name"], row["sst_size_bytes"]))
        for kind, product_id, name, size in pairs:
            if product_id in seen:
                continue
            seen.add(product_id)
            wanted.append(
                {
                    "kind": kind,
                    "product_id": product_id,
                    "product_name": name,
                    "size_bytes": float(size or 0),
                    "dest": args.products_dir / kind / f"{name}.zip",
                }
            )
    return wanted


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest-dir", type=Path, default=MANIFEST_DIR)
    parser.add_argument("--products-dir", type=Path, default=REPO_ROOT / "data" / "products")
    parser.add_argument("--max-products", type=int, default=0, help="Only the first N OLCI scenes (0 = all).")
    parser.add_argument("--no-sst", action="store_true", help="Skip the SST orbit files.")
    parser.add_argument("--dry-run", action="store_true", help="Report counts and sizes; download nothing.")
    parser.add_argument("--connections", type=int, default=4, help="Parallel downloads (CDSE allows 4 per account).")
    parser.add_argument("--timeout", type=float, default=1800.0, help="Socket timeout per download, seconds.")
    parser.add_argument("--retries", type=int, default=6)
    parser.add_argument("--retry-sleep", type=float, default=30.0)
    return parser.parse_args()


if __name__ == "__main__":
    main()
