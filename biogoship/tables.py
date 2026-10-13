"""Small CSV, time and geometry helpers shared by the pipeline scripts."""

from __future__ import annotations

import csv
import math
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
MANIFEST_DIR = REPO_ROOT / "data" / "manifests"
EARTH_RADIUS_KM = 6371.0088


def read_csv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = fieldnames or (list(rows[0].keys()) if rows else [])
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def append_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    """Append rows, writing the header only when the file is new. Flushes so a killed job loses nothing."""
    if not rows:
        return
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    is_new = not path.exists() or path.stat().st_size == 0
    with path.open("a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        if is_new:
            writer.writeheader()
        writer.writerows(rows)
        f.flush()


TILE_LOG_FIELDS = [
    "site_id",
    "product_id",
    "product_name",
    "position",  # fixed | advected
    "target_lat",
    "target_lon",
    "status",  # valid | cloudy | edge | out_of_swath | no_position
    "offset_km",  # distance from the target position to the tile centre; ~0 unless a nearby clear window was used
    "valid_fraction",
    "reject_flags",
    "time_gap_hours",
    "pixel_km",
    "row",
    "col",
    "tile_lat",
    "tile_lon",
    "sst_valid_fraction",
    "tile_path",
]
# Column layouts the log has had. Rows are recognised by their number of fields, so a
# file holding rows from more than one version still reads correctly.
_TILE_LOG_LAYOUTS = [
    TILE_LOG_FIELDS,
    [name for name in TILE_LOG_FIELDS if name != "offset_km"],
    [name for name in TILE_LOG_FIELDS if name not in ("offset_km", "position", "target_lat", "target_lon")],
]
_TILE_LOG_BY_WIDTH = {len(layout): layout for layout in _TILE_LOG_LAYOUTS}


def read_tile_log(path: Path) -> list[dict[str, str]]:
    """All rows of data/tile_log.csv in the current column layout, whatever version wrote them."""
    rows = []
    with Path(path).open(newline="") as f:
        reader = csv.reader(f)
        next(reader, None)  # the header only describes the rows of the version that created the file
        for line_number, values in enumerate(reader, start=2):
            if not values:
                continue
            layout = _TILE_LOG_BY_WIDTH.get(len(values))
            if layout is None:
                raise ValueError(f"{path}, line {line_number}: {len(values)} fields; expected one of {sorted(_TILE_LOG_BY_WIDTH)}")
            row = {name: "" for name in TILE_LOG_FIELDS}
            row.update(zip(layout, values))
            row["position"] = row["position"] or "fixed"  # logs from before positions existed
            rows.append(row)
    return rows


def append_tile_log(path: Path, rows: list[dict]) -> None:
    """Append rows, first bringing an older-layout file up to date so header and rows always agree."""
    repair_tile_log(path)
    append_csv(path, rows, TILE_LOG_FIELDS)


def repair_tile_log(path: Path) -> bool:
    """Rewrite the log in the current layout if any of it is not; returns True if it was rewritten.

    The file as it was is kept once, as tile_log.before_layout_update.csv.
    """
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return False
    with path.open(newline="") as f:
        reader = csv.reader(f)
        header = next(reader, [])
        widths = {len(values) for values in reader if values}
    if header == TILE_LOG_FIELDS and widths <= {len(TILE_LOG_FIELDS)}:
        return False
    rows = read_tile_log(path)
    backup = path.with_name(f"{path.stem}.before_layout_update.csv")
    if backup.exists():
        path.unlink()
    else:
        path.replace(backup)
    write_csv(path, rows, TILE_LOG_FIELDS)
    return True


def latest_tile_log_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """One row per (site, scene, position): the most recent, since a re-check appends a new row."""
    latest = {(row["site_id"], row["product_id"], row["position"]): row for row in rows}
    return list(latest.values())


def parse_utc(value: str) -> datetime:
    return datetime.fromisoformat(value.strip().replace("Z", "+00:00")).astimezone(timezone.utc)


def format_utc(value: datetime) -> str:
    return value.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def wrap_lon_delta(delta_deg):
    """Longitude difference folded into [-180, 180), so the antimeridian is not a 360 degree jump."""
    return (delta_deg + 180.0) % 360.0 - 180.0


def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = phi2 - phi1
    dlam = math.radians(wrap_lon_delta(lon2 - lon1))
    a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlam / 2) ** 2
    return 2 * EARTH_RADIUS_KM * math.asin(min(1.0, math.sqrt(a)))


def repo_relative(path: Path) -> str:
    """Path as stored in manifests: relative to the repo root, so manifests move between machines."""
    path = Path(path).resolve()
    try:
        return path.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(path)


def from_repo(path_text: str) -> Path:
    path = Path(path_text)
    return path if path.is_absolute() else REPO_ROOT / path
