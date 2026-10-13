"""Cut Granite-ready 42x42 tiles out of Sentinel-3 products, centred on a sample site.

A tile is a float32 array of shape (17, 42, 42): 16 OLCI water-leaving reflectance
bands (Oa01-Oa12, Oa16, Oa17, Oa18, Oa21) plus SLSTR sea surface temperature, each
z-scored with the per-channel statistics of the ibm-granite/granite-geospatial-ocean
pretraining data. That layout and normalisation are what the Granite encoder expects.
"""

from __future__ import annotations

import re
from pathlib import Path

import h5py
import numpy as np

from .tables import wrap_lon_delta

TILE_SIZE = 42  # pixels; 12.6 km at OLCI full resolution
PIXEL_KM = 0.3  # nominal OLCI full-resolution pixel spacing
HALF = TILE_SIZE // 2

OLCI_BANDS = [
    "Oa01", "Oa02", "Oa03", "Oa04", "Oa05", "Oa06", "Oa07", "Oa08",
    "Oa09", "Oa10", "Oa11", "Oa12", "Oa16", "Oa17", "Oa18", "Oa21",
]  # fmt: skip
OLCI_BAND_FILES = [f"{band}_reflectance.nc" for band in OLCI_BANDS]
OLCI_CHECK_FILES = ["geo_coordinates.nc", "wqsf.nc"]  # enough to tell whether a site is cloud-free

# Granite Geospatial Ocean pretraining statistics, physical units, in OLCI_BANDS order.
OLCI_MEANS = np.array([
    0.0235427398,   0.0226303495,   0.0199877248,   0.0166938124,
    0.0119924026,   0.00767917988,  0.00251636861,  0.00189688827,
    0.00192718330,  0.00190564570,  0.00103529217,  0.000566890440,
    0.000595696267, 0.000402757423, 0.000423631744, 0.000105166233,
], dtype=np.float32)  # fmt: skip
OLCI_STDS = np.array([
    0.00776708,  0.00733259,  0.00633057,  0.00615707,
    0.00610327,  0.00663780,  0.00539699,  0.00511585,
    0.00507850,  0.00507704,  0.00484563,  0.00415998,
    0.00441236,  0.00408463,  0.00400387,  0.00370793,
], dtype=np.float32)  # fmt: skip
SST_MEAN_K = np.float32(293.908469)
SST_STD_K = np.float32(2.67577808)

# OLCI water quality and science flags (wqsf.nc, one uint64 per pixel).
WQSF_WATER = np.uint64(2)
WQSF_BITS = {
    "INVALID": 1,
    "LAND": 4,
    "CLOUD": 8,
    "SNOW_ICE": 16,
    "INLAND_WATER": 32,
    "COASTLINE": 64,
    "SUSPECT": 512,
    "SATURATED": 2048,
    "HIGHGLINT": 8192,
    "AC_FAIL": 262144,
    "OCNN_FAIL": 1048576,
    "CLOUD_AMBIGUOUS": 8388608,
    "CLOUD_MARGIN": 16777216,
}
# A pixel is usable open ocean when WATER is set and none of these are.
REJECT_FLAGS = [
    "INVALID", "LAND", "CLOUD", "SNOW_ICE", "INLAND_WATER", "COASTLINE",
    "SUSPECT", "SATURATED", "AC_FAIL", "CLOUD_AMBIGUOUS", "CLOUD_MARGIN",
]  # fmt: skip
WQSF_REJECT = np.uint64(sum(WQSF_BITS[name] for name in REJECT_FLAGS))

SLSTR_MIN_QUALITY = 4  # GHRSST quality_level: 4 acceptable, 5 best
SST_SEARCH_PAD_DEG = 0.1


def _attr(dataset: h5py.Dataset, name: str, default: float) -> float:
    return float(np.ravel(dataset.attrs.get(name, default))[0])


# ----------------------------------------------------------------------------
# OLCI
# ----------------------------------------------------------------------------
def load_geo(geo_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Per-pixel latitude and longitude (float32, NaN where missing)."""
    with h5py.File(geo_path, "r") as f:
        out = []
        for name in ("latitude", "longitude"):
            ds = f[name]
            raw = ds[:]
            values = raw.astype(np.float32) * np.float32(_attr(ds, "scale_factor", 1.0))
            values[raw == int(_attr(ds, "_FillValue", -2147483648))] = np.nan
            out.append(values)
    return out[0], out[1]


def load_wqsf(wqsf_path: Path) -> np.ndarray:
    with h5py.File(wqsf_path, "r") as f:
        return f["WQSF"][:]


def ocean_mask(wqsf: np.ndarray) -> np.ndarray:
    return ((wqsf & WQSF_WATER) != 0) & ((wqsf & WQSF_REJECT) == 0)


def dominant_reject_flags(wqsf_window: np.ndarray, top: int = 3) -> str:
    """Which flags spoiled a window, e.g. 'CLOUD=0.81;CLOUD_MARGIN=0.12'."""
    fractions = {name: float(((wqsf_window & np.uint64(bit)) != 0).mean()) for name, bit in WQSF_BITS.items()}
    ranked = sorted(((frac, name) for name, frac in fractions.items() if frac > 0), reverse=True)[:top]
    return ";".join(f"{name}={frac:.2f}" for frac, name in ranked)


def locate_pixel(lat: np.ndarray, lon: np.ndarray, site_lat: float, site_lon: float) -> tuple[int, int, float]:
    """Row, column and distance (km) of the pixel nearest to a site."""
    cos_lat = np.float32(np.cos(np.radians(site_lat)))
    d2 = (lat - np.float32(site_lat)) ** 2 + (wrap_lon_delta(lon - np.float32(site_lon)) * cos_lat) ** 2
    if not np.isfinite(d2).any():
        return -1, -1, float("inf")
    row, col = np.unravel_index(int(np.nanargmin(d2)), d2.shape)
    return int(row), int(col), float(np.sqrt(d2[row, col])) * 111.195


def window_origin(row: int, col: int, shape: tuple[int, int]) -> tuple[int, int] | None:
    """Top-left corner of the TILE_SIZE window centred on (row, col); None if it leaves the swath."""
    r0, c0 = row - HALF, col - HALF
    if r0 < 0 or c0 < 0 or r0 + TILE_SIZE > shape[0] or c0 + TILE_SIZE > shape[1]:
        return None
    return r0, c0


def find_clear_window(
    valid: np.ndarray, row: int, col: int, max_offset_px: int, min_valid: float
) -> tuple[int, int, float] | None:
    """Nearest TILE_SIZE window to (row, col) that is at least min_valid clear.

    Searches every window whose centre is within max_offset_px pixels and lies fully
    inside the scene. Returns (r0, c0, clear fraction) of the closest one, or None.
    """
    centre_r0, centre_c0 = row - HALF, col - HALF
    r_lo, r_hi = max(centre_r0 - max_offset_px, 0), min(centre_r0 + max_offset_px, valid.shape[0] - TILE_SIZE)
    c_lo, c_hi = max(centre_c0 - max_offset_px, 0), min(centre_c0 + max_offset_px, valid.shape[1] - TILE_SIZE)
    if r_lo > r_hi or c_lo > c_hi:
        return None
    region = valid[r_lo : r_hi + TILE_SIZE, c_lo : c_hi + TILE_SIZE]
    total = np.pad(region.astype(np.int32).cumsum(axis=0).cumsum(axis=1), ((1, 0), (1, 0)))  # summed-area table
    sums = total[TILE_SIZE:, TILE_SIZE:] - total[:-TILE_SIZE, TILE_SIZE:] - total[TILE_SIZE:, :-TILE_SIZE] + total[:-TILE_SIZE, :-TILE_SIZE]
    fraction = sums / float(TILE_SIZE * TILE_SIZE)
    rows, cols = np.mgrid[r_lo : r_hi + 1, c_lo : c_hi + 1]
    distance2 = (rows - centre_r0) ** 2 + (cols - centre_c0) ** 2
    usable = (fraction >= min_valid) & (distance2 <= max_offset_px**2)
    if not usable.any():
        return None
    # Closest first; among equally close windows, the clearest.
    order = np.lexsort((-fraction[usable], distance2[usable]))
    best = np.flatnonzero(usable)[order[0]]
    i, j = np.unravel_index(best, fraction.shape)
    return int(r_lo + i), int(c_lo + j), float(fraction[i, j])


def read_olci_window(product_dir: Path, r0: int, c0: int) -> np.ndarray:
    """Normalised (16, 42, 42) reflectance window. Fill pixels become 0, i.e. the channel mean."""
    bands = np.empty((len(OLCI_BANDS), TILE_SIZE, TILE_SIZE), dtype=np.float32)
    for i, (band, filename) in enumerate(zip(OLCI_BANDS, OLCI_BAND_FILES)):
        with h5py.File(product_dir / filename, "r") as f:
            ds = f[f"{band}_reflectance"]
            raw = ds[r0 : r0 + TILE_SIZE, c0 : c0 + TILE_SIZE]
            physical = raw.astype(np.float32) * np.float32(_attr(ds, "scale_factor", 1.0)) + np.float32(
                _attr(ds, "add_offset", 0.0)
            )
            normalised = (physical - OLCI_MEANS[i]) / OLCI_STDS[i]
            normalised[raw == int(_attr(ds, "_FillValue", 65535))] = 0.0
        bands[i] = normalised
    return bands


# ----------------------------------------------------------------------------
# SLSTR sea surface temperature
# ----------------------------------------------------------------------------
def is_sst_l2p_file(filename: str) -> bool:
    """The GHRSST L2P NetCDF is the one data file in an SL_2_WST product."""
    return filename.endswith(".nc") and re.search(r"\d{14}-\w+-L2P", filename) is not None


class SstSwath:
    """One orbit of SLSTR SST. Coordinates are held in memory; SST is read only where needed."""

    def __init__(self, l2p_path: Path) -> None:
        self.path = Path(l2p_path)
        with h5py.File(self.path, "r") as f:
            self.lat = f["lat"][:].astype(np.float32)
            self.lon = f["lon"][:].astype(np.float32)

    def resample(self, tile_lat: np.ndarray, tile_lon: np.ndarray) -> np.ndarray:
        """Nearest-neighbour SST in kelvin on the OLCI tile grid; NaN where no good SST pixel is near.

        SLSTR is ~1 km against OLCI's 300 m, so several OLCI pixels share one SST value.
        Nearest neighbour keeps quality boundaries sharp instead of blending across bad pixels.
        """
        empty = np.full(tile_lat.shape, np.nan, dtype=np.float32)
        centre_lat = float(np.nanmean(tile_lat))
        centre_lon = float(tile_lon[tile_lon.shape[0] // 2, tile_lon.shape[1] // 2])
        cos_lat = max(np.cos(np.radians(centre_lat)), 0.05)
        half_lat = float(np.nanmax(tile_lat) - np.nanmin(tile_lat)) / 2 + SST_SEARCH_PAD_DEG
        half_lon = float(np.nanmax(np.abs(wrap_lon_delta(tile_lon - centre_lon)))) + SST_SEARCH_PAD_DEG / cos_lat

        near = (np.abs(self.lat - centre_lat) <= half_lat) & (
            np.abs(wrap_lon_delta(self.lon - centre_lon)) <= half_lon
        )
        rows = np.flatnonzero(near.any(axis=1))
        if rows.size == 0:
            return empty
        r0, r1 = int(rows[0]), int(rows[-1]) + 1

        with h5py.File(self.path, "r") as f:
            ds = f["sea_surface_temperature"]
            raw = ds[0, r0:r1, :]
            sst = raw.astype(np.float32) * np.float32(_attr(ds, "scale_factor", 0.01)) + np.float32(
                _attr(ds, "add_offset", 273.15)
            )
            sst[raw == int(_attr(ds, "_FillValue", -32768))] = np.nan
            quality = f["quality_level"][0, r0:r1, :]

        good = near[r0:r1] & (quality >= SLSTR_MIN_QUALITY) & np.isfinite(sst)
        if not good.any():
            return empty
        cand_lat = self.lat[r0:r1][good]
        cand_lon = self.lon[r0:r1][good]
        cand_sst = sst[good]

        d2 = (tile_lat.reshape(-1, 1) - cand_lat) ** 2 + (
            wrap_lon_delta(tile_lon.reshape(-1, 1) - cand_lon) * np.float32(cos_lat)
        ) ** 2
        return cand_sst[np.argmin(d2, axis=1)].reshape(tile_lat.shape).astype(np.float32)


def normalise_sst(sst_kelvin: np.ndarray) -> np.ndarray:
    return ((sst_kelvin - SST_MEAN_K) / SST_STD_K).astype(np.float32)


def stack_tile(olci: np.ndarray, sst_kelvin: np.ndarray | None) -> np.ndarray:
    """(17, 42, 42) Granite input. A missing SST channel is NaN; the embedding step fills it."""
    if sst_kelvin is None:
        sst = np.full((1, TILE_SIZE, TILE_SIZE), np.nan, dtype=np.float32)
    else:
        sst = normalise_sst(sst_kelvin)[np.newaxis]
    return np.concatenate([olci, sst], axis=0).astype(np.float32)
