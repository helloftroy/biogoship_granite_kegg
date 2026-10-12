"""Offline tests: synthetic Sentinel-3 files stand in for real downloads.

Run with:  python -m unittest discover -s tests
"""

from __future__ import annotations

import argparse
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import download_tiles  # noqa: E402
import search_products  # noqa: E402
from biogoship import cdse, tiles  # noqa: E402
from biogoship.tables import haversine_km  # noqa: E402

SHAPE = (120, 160)
LAT0, LON0, STEP = 10.0, -30.0, 0.003  # ~300 m pixels
OLCI_NAME = "S3B_OL_2_WFR____20200329T114129_20200329T114429_20210701T233619_0180_037_137______MAR_R_NT_003.SEN3"
SST_FILE = "20200329103000-MAR-L2P_GHRSST-SSTskin-SLSTRB-20200330090904-v02.0-fv01.0.nc"


def grid() -> tuple[np.ndarray, np.ndarray]:
    rows, cols = np.mgrid[0 : SHAPE[0], 0 : SHAPE[1]]
    return LAT0 - rows * STEP, LON0 + cols * STEP  # north-up, like a descending daytime pass


def write_olci(product_dir: Path, cloudy_cols: slice | None = None) -> None:
    product_dir.mkdir(parents=True)
    lat, lon = grid()
    with h5py.File(product_dir / "geo_coordinates.nc", "w") as f:
        for name, values in (("latitude", lat), ("longitude", lon)):
            ds = f.create_dataset(name, data=np.round(values * 1e6).astype(np.int32))
            ds.attrs["scale_factor"] = np.array([1e-6])
            ds.attrs["_FillValue"] = np.array([-2147483648], dtype=np.int32)
    wqsf = np.full(SHAPE, tiles.WQSF_WATER, dtype=np.uint64)
    if cloudy_cols is not None:
        wqsf[:, cloudy_cols] |= np.uint64(tiles.WQSF_BITS["CLOUD"])
    with h5py.File(product_dir / "wqsf.nc", "w") as f:
        f.create_dataset("WQSF", data=wqsf)
    for i, (band, filename) in enumerate(zip(tiles.OLCI_BANDS, tiles.OLCI_BAND_FILES)):
        with h5py.File(product_dir / filename, "w") as f:
            # Physical reflectance equal to the channel mean, so the normalised value is 0...
            raw = np.full(SHAPE, round(float(tiles.OLCI_MEANS[i]) / 2e-5), dtype=np.uint16)
            raw[0, 0] = 65535
            ds = f.create_dataset(f"{band}_reflectance", data=raw)
            ds.attrs["scale_factor"] = np.array([2e-5], dtype=np.float32)
            ds.attrs["add_offset"] = np.array([0.0], dtype=np.float32)
            ds.attrs["_FillValue"] = np.array([65535], dtype=np.uint16)


def write_sst(path: Path) -> None:
    """A coarser (1 km) grid over the same area; SST in kelvin equals 290 + latitude."""
    path.parent.mkdir(parents=True, exist_ok=True)
    rows, cols = np.mgrid[0:60, 0:80]
    lat = (LAT0 + 0.05 - rows * 0.009).astype(np.float32)
    lon = (LON0 - 0.05 + cols * 0.009).astype(np.float32)
    kelvin = 290.0 + lat
    with h5py.File(path, "w") as f:
        f.create_dataset("lat", data=lat)
        f.create_dataset("lon", data=lon)
        ds = f.create_dataset("sea_surface_temperature", data=np.round((kelvin - 273.15) / 0.01).astype(np.int16)[None])
        ds.attrs["scale_factor"] = np.array([0.01], dtype=np.float32)
        ds.attrs["add_offset"] = np.array([273.15], dtype=np.float32)
        ds.attrs["_FillValue"] = np.array([-32768], dtype=np.int16)
        f.create_dataset("quality_level", data=np.full((1, 60, 80), 5, dtype=np.int8))


class TileGeometryTests(unittest.TestCase):
    def test_locate_pixel_and_window(self):
        lat, lon = (a.astype(np.float32) for a in grid())
        row, col, km = tiles.locate_pixel(lat, lon, LAT0 - 60 * STEP, LON0 + 80 * STEP)
        self.assertEqual((row, col), (60, 80))
        self.assertLess(km, 0.05)
        self.assertEqual(tiles.window_origin(60, 80, SHAPE), (39, 59))
        self.assertIsNone(tiles.window_origin(5, 80, SHAPE))
        self.assertGreater(tiles.locate_pixel(lat, lon, LAT0 + 2.0, LON0)[2], 100)

    def test_locate_pixel_across_antimeridian(self):
        lon = np.array([[179.998, -179.999]], dtype=np.float32)
        lat = np.zeros_like(lon)
        self.assertEqual(tiles.locate_pixel(lat, lon, 0.0, -179.9995)[:2], (0, 1))
        self.assertLess(haversine_km(0, 179.999, 0, -179.999), 0.3)

    def test_reject_mask_matches_florida_pipeline(self):
        self.assertEqual(
            int(tiles.WQSF_REJECT), 1 | 4 | 8 | 16 | 32 | 64 | 512 | 2048 | 8388608 | 16777216 | 262144
        )
        flags = np.array([2, 2 | 8, 2 | 8192, 4], dtype=np.uint64)  # water, cloud, glint (kept), land
        self.assertEqual(tiles.ocean_mask(flags).tolist(), [True, False, True, False])

    def test_product_name_fields(self):
        row = cdse.product_row({"Id": "x", "Name": OLCI_NAME, "ContentDate": {"Start": "s", "End": "e"}})
        self.assertEqual(
            (row["satellite"], row["product_type"], row["timeliness"], row["created"]),
            ("S3B", "OL_2_WFR___", "NT", "20210701T233619"),
        )
        older = {**row, "product_id": "old", "created": "20200401T000000"}
        self.assertEqual([r["product_id"] for r in cdse.newest_versions([older, row])], ["x"])


class SstPairingTests(unittest.TestCase):
    def test_prefers_same_satellite_orbit_containing_the_scene(self):
        olci = {"satellite": "S3B", "start_utc": "2020-03-29T11:41:29Z"}
        sst = [
            {"product_id": "a", "product_name": "A", "satellite": "S3A", "size_bytes": "1",
             "start_utc": "2020-03-29T11:00:00Z", "end_utc": "2020-03-29T12:40:00Z"},
            {"product_id": "b", "product_name": "B", "satellite": "S3B", "size_bytes": "1",
             "start_utc": "2020-03-29T10:30:00Z", "end_utc": "2020-03-29T12:10:00Z"},
        ]  # fmt: skip
        self.assertEqual(search_products.pair_sst(olci, sst)["sst_product_id"], "b")
        self.assertEqual(search_products.pair_sst(olci, sst[:1])["sst_relation"], "cross_satellite")
        far = [{**sst[1], "start_utc": "2020-03-30T10:30:00Z", "end_utc": "2020-03-30T12:10:00Z"}]
        self.assertEqual(search_products.pair_sst(olci, far)["sst_product_id"], "")


class EndToEndTests(unittest.TestCase):
    """download_tiles on a fake scene: fetch() is a no-op because the files are already in place."""

    def test_clear_cloudy_edge_and_out_of_swath_sites(self):
        with tempfile.TemporaryDirectory() as tmp:
            data_dir = Path(tmp)
            args = argparse.Namespace(
                data_dir=data_dir, download_mode="files", no_sst=False, min_valid=0.8, keep_products=False,
                sst_cache=2, connections=1, timeout=1, retries=1, retry_sleep=0,
            )  # fmt: skip
            worker = download_tiles.TileWorker(args)
            write_olci(worker.product_dir(OLCI_NAME), cloudy_cols=slice(110, 160))
            write_sst(worker.sst_dir / "SST_PRODUCT" / SST_FILE)
            product = {
                "product_id": "p1", "product_name": OLCI_NAME, "satellite": "S3B",
                "start_utc": "2020-03-29T11:41:29.000000Z", "sst_product_id": "s1", "sst_product_name": "SST_PRODUCT",
            }  # fmt: skip

            def site(name: str, row: float, col: float) -> tuple[dict, dict]:
                return (
                    {"site_id": name, "lat": str(LAT0 - row * STEP), "lon": str(LON0 + col * STEP),
                     "time_utc": "2020-03-29T06:00:00Z"},
                    {"site_id": name, "time_gap_hours": "5.69"},
                )  # fmt: skip

            sites = [site("clear", 60, 60), site("cloudy", 60, 125), site("edge", 5, 60), site("outside", 60, 900)]
            rows = {row["site_id"]: row for row in worker.process(product, sites)}

            self.assertEqual(
                {name: row["status"] for name, row in rows.items()},
                {"clear": "valid", "cloudy": "cloudy", "edge": "edge", "outside": "out_of_swath"},
            )
            self.assertIn("CLOUD", rows["cloudy"]["reject_flags"])
            tile_file = data_dir / "tiles" / "2020" / "03" / "tile_S3B_20200329T114129_clear.npz"
            self.assertTrue(tile_file.exists())
            with np.load(tile_file) as tile:
                data = tile["data"]
                self.assertEqual(data.shape, (17, 42, 42))
                self.assertEqual(data.dtype, np.float32)
                np.testing.assert_allclose(data[:16], 0.0, atol=0.01)  # reflectance == channel mean
                expected_sst = (290.0 + (LAT0 - 60 * STEP) - tiles.SST_MEAN_K) / tiles.SST_STD_K
                np.testing.assert_allclose(data[16].mean(), expected_sst, atol=0.01)
                self.assertTrue(tile["valid_mask"].all())
                self.assertAlmostEqual(float(tile["lat_center"]), LAT0 - 60 * STEP, places=4)
                self.assertEqual(str(tile["site_id"]), "clear")
                self.assertEqual(str(tile["date"]), "2020-03-29T11:41:29Z")


if __name__ == "__main__":
    unittest.main()
