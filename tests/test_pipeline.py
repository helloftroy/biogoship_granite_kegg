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
from biogoship import advection, cdse, tiles  # noqa: E402
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

    def test_find_nearest_clear_window(self):
        valid = np.ones(SHAPE, dtype=bool)
        valid[:, 110:] = False  # cloud bank over the eastern columns
        # Target at column 125 is under cloud; a window is 80% clear once it starts at column 76 or less.
        self.assertEqual(tiles.find_clear_window(valid, 60, 125, 60, 0.8)[:2], (39, 76))
        self.assertIsNone(tiles.find_clear_window(valid, 60, 125, 20, 0.8))  # nothing clear within 20 pixels
        self.assertEqual(tiles.find_clear_window(valid, 60, 60, 60, 0.8), (39, 39, 1.0))  # already clear: stays centred
        self.assertEqual(tiles.find_clear_window(valid, 5, 60, 60, 0.8)[:2], (0, 39))  # near the top edge: shifted inside

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


class AdvectionTests(unittest.TestCase):
    def test_position_along_track(self):
        track = np.array([[-2.0, 0.0, 179.8], [0.0, 0.0, 180.0], [2.0, 1.0, -179.8]])  # crosses the antimeridian
        lat, lon = advection.position_at(track, 1.0)
        self.assertAlmostEqual(lat, 0.5)
        self.assertAlmostEqual(abs(lon), 179.9, places=6)
        self.assertAlmostEqual(advection.position_at(track, -1.0)[1], 179.9, places=6)
        self.assertIsNone(advection.position_at(track, 3.0))
        self.assertAlmostEqual(advection.path_length_km(track, -2.0), 22.24, delta=0.05)

    def test_uniform_current_gives_known_drift(self):
        try:
            import parcels  # noqa: F401
            import xarray as xr
        except ImportError:
            self.skipTest("parcels/xarray not installed")
        import logging
        import warnings
        from datetime import datetime, timezone

        warnings.simplefilter("ignore")
        logging.getLogger("parcels").setLevel(logging.ERROR)
        lat, lon = np.arange(-5, 5.01, 0.25), np.arange(-35, -24.99, 0.25)
        times = np.arange("2020-03-25", "2020-04-02", dtype="datetime64[D]").astype("datetime64[ns]")
        shape = (len(times), len(lat), len(lon))
        currents = xr.Dataset(
            {"uo": (("time", "latitude", "longitude"), np.full(shape, 0.5, "float32")),
             "vo": (("time", "latitude", "longitude"), np.zeros(shape, "float32"))},
            coords={"time": times, "latitude": lat, "longitude": lon},
        )  # fmt: skip
        release = {"id": "s", "lat": 0.0, "lon": -30.0, "time": datetime(2020, 3, 28, 6, 30, tzinfo=timezone.utc)}
        track = advection.advect(currents, [release], max_hours=24)["s"]
        for hours, sign in ((24.0, 1), (-24.0, -1)):  # 0.5 m/s east for a day is 43.2 km
            where = advection.position_at(track, hours)
            self.assertAlmostEqual(haversine_km(0.0, -30.0, *where), 43.2, delta=0.2)
            self.assertGreater(sign * (where[1] + 30.0), 0)


class EndToEndTests(unittest.TestCase):
    """download_tiles on a fake scene: fetch() is a no-op because the files are already in place."""

    def test_clear_cloudy_edge_and_out_of_swath_sites(self):
        with tempfile.TemporaryDirectory() as tmp:
            data_dir = Path(tmp)
            args = argparse.Namespace(
                data_dir=data_dir, download_mode="files", no_sst=False, min_valid=0.8, keep_products=False,
                sst_cache=2, connections=1, products_dir=None, max_offset_km=0.0, timeout=1, retries=1, retry_sleep=0,
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
            # Water sampled at "clear" was 30 pixels east at overpass time; "cloudy" has no back-tracked position.
            advected = {("clear", "p1"): (LAT0 - 60 * STEP, LON0 + 90 * STEP)}
            targets = download_tiles.targets_for(sites, "p1", ["fixed", "advected"], advected)
            all_rows = worker.process(product, targets)
            rows = {row["site_id"]: row for row in all_rows if row["position"] == "fixed"}
            moved = {row["site_id"]: row for row in all_rows if row["position"] == "advected"}
            self.assertEqual(moved["clear"]["status"], "valid")
            self.assertEqual(moved["cloudy"]["status"], "no_position")
            with np.load(data_dir / "tiles" / "2020" / "03" / "tile_S3B_20200329T114129_clear_advected.npz") as tile:
                self.assertAlmostEqual(float(tile["lon_center"]), LON0 + 90 * STEP, places=4)
                self.assertEqual(str(tile["position"]), "advected")

            self.assertEqual(
                {name: row["status"] for name, row in rows.items()},
                {"clear": "valid", "cloudy": "cloudy", "edge": "edge", "outside": "out_of_swath"},
            )
            self.assertIn("CLOUD", rows["cloudy"]["reject_flags"])
            tile_file = data_dir / "tiles" / "2020" / "03" / "tile_S3B_20200329T114129_clear_fixed.npz"
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



class LocalArchiveRunTests(unittest.TestCase):
    """The scripts run as on the cluster: scenes already downloaded as archives, no network."""

    def run_script(self, module, argv: list[str]) -> None:
        old = sys.argv
        sys.argv = [module.__name__, *argv]
        try:
            module.main()
        finally:
            sys.argv = old

    def test_tile_manifest_and_embeddings(self):
        import contextlib
        import io
        import zipfile

        import build_tile_manifest
        from biogoship.tables import latest_tile_log_rows, read_csv, write_csv

        sst_name = "S3B_SL_2_WST____20200329T103000_20200329T121000_20200330T090904_6059_037_137______MAR_O_NT_003.SEN3"
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifests, data = root / "manifests", root / "data"
            # Scene archives, laid out like real products: <name>.zip containing <name>/<files>.
            write_olci(root / "src" / OLCI_NAME, cloudy_cols=slice(110, 160))
            write_sst(root / "src" / sst_name / SST_FILE)
            for kind, name in (("olci", OLCI_NAME), ("sst", sst_name)):
                (data / "products" / kind).mkdir(parents=True)
                with zipfile.ZipFile(data / "products" / kind / f"{name}.zip", "w") as archive:
                    for path in (root / "src" / name).iterdir():
                        archive.write(path, arcname=f"{name}/{path.name}")

            def site(name: str, col: int) -> dict:
                return {"site_id": name, "time_utc": "2020-03-29T06:00:00Z", "lat": LAT0 - 60 * STEP,
                        "lon": LON0 + col * STEP, "section": "TEST"}  # fmt: skip

            write_csv(manifests / "biogoship_sites.csv", [site("clear", 60), site("cloudy", 125)])
            write_csv(manifests / "biogoship_samples.csv", [
                {"sample_id": f"{s}_{rep}", "site_id": s, "collection_time_utc": "2020-03-29T06:00:00Z",
                 "lat": LAT0 - 60 * STEP, "lon": LON0 + col * STEP, "section": "TEST", "metagenome_runs": "SRR1"}
                for s, col in (("clear", 60), ("cloudy", 125)) for rep in ("a", "b")
            ])  # fmt: skip
            write_csv(manifests / "products.csv", [{
                "product_id": "p1", "product_name": OLCI_NAME, "satellite": "S3B", "start_utc": "2020-03-29T11:41:29.000000Z",
                "sst_product_id": "s1", "sst_product_name": sst_name,
            }])  # fmt: skip
            write_csv(manifests / "site_products.csv", [
                {"site_id": s, "product_id": "p1", "time_gap_hours": "5.69"} for s in ("clear", "cloudy")
            ])  # fmt: skip
            paths = ["--manifest-dir", str(manifests), "--data-dir", str(data), "--products-dir", str(data / "products")]
            common = [*paths, "--max-offset-km", "0"]  # centred windows only, to begin with

            with contextlib.redirect_stdout(io.StringIO()):
                # First the sampling positions only, as when back-tracking has not run yet.
                self.run_script(download_tiles, [*common, "--positions", "fixed"])
                self.assertEqual(len(read_csv(data / "tile_log.csv")), 2)
                # A log from before the position column existed is converted, not tripped over.
                legacy = [{k: v for k, v in row.items() if k not in ("position", "target_lat", "target_lon")}
                          for row in read_csv(data / "tile_log.csv")]  # fmt: skip
                write_csv(data / "tile_log.csv", legacy)
                self.run_script(download_tiles, [*common, "--positions", "fixed"])
                self.assertEqual([r["position"] for r in read_csv(data / "tile_log.csv")], ["fixed", "fixed"])
                self.assertTrue((data / "tile_log.legacy.csv").exists())
                # Then both: only the advected positions are new work. "cloudy" drifted west into clear water.
                write_csv(manifests / "advected_positions.csv", [
                    {"site_id": "clear", "product_id": "p1", "advected_lat": LAT0 - 60 * STEP,
                     "advected_lon": LON0 + 500 * STEP, "status": "ok"},
                    {"site_id": "cloudy", "product_id": "p1", "advected_lat": LAT0 - 60 * STEP,
                     "advected_lon": LON0 + 70 * STEP, "status": "ok"},
                ])  # fmt: skip
                self.run_script(download_tiles, common)
                log = read_csv(data / "tile_log.csv")
                self.assertEqual(
                    {(r["site_id"], r["position"]): r["status"] for r in log},
                    {("clear", "fixed"): "valid", ("cloudy", "fixed"): "cloudy",
                     ("clear", "advected"): "out_of_swath", ("cloudy", "advected"): "valid"},
                )  # fmt: skip
                self.run_script(download_tiles, common)  # nothing left to do
                self.assertEqual(len(read_csv(data / "tile_log.csv")), 4)
                self.assertTrue((data / "products" / "olci" / f"{OLCI_NAME}.zip").exists())  # archives are kept

                out_csv = manifests / "biogoship_tile_manifest.csv"
                self.run_script(build_tile_manifest, ["--manifest-dir", str(manifests), "--tile-log",
                                                      str(data / "tile_log.csv"), "--out-csv", str(out_csv)])  # fmt: skip
            manifest = read_csv(out_csv)
            self.assertEqual(sorted((r["sample_id"], r["position"]) for r in manifest),
                             [("clear_a", "fixed"), ("clear_b", "fixed"), ("cloudy_a", "advected"), ("cloudy_b", "advected")])  # fmt: skip
            flags = {(r["site_id"], r["position"]): r["flag"] for r in read_csv(manifests / "site_tile_status.csv")}
            self.assertEqual(flags, {("clear", "fixed"): "ok", ("clear", "advected"): "outside_all_scenes",
                                     ("cloudy", "fixed"): "all_cloudy", ("cloudy", "advected"): "ok"})  # fmt: skip

            # Allowing a nearby clear window rescues the cloudy sampling position, about 9 km west of it.
            import cloud_report

            with contextlib.redirect_stdout(io.StringIO()) as printed:
                self.run_script(download_tiles, [*paths, "--recheck", "cloudy", "edge"])  # default: within 20 km
                self.run_script(cloud_report, ["--manifest-dir", str(manifests), "--tile-log", str(data / "tile_log.csv")])
            latest = {(r["site_id"], r["position"]): r for r in latest_tile_log_rows(read_csv(data / "tile_log.csv"))}
            self.assertEqual(latest[("cloudy", "fixed")]["status"], "valid")
            self.assertTrue(7.0 < float(latest[("cloudy", "fixed")]["offset_km"]) < 11.0)
            self.assertEqual(latest[("clear", "fixed")]["offset_km"], "0.000")
            self.assertIn("Tiles taken from a nearby clear window: 1 of 2", printed.getvalue())

            try:
                import torch

                import embed_tiles
                from biogoship import granite
            except ImportError:
                self.skipTest("torch not installed; embedding step not exercised")

            class StubEncoder(torch.nn.Module):  # stands in for Granite: 441 tokens of width 512 per tile
                def forward(self, x):
                    self.seen_shape = tuple(x.shape)
                    return [x.mean(dim=(1, 2, 3, 4)).reshape(-1, 1, 1).expand(-1, 441, granite.EMBED_DIM)]

            stub, real = StubEncoder(), granite.build_encoder
            granite.build_encoder = lambda *args, **kwargs: (stub, "stub encoder")
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    self.run_script(embed_tiles, ["--manifest-csv", str(out_csv), "--out-dir", str(data / "embeddings")])
            finally:
                granite.build_encoder = real
            embeddings = np.load(data / "embeddings" / "tile_embeddings.npy")
            self.assertEqual(embeddings.shape, (2, granite.EMBED_DIM))  # two distinct tiles, four sample rows
            self.assertEqual(stub.seen_shape, (2, 17, 1, 42, 42))
            sample_rows = read_csv(data / "embeddings" / "sample_embeddings.csv")
            self.assertEqual(len(sample_rows), 4)
            self.assertEqual({r["embedding_row"] for r in sample_rows}, {"0", "1"})
            self.assertTrue((data / "embeddings" / "granite_embeddings_bundle.zip").exists())


if __name__ == "__main__":
    unittest.main()
