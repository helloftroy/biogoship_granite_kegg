"""Deterministic advection of sampled water with surface currents, using OceanParcels.

One particle is released at each site at its sampling time and advected both
backward and forward in time with fourth-order Runge-Kutta. The resulting
trajectory says where that parcel of water was (or will be) at the moment of each
satellite overpass, so a tile can be cut there instead of at the sampling position.

Currents: Copernicus Marine MULTIOBS_GLO_PHY_MYNRT_015_003 (GlobCurrent total
velocity uo/vo, 0.25 degree, at 0 m and 15 m).

No diffusion and no ensemble: one trajectory per site. Requires parcels 3.x and xarray.
"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

from .tables import haversine_km, parse_utc, wrap_lon_delta

CURRENTS_PRODUCT = "MULTIOBS_GLO_PHY_MYNRT_015_003"
CURRENTS_DATASET = "cmems_obs-mob_glo_phy-cur_my_0.25deg_P1D-m"  # multi-year, daily means
CURRENT_VARIABLES = ["uo", "vo"]  # total eastward / northward velocity, m/s


def section_boxes(
    sites: list[dict[str, str]],
    site_products: list[dict[str, str]],
    scene_times: dict[str, datetime],
    margin_deg: float,
    pad_days: float,
) -> list[dict]:
    """Per section, the lon/lat box and time span the currents must cover.

    The box is the sites' extent plus margin_deg, so water can drift without leaving
    it. The time span runs from the earliest to the latest of all sampling and
    overpass times, plus pad_days at each end.
    """
    by_site = {site["site_id"]: site for site in sites}
    boxes: dict[str, dict] = {}
    for link in site_products:
        site = by_site[link["site_id"]]
        box = boxes.setdefault(site["section"], {"section": site["section"], "sites": {}, "times": []})
        box["sites"][site["site_id"]] = site
        box["times"] += [parse_utc(site["time_utc"]), scene_times[link["product_id"]]]
    out = []
    for box in boxes.values():
        lats = [float(site["lat"]) for site in box["sites"].values()]
        lons = [float(site["lon"]) for site in box["sites"].values()]
        out.append(
            {
                "section": box["section"],
                "sites": list(box["sites"].values()),
                "lat_min": max(min(lats) - margin_deg, -89.5),
                "lat_max": min(max(lats) + margin_deg, 89.5),
                "lon_min": min(lons) - margin_deg,
                "lon_max": max(lons) + margin_deg,
                "start": min(box["times"]) - timedelta(days=pad_days),
                "end": max(box["times"]) + timedelta(days=pad_days),
            }
        )
    return sorted(out, key=lambda box: box["start"])


def section_filename(section: str) -> str:
    return "currents_" + "".join(ch if ch.isalnum() else "_" for ch in section) + ".nc"


def load_currents(path: Path, depth_m: float):
    """Open a currents subset as an xarray Dataset of (time, latitude, longitude) at one depth level.

    Land and missing values become zero velocity, so a particle that reaches the
    coast stops there instead of raising an error.
    """
    import xarray as xr

    ds = xr.open_dataset(path)[CURRENT_VARIABLES]
    for name in ("depth", "elevation"):
        if name in ds.dims:
            level = int(np.argmin(np.abs(np.abs(ds[name].values) - abs(depth_m))))
            ds = ds.isel({name: level}, drop=True)
    return ds.transpose("time", "latitude", "longitude").fillna(0.0).load()


def advect(
    currents,
    releases: list[dict],
    max_hours: float,
    dt_minutes: float = 20.0,
    output_minutes: float = 60.0,
) -> dict[str, np.ndarray]:
    """Trajectories for releases [{id, lat, lon, time (datetime)}], max_hours each way.

    Returns id -> array of rows (hours_from_release, lat, lon), sorted by time and
    including the release point at 0. A trajectory is cut short where the particle
    leaves the currents domain.
    """
    import parcels

    fieldset = parcels.FieldSet.from_xarray_dataset(
        currents,
        variables={"U": "uo", "V": "vo"},
        dimensions={"lon": "longitude", "lat": "latitude", "time": "time"},
        mesh="spherical",
    )
    origin = np.datetime64(currents["time"].values[0], "ns")  # Parcels counts seconds from the first field time
    release_s = np.array(
        [(np.datetime64(r["time"].astimezone(timezone.utc).replace(tzinfo=None), "ns") - origin) / np.timedelta64(1, "s") for r in releases]
    )

    def delete_lost(particle, fieldset, time):  # noqa: ARG001 - Parcels kernel signature
        if particle.state >= 50:  # any error state, e.g. out of the currents domain
            particle.delete()

    step = timedelta(minutes=output_minutes)
    tracks: dict[int, list[tuple[float, float, float]]] = {index: [] for index in range(len(releases))}
    for direction in (-1, 1):
        pset = parcels.ParticleSet(
            fieldset=fieldset,
            pclass=parcels.ScipyParticle,
            lon=[r["lon"] for r in releases],
            lat=[r["lat"] for r in releases],
            time=release_s,
        )
        index_of = {particle.id: index for index, particle in enumerate(pset)}
        # Particles start at their own release time; the clock runs until the last one has had max_hours.
        span_hours = (release_s.max() - release_s.min()) / 3600.0 + max_hours
        for _ in range(int(np.ceil(span_hours * 60.0 / output_minutes)) + 1):  # +1: reach max_hours for all
            if len(pset) == 0:
                break
            pset.execute(
                [parcels.AdvectionRK4, delete_lost],
                runtime=step,
                dt=direction * timedelta(minutes=dt_minutes),
                verbose_progress=False,
            )
            finished = []
            for particle in pset:
                index = index_of[particle.id]
                hours = (float(particle.time) - release_s[index]) / 3600.0
                tracks[index].append((hours, float(particle.lat), float(particle.lon)))
                finished.append(abs(hours) >= max_hours)
            if any(finished):  # stop integrating particles that have covered their window
                pset.remove_booleanvector(np.array(finished))

    out = {}
    for index, release in enumerate(releases):
        rows = np.array([(0.0, release["lat"], release["lon"]), *tracks[index]], dtype=float)
        rows = rows[np.unique(np.round(rows[:, 0], 6), return_index=True)[1]]  # sorted, one row per time
        out[release["id"]] = rows
    return out


def position_at(track: np.ndarray, hours: float) -> tuple[float, float] | None:
    """Linear interpolation of a track to a time offset; None if the track does not reach it."""
    offsets = track[:, 0]
    if hours < offsets[0] - 1e-6 or hours > offsets[-1] + 1e-6:
        return None
    lat = float(np.interp(hours, offsets, track[:, 1]))
    unwrapped = track[0, 2] + np.cumsum(np.r_[0.0, wrap_lon_delta(np.diff(track[:, 2]))])
    lon = float(wrap_lon_delta(np.interp(hours, offsets, unwrapped)))
    return lat, lon


def path_length_km(track: np.ndarray, hours: float) -> float:
    """Distance travelled along the track between release and a time offset."""
    lo, hi = sorted((0.0, hours))
    inside = track[(track[:, 0] >= lo) & (track[:, 0] <= hi)]
    ends = [position_at(track, lo), position_at(track, hi)]
    if ends[0] is None or ends[1] is None:
        return float("nan")
    points = [ends[0], *[(row[1], row[2]) for row in inside], ends[1]]
    return sum(haversine_km(a[0], a[1], b[0], b[1]) for a, b in zip(points[:-1], points[1:]))


# ----------------------------------------------------------------------------
# Tracks on disk, and positions at overpass times
# ----------------------------------------------------------------------------
POSITION_FIELDS = [
    "site_id",
    "product_id",
    "time_gap_hours",
    "site_lat",
    "site_lon",
    "advected_lat",
    "advected_lon",
    "displacement_km",
    "path_km",
    "status",
]


def load_tracks(trajectories_csv: Path) -> dict[str, np.ndarray]:
    """site_id -> rows (hours_from_sample, lat, lon), from the file backtrack_sites.py writes."""
    import csv

    rows: dict[str, list[tuple[float, float, float]]] = {}
    with Path(trajectories_csv).open(newline="") as f:
        for row in csv.DictReader(f):
            rows.setdefault(row["site_id"], []).append(
                (float(row["hours_from_sample"]), float(row["lat"]), float(row["lon"]))
            )
    return {site_id: np.array(sorted(track), dtype=float) for site_id, track in rows.items()}


def advected_position_row(site: dict[str, str], link: dict[str, str], track: np.ndarray | None) -> dict:
    """Where the site's water was at one scene's sensing time.

    status: ok | stalled (never moved: no current data at the site) | left_domain (track
    ends before the overpass time) | no_track (the site was not tracked).
    """
    lat, lon, hours = float(site["lat"]), float(site["lon"]), float(link["time_gap_hours"])
    row = {
        "site_id": site["site_id"],
        "product_id": link["product_id"],
        "time_gap_hours": link["time_gap_hours"],
        "site_lat": site["lat"],
        "site_lon": site["lon"],
        "status": "no_track",
    }
    if track is None:
        return row
    where = position_at(track, hours)
    if where is None:
        row["status"] = "left_domain"
        return row
    displacement = haversine_km(lat, lon, where[0], where[1])
    row.update(
        advected_lat=f"{where[0]:.6f}",
        advected_lon=f"{where[1]:.6f}",
        displacement_km=f"{displacement:.3f}",
        path_km=f"{path_length_km(track, hours):.3f}",
        status="stalled" if displacement == 0.0 and abs(hours) >= 1.0 else "ok",
    )
    return row


def advected_position_rows(
    sites: list[dict[str, str]], site_products: list[dict[str, str]], tracks: dict[str, np.ndarray]
) -> list[dict]:
    by_site = {site["site_id"]: site for site in sites}
    return [advected_position_row(by_site[link["site_id"]], link, tracks.get(link["site_id"])) for link in site_products]
