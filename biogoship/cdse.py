"""Copernicus Data Space Ecosystem (CDSE): authentication, catalogue search and downloads.

Catalogue searches and file listings need no account. Downloads need credentials in
the environment (see cdse_credentials.sh.example):

    CDSE_USERNAME / CDSE_PASSWORD     account login, or
    CDSE_ACCESS_TOKEN                 a pasted short-lived token

Optional knobs: CDSE_TOKEN_METHOD=urllib|curl, CDSE_TOKEN_TIMEOUT, CDSE_CONNECT_TIMEOUT,
CDSE_TOTP, and CDSE_SSL_NO_VERIFY=1 as a last resort when the Python certificate store
on an HPC node cannot verify CDSE.
"""

from __future__ import annotations

import json
import os
import socket
import ssl
import subprocess
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

CATALOGUE_URL = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"
DOWNLOAD_URL = "https://download.dataspace.copernicus.eu/odata/v1/Products"
TOKEN_URL = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"
CLIENT_ID = "cdse-public"

OLCI_WFR = "OL_2_WFR___"  # OLCI level-2 water full resolution (300 m)
SLSTR_WST = "SL_2_WST___"  # SLSTR level-2 sea surface temperature (GHRSST L2P, one file per orbit)

RETRY_STATUS = {429, 500, 502, 503, 504}
PAGE_SIZE = 1000
TOKEN_REFRESH_MARGIN_S = 90.0


# ----------------------------------------------------------------------------
# Authentication
# ----------------------------------------------------------------------------
def ssl_context() -> ssl.SSLContext:
    if os.getenv("CDSE_SSL_NO_VERIFY", "").lower() in {"1", "true", "yes"}:
        return ssl._create_unverified_context()
    try:
        import certifi

        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        return ssl.create_default_context()


class TokenProvider:
    """Hands out a valid access token, renewing it shortly before it expires.

    CDSE access tokens live about ten minutes. Renewal uses the refresh token when
    one is still valid, so a long job holds one login session instead of opening a
    new one per product. Safe to share between download threads.
    """

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._access = ""
        self._access_expiry = 0.0
        self._refresh = ""
        self._refresh_expiry = 0.0

    def get(self, force: bool = False) -> str:
        pasted = os.getenv("CDSE_ACCESS_TOKEN", "").strip().strip("\"'`")
        if pasted:
            return pasted
        with self._lock:
            now = time.time()
            if force or not self._access or now >= self._access_expiry - TOKEN_REFRESH_MARGIN_S:
                self._renew(now)
            return self._access

    def _renew(self, now: float) -> None:
        payload = None
        if self._refresh and now < self._refresh_expiry - TOKEN_REFRESH_MARGIN_S:
            try:
                payload = request_token({"grant_type": "refresh_token", "refresh_token": self._refresh})
            except RuntimeError:
                payload = None  # refresh token rejected; fall back to a full login
        if payload is None:
            payload = request_token(password_form())
        self._access = payload["access_token"]
        self._access_expiry = now + float(payload.get("expires_in", 600))
        self._refresh = payload.get("refresh_token", "")
        self._refresh_expiry = now + float(payload.get("refresh_expires_in", 0))


def password_form() -> dict[str, str]:
    username = os.getenv("CDSE_USERNAME")
    password = os.getenv("CDSE_PASSWORD")
    if not username or not password:
        raise RuntimeError("Set CDSE_USERNAME and CDSE_PASSWORD (or CDSE_ACCESS_TOKEN) before downloading.")
    form = {"grant_type": "password", "username": username, "password": password}
    if os.getenv("CDSE_TOTP"):
        form["totp"] = os.environ["CDSE_TOTP"]
    return form


def request_token(form: dict[str, str]) -> dict[str, Any]:
    data = urllib.parse.urlencode({"client_id": CLIENT_ID, **form}).encode("utf-8")
    timeout = float(os.getenv("CDSE_TOKEN_TIMEOUT", "300"))
    method = os.getenv("CDSE_TOKEN_METHOD", "urllib").lower()
    if method == "curl":
        return _request_token_curl(data, timeout)
    if method not in {"urllib", "python"}:
        raise RuntimeError("CDSE_TOKEN_METHOD must be 'urllib' or 'curl'.")
    request = urllib.request.Request(
        TOKEN_URL, data=data, headers={"Content-Type": "application/x-www-form-urlencoded"}, method="POST"
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout, context=ssl_context()) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")[:500]
        raise RuntimeError(f"CDSE token request failed HTTP {exc.code}: {body}") from exc
    except (urllib.error.URLError, TimeoutError, socket.timeout) as exc:
        raise RuntimeError(f"CDSE token request could not reach {TOKEN_URL} within {timeout:g}s: {exc}") from exc


def _request_token_curl(data: bytes, timeout: float) -> dict[str, Any]:
    connect_timeout = float(os.getenv("CDSE_CONNECT_TIMEOUT", str(min(timeout, 300.0))))
    cmd = [
        "curl", "-sS", "--fail", "--show-error",
        "--connect-timeout", f"{connect_timeout:g}",
        "--max-time", f"{timeout:g}",
        "-H", "Content-Type: application/x-www-form-urlencoded",
        "--data-binary", "@-",
        TOKEN_URL,
    ]  # fmt: skip
    try:
        result = subprocess.run(cmd, input=data, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError("CDSE_TOKEN_METHOD=curl requested, but curl is not available.") from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"CDSE token curl request failed: {exc.stderr.decode('utf-8', 'replace').strip()}") from exc
    try:
        return json.loads(result.stdout.decode("utf-8"))
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"CDSE token curl response was not JSON: {result.stdout[:300]!r}") from exc


# ----------------------------------------------------------------------------
# Catalogue
# ----------------------------------------------------------------------------
def http_json(url: str, timeout: float = 120.0, retries: int = 8, base_sleep: float = 15.0) -> dict:
    """GET a JSON document, backing off on rate limits, server errors and dropped connections."""
    request = urllib.request.Request(url, headers={"Accept": "application/json"})
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=timeout, context=ssl_context()) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            if exc.code not in RETRY_STATUS or attempt == retries:
                raise
            sleep_s = _retry_after(exc) or min(base_sleep * 2**attempt, 900.0)
        except (urllib.error.URLError, TimeoutError, socket.timeout, ConnectionError):
            if attempt == retries:
                raise
            sleep_s = min(base_sleep * 2**attempt, 900.0)
        time.sleep(sleep_s)
    raise AssertionError("unreachable")


def _retry_after(exc: urllib.error.HTTPError) -> float | None:
    try:
        return float(exc.headers.get("Retry-After", ""))
    except ValueError:
        return None


def odata_time(dt: datetime) -> str:
    return dt.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.000Z")


def search_products(
    product_type: str,
    start: datetime,
    end: datetime,
    point: tuple[float, float] | None = None,
    **http_kwargs,
) -> list[dict[str, str]]:
    """Sentinel-3 products of one type sensed in [start, end), optionally covering a (lon, lat) point.

    Only non-time-critical (NT) products are returned; where a scene was reprocessed,
    only its newest version is kept.
    """
    clauses = [
        "Collection/Name eq 'SENTINEL-3'",
        "Attributes/OData.CSC.StringAttribute/any(att:att/Name eq 'productType' "
        f"and att/OData.CSC.StringAttribute/Value eq '{product_type}')",
        f"ContentDate/Start ge {odata_time(start)}",
        f"ContentDate/Start lt {odata_time(end)}",
    ]
    if point is not None:
        lon, lat = point
        clauses.append(f"OData.CSC.Intersects(area=geography'SRID=4326;POINT({lon:.6f} {lat:.6f})')")
    rows: list[dict[str, str]] = []
    skip = 0
    while True:
        params = {
            "$filter": " and ".join(clauses),
            "$orderby": "ContentDate/Start",
            "$top": str(PAGE_SIZE),
            "$skip": str(skip),
            "$select": "Id,Name,ContentDate,ContentLength,Online",
        }
        page = http_json(f"{CATALOGUE_URL}?{urllib.parse.urlencode(params)}", **http_kwargs).get("value", [])
        rows.extend(product_row(item) for item in page)
        if len(page) < PAGE_SIZE:
            break
        skip += PAGE_SIZE
    return newest_versions([row for row in rows if row["timeliness"] == "NT"])


def product_row(item: dict) -> dict[str, str]:
    name = item["Name"]
    dates = item.get("ContentDate") or {}
    return {
        "product_id": item["Id"],
        "product_name": name,
        "product_type": name[4:15],
        "satellite": name[:3],
        "start_utc": dates.get("Start", ""),
        "end_utc": dates.get("End", ""),
        "size_bytes": str(item.get("ContentLength", "")),
        "online": str(item.get("Online", "")),
        "timeliness": name[88:90],
        "created": name[48:63],
    }


def newest_versions(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Drop older processings of the same scene (same satellite, type and sensing start)."""
    best: dict[tuple[str, str, str], dict[str, str]] = {}
    for row in rows:
        key = (row["satellite"], row["product_type"], row["product_name"][16:31])
        if key not in best or row["created"] > best[key]["created"]:
            best[key] = row
    return sorted(best.values(), key=lambda row: (row["start_utc"], row["product_name"]))


def list_files(product_id: str, product_name: str, **http_kwargs) -> dict[str, int]:
    """File name -> size in bytes for the files inside a product's .SEN3 folder."""
    url = f"{CATALOGUE_URL}({product_id})/Nodes({urllib.parse.quote(product_name, safe='')})/Nodes"
    nodes = http_json(url, **http_kwargs).get("result", [])
    return {node["Name"]: int(node.get("ContentLength") or 0) for node in nodes}


# ----------------------------------------------------------------------------
# Downloads
# ----------------------------------------------------------------------------
def file_url(product_id: str, product_name: str, filename: str) -> str:
    """One file from inside a product, so the whole 300-600 MB archive is not needed."""
    quote = lambda text: urllib.parse.quote(text, safe="")  # noqa: E731
    return f"{DOWNLOAD_URL}({product_id})/Nodes({quote(product_name)})/Nodes({quote(filename)})/$value"


def product_zip_url(product_id: str) -> str:
    return f"{DOWNLOAD_URL}({product_id})/$value"


def download(
    url: str,
    dest: Path,
    tokens: TokenProvider,
    expected_bytes: int = 0,
    timeout: float = 1800.0,
    retries: int = 6,
    retry_sleep: float = 30.0,
) -> None:
    """Download url to dest (via a .part file), retrying with a fresh token on failure."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    error = ""
    for attempt in range(1, retries + 1):
        try:
            token = tokens.get(force=attempt > 1)
            request = urllib.request.Request(url, headers={"Authorization": f"Bearer {token}"})
            with urllib.request.urlopen(request, timeout=timeout, context=ssl_context()) as resp, part.open("wb") as out:
                while chunk := resp.read(1024 * 1024):
                    out.write(chunk)
            size = part.stat().st_size
            if expected_bytes and size != expected_bytes:
                raise RuntimeError(f"incomplete download: {size} of {expected_bytes} bytes")
            part.replace(dest)
            return
        except urllib.error.HTTPError as exc:
            error = f"HTTP {exc.code}: {exc.read().decode('utf-8', 'replace')[:300]}"
        except Exception as exc:  # noqa: BLE001 - network failures come in many types; all are retried
            error = f"{type(exc).__name__}: {exc}"
        part.unlink(missing_ok=True)
        if attempt < retries:
            print(f"    download attempt {attempt}/{retries} failed ({error}); retrying in {retry_sleep:g}s", flush=True)
            time.sleep(retry_sleep)
    raise RuntimeError(f"download failed after {retries} attempts: {dest.name}: {error}")
