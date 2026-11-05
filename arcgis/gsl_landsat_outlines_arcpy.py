"""
gsl_landsat_outlines_arcpy.py  (version 2: per-pixel cloud compositing)

Great Salt Lake outlines from Landsat imagery (1984 to present) in ArcGIS Pro.

Imagery source: Esri Living Atlas "Landsat Level-2" imagery layer
(Landsat 4, 5, 7, 8, 9 Collection 2 Level-2 surface reflectance, 1982-present).
It needs an ArcGIS Online organizational account but does NOT use credits.
ArcGIS Pro must be signed in to that account when you run this.

Why version 2: the scene-level cloud-cover number describes the whole
185-km scene, not the lake, and around Great Salt Lake it is often inflated
because the cloud mask flags bright salt flats and playas as cloud. Version 1
also needed one low-cloud scene for every part of the lake. Now clouds are
removed pixel by pixel from every scene, and each window's outline is built
from all the clear looks at each pixel. A cloudy scene still contributes its
clear pixels, so the scene-level limit can be high (CLOUD_MAX = 80).
Quality is judged where it matters, along the shoreline (edge_gap field).

For each time window the script:
  1. finds Tier-1 scenes over the lake and orders them least-cloudy first,
  2. downloads green, SWIR-1 and Pixel QA for each scene (30 m, UTM 12N),
  3. masks fill, cloud, cloud shadow, cirrus (and dilated cloud) from Pixel QA,
  4. computes MNDWI = (green - SWIR1) / (green + SWIR1) and classifies water
     with an Otsu threshold set separately for each scene,
  5. adds scenes until 98% of the area has at least 2 clear looks,
  6. calls a pixel water if it was water in at least half of its clear looks,
  7. keeps water connected to the seed points, fills interior holes that are
     mostly unobserved (real islands stay), and vectorizes the result.

Output (in OUT_FOLDER):
  GSL_landsat.gdb/GSL_outlines_polygons   all windows, largest lake first
  GSL_landsat.gdb/GSL_outlines_lines      same, as outlines
  GSL_Landsat_outlines.shp / _lines.shp   shapefile copies
  GSL_Landsat_scenes.csv                  scenes used for each window

Quality fields:
  edge_gap   % of shoreline pixels that border unobserved (cloud) pixels.
             This is where the outline follows a cloud edge, not the shore.
  lake_gap   % of the lake area filled in because it was never observed
  quality    good (edge_gap <= 2%), check (<= 10%), poor (> 10%)

Run from the ArcGIS Pro Python window, a Pro notebook, or propy.bat.
No Spatial Analyst needed (uses numpy and scipy, both included with Pro).
"""

import csv
import datetime as dt
import json
import math
import os
import re
import shutil
import time

import numpy as np
import requests
from scipy import ndimage
import arcpy

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------
# Work in a LOCAL folder. OneDrive locks geodatabase files while it syncs them,
# which makes ArcGIS fail with ERROR 000210 ("Cannot create output").
OUT_FOLDER = r"C:\GSL\landsat_outlines"
# When the run finishes, the shapefiles and scene list are copied here
# (shapefiles sync safely; the geodatabase stays local). None = no copy.
COPY_RESULTS_TO = None   # e.g. r"C:\Users\you\OneDrive\GSL\landsat_outlines"
# Reuse years already in the geodatabase instead of downloading them again.
SKIP_EXISTING = True
# Years to redo even if they are already in the geodatabase, e.g. [2016, 2022].
RERUN_YEARS = []

FIRST_YEAR, LAST_YEAR = 1984, 2026

# One window per year at the seasonal low (the lake usually bottoms out in
# October-November). Month-day strings.
WINDOW_START, WINDOW_END = "09-01", "11-15"

# Extra windows, e.g. the 1986-87 record-high stand. Add or remove freely.
EXTRA_WINDOWS = [
    ("1986 record-high stand", "1986-05-15", "1986-07-31"),
    ("1987 spring high",       "1987-04-15", "1987-06-30"),
]

CLOUD_MAX = 80          # % scene cloud cover; only drops scenes that are nearly all cloud
MAX_SCENES = 12         # most scenes used for one window
MIN_CLEAR_OBS = 2       # stop adding scenes once STOP_FRAC of the area
STOP_FRAC = 0.98        #   has at least this many clear looks
USE_L7_FILL = True      # add Landsat 7 after the other sensors if gaps remain
TIER_1_ONLY = True      # Tier-1 scenes have the best geolocation
MASK_DILATED_CLOUD = True   # also mask the buffer the QA band draws around clouds

WATER_FRAC = 0.5        # water if water in >= this share of a pixel's clear looks
THRESHOLD = "otsu"      # or a number such as 0.0
OTSU_MIN, OTSU_MAX = -0.15, 0.30

CELL = 30               # m
# Area searched (lon/lat); wide enough for the 1986-87 high stand.
REGION = (-113.35, 40.55, -111.85, 41.80)

# Points that are open water in every year. Any water body touching one is
# kept. Add the bays if you want them kept when they separate at low levels.
SEEDS = [
    (-112.45, 41.05),   # south arm (Gilbert Bay)
    (-112.683, 41.323), # north arm (Gunnison Bay), >= 10 km from shore even in low years
    # (-112.05, 40.97), # Farmington Bay
    # (-112.30, 41.35), # Bear River Bay
]

FILL_CLOUD_HOLES = True  # fill lake holes that are mostly unobserved
MIN_HOLE_KM2 = 0.05      # remove smaller interior holes (0 = keep all)
EDGE_GAP_GOOD, EDGE_GAP_CHECK = 2.0, 10.0
KEEP_RASTERS = False     # keep lake, water-frequency and clear-count rasters

SERVICE = ("https://landsat.imagery1.arcgis.com/arcgis/rest/services/"
           "LandsatC2L2/ImageServer")
SR_OUT = arcpy.SpatialReference(32612)   # WGS 84 / UTM 12N
MAX_TILE_PX = 1000     # pixels per exportImage request (~12 MB)
MIN_TILE_PX = 250      # smallest tile when splitting after server errors

# Service bands (1-based): 3 Green, 6 SWIR 1, 8 Pixel QA -> 0-based ids
BAND_IDS = "2,5,7"
NODATA = -9999.0
# Pixel QA bits: 0 fill, 1 dilated cloud, 2 cirrus, 3 cloud, 4 cloud shadow
QA_BAD = 0b11101 | (0b00010 if MASK_DILATED_CLOUD else 0)
# ---------------------------------------------------------------------------

CROSS = ndimage.generate_binary_structure(2, 1)   # 4-connectivity


def log(msg):
    print(msg, flush=True)


def remove(path):
    try:
        if path.lower().endswith(".tif") and arcpy.Exists(path):
            arcpy.management.Delete(path)
        elif os.path.exists(path):
            os.remove(path)
    except Exception:
        pass            # a locked temp file is harmless


def session():
    tok = arcpy.GetSigninToken()
    if not tok or not tok.get("token"):
        raise RuntimeError("ArcGIS Pro is not signed in to an ArcGIS Online "
                           "organizational account.")
    s = requests.Session()
    if tok.get("referer"):
        s.headers["Referer"] = tok["referer"]
    s.token = tok["token"]          # sent in the POST body, not the URL
    return s


def epoch_ms(day, end_of_day=False):
    d = dt.datetime.strptime(day, "%Y-%m-%d").replace(tzinfo=dt.timezone.utc)
    if end_of_day:
        d += dt.timedelta(days=1) - dt.timedelta(milliseconds=1)
    return int(d.timestamp() * 1000)


def sensors_for(year):
    if year <= 2011:
        return ["LT05", "LT04"]
    if year == 2012:
        return ["LE07"]             # only Landsat 7 was imaging in 2012
    if year <= 2020:
        return ["LC08"]
    return ["LC08", "LC09"]


def wgs_polygon(xmin, ymin, xmax, ymax):
    pts = [arcpy.Point(xmin, ymin), arcpy.Point(xmin, ymax),
           arcpy.Point(xmax, ymax), arcpy.Point(xmax, ymin)]
    return arcpy.Polygon(arcpy.Array(pts), arcpy.SpatialReference(4326))


AOI = wgs_polygon(*REGION)


def make_grid():
    ext = AOI.projectAs(SR_OUT).extent
    xmin = math.floor(ext.XMin / CELL) * CELL
    ymax = math.ceil(ext.YMax / CELL) * CELL
    ncol = int(math.ceil((ext.XMax - xmin) / CELL))
    nrow = int(math.ceil((ymax - ext.YMin) / CELL))
    return xmin, ymax, ncol, nrow


# ------------------------------ catalog ------------------------------------
PID = re.compile(r"^L[CETOM]0\d_L\w{3}_\d{6}_\d{8}_\d{8}_\d{2}_(T1|T2|RT)$")


def query_scenes(s, start, end):
    """Scenes over the lake in a date window. Matches the query Esri's own
    Landsat Explorer sends: category = 1 keeps real scenes (not overview
    tiles), dates go in the where clause, the product ID is in 'name', and
    cloudcover is a fraction (0-1)."""
    where = (f"category = 1 AND acquisitiondate BETWEEN "
             f"timestamp '{start} 00:00:00' AND timestamp '{end} 23:59:59'")
    params = {
        "f": "json", "where": where, "token": s.token,
        "geometry": json.dumps({"xmin": REGION[0], "ymin": REGION[1],
                                "xmax": REGION[2], "ymax": REGION[3],
                                "spatialReference": {"wkid": 4326}}),
        "geometryType": "esriGeometryEnvelope",
        "spatialRel": "esriSpatialRelIntersects",
        "outFields": ("objectid,name,landsat_product_id,acquisitiondate,"
                      "cloudcover,category"),
        "returnGeometry": "true", "outSR": 4326,
        "orderByFields": "acquisitiondate", "resultRecordCount": 1000,
    }
    feats, offset = [], 0
    while True:
        r = s.post(SERVICE + "/query", data={**params, "resultOffset": offset},
                   timeout=180)
        r.raise_for_status()
        js = r.json()
        if "error" in js:
            raise RuntimeError(f"Catalog query failed: {js['error']}")
        page = js.get("features", [])
        feats += page
        if not js.get("exceededTransferLimit") or not page:
            break
        offset += len(page)
    if not feats:
        log("  Catalog returned no scenes. Raw response: "
            + json.dumps(js)[:300])

    scenes, bad_id = [], []
    for f in feats:
        a = {k.lower(): v for k, v in f["attributes"].items()}
        pid = next((str(a.get(k) or "") for k in ("name", "landsat_product_id")
                    if PID.match(str(a.get(k) or ""))), "")
        rings = (f.get("geometry") or {}).get("rings")
        if not pid or not rings:
            bad_id.append(str(a.get("name")))
            continue
        a["landsat_product_id"] = pid
        arr = arcpy.Array([arcpy.Array([arcpy.Point(x, y) for x, y in ring])
                           for ring in rings])
        a["geom"] = arcpy.Polygon(arr, arcpy.SpatialReference(4326))
        cc = a.get("cloudcover")
        a["cloudcover"] = 100.0 if cc is None or cc < 0 else float(cc)
        scenes.append(a)
    if bad_id:
        log(f"  {len(bad_id)} records ignored (no product ID or footprint), "
            f"e.g. {bad_id[:3]}")
    # cloudcover is a fraction on this service; convert to percent
    if scenes and max(sc["cloudcover"] for sc in scenes) <= 1.0:
        for sc in scenes:
            sc["cloudcover"] *= 100.0
    return scenes


def order_scenes(scenes, year):
    prim = sensors_for(year)
    def ok(sc, prefixes):
        pid = sc["landsat_product_id"]
        return (pid[:4] in prefixes and sc["cloudcover"] <= CLOUD_MAX
                and (not TIER_1_ONLY or pid.endswith("_T1")))
    first = sorted([sc for sc in scenes if ok(sc, prim)],
                   key=lambda sc: sc["cloudcover"])
    extra = []
    if USE_L7_FILL and "LE07" not in prim:
        extra = sorted([sc for sc in scenes if ok(sc, ["LE07"])],
                       key=lambda sc: sc["cloudcover"])
    return first + extra


# ------------------------------ download -----------------------------------
def scene_bounds(sc, grid):
    xmin0, ymax0, ncol, nrow = grid
    g = sc["geom"].intersect(AOI, 4)
    if g is None or g.area == 0:
        return None
    ext = g.projectAs(SR_OUT).extent
    c0 = max(0, int(math.floor((ext.XMin - xmin0) / CELL)))
    c1 = min(ncol, int(math.ceil((ext.XMax - xmin0) / CELL)))
    r0 = max(0, int(math.floor((ymax0 - ext.YMax) / CELL)))
    r1 = min(nrow, int(math.ceil((ymax0 - ext.YMin) / CELL)))
    return None if (r1 <= r0 or c1 <= c0) else (r0, r1, c0, c1)


_SIMPLE_PARAMS = False      # set once the full parameter set has failed


class TileError(RuntimeError):
    pass


def _export(s, sc, x0, y1, w, h, simple):
    """One exportImage request; returns a (3, h, w) array or raises TileError."""
    oid = int(sc["objectid"])
    rule = {"mosaicMethod": "esriMosaicLockRaster", "lockRasterIds": [oid],
            "where": f"objectid in ({oid})", "ascending": False}
    params = {
        "f": "image", "format": "tiff", "pixelType": "F32",
        "bbox": f"{x0},{y1 - h * CELL},{x0 + w * CELL},{y1}",
        "bboxSR": SR_OUT.factoryCode, "imageSR": SR_OUT.factoryCode,
        "size": f"{w},{h}", "bandIds": BAND_IDS,
        "renderingRule": json.dumps({"rasterFunction": "None"}),
        "mosaicRule": json.dumps(rule), "token": s.token,
    }
    if not simple:
        params.update(noData=NODATA, interpolation="RSP_NearestNeighbor",
                      compression="LZ77")
    try:
        r = s.post(SERVICE + "/exportImage", data=params, timeout=300)
    except requests.RequestException as e:
        raise TileError(f"request failed: {e}")
    if r.status_code != 200 or r.content[:4] not in (b"II*\x00", b"MM\x00*"):
        raise TileError(f"HTTP {r.status_code}: {r.text[:300].strip()}")
    tif = os.path.join(OUT_FOLDER, "rasters", f"_tile_{os.getpid()}.tif")
    remove(tif)
    with open(tif, "wb") as fh:
        fh.write(r.content)
    try:
        a = arcpy.RasterToNumPyArray(tif, nodata_to_value=NODATA)
    finally:
        remove(tif)
    if a.ndim != 3 or a.shape[0] < 3 or a.shape[1:] != (h, w):
        raise TileError(f"unexpected tile shape {a.shape}")
    return a[:3].astype(np.float32)


def _fetch(s, sc, out, grid, rr, cc, h, w, origin):
    """Fill out[:, rr:rr+h, cc:cc+w] (grid offsets origin), retrying and
    splitting the tile on server errors."""
    global _SIMPLE_PARAMS
    xmin0, ymax0, _, _ = grid
    x0, y1 = xmin0 + cc * CELL, ymax0 - rr * CELL
    err = None
    for attempt in range(2):
        try:
            a = _export(s, sc, x0, y1, w, h, _SIMPLE_PARAMS)
            r0, c0 = origin
            out[:, rr - r0:rr - r0 + h, cc - c0:cc - c0 + w] = a
            return
        except TileError as e:
            err = e
            time.sleep(3 * (attempt + 1))
    if max(h, w) > MIN_TILE_PX:              # split into quarters and retry
        h2, w2 = (h + 1) // 2, (w + 1) // 2
        for dr, hh in ((0, h2), (h2, h - h2)):
            for dc, ww in ((0, w2), (w2, w - w2)):
                if hh > 0 and ww > 0:
                    _fetch(s, sc, out, grid, rr + dr, cc + dc, hh, ww, origin)
        return
    if not _SIMPLE_PARAMS:                   # last resort: fewer parameters
        log(f"    small tiles still fail ({err}); retrying with simpler request")
        _SIMPLE_PARAMS = True
        _fetch(s, sc, out, grid, rr, cc, h, w, origin)
        return
    raise RuntimeError(f"exportImage keeps failing for "
                       f"{sc['landsat_product_id']}: {err}")


def download_scene(s, sc, bounds, grid, workdir):
    r0, r1, c0, c1 = bounds
    out = np.full((3, r1 - r0, c1 - c0), NODATA, np.float32)
    for rr in range(r0, r1, MAX_TILE_PX):
        for cc in range(c0, c1, MAX_TILE_PX):
            _fetch(s, sc, out, grid, rr, cc, min(MAX_TILE_PX, r1 - rr),
                   min(MAX_TILE_PX, c1 - cc), (r0, c0))
    return out


# --------------------------- classification --------------------------------
def otsu(values):
    hist, edges = np.histogram(values, bins=256, range=(-1, 1))
    centers = (edges[:-1] + edges[1:]) / 2
    w0 = np.cumsum(hist).astype(float)
    w1 = w0[-1] - w0
    m0 = np.cumsum(hist * centers)
    with np.errstate(divide="ignore", invalid="ignore"):
        between = w0 * w1 * (m0 / w0 - (m0[-1] - m0) / w1) ** 2
    between[~np.isfinite(between)] = 0
    return float(edges[np.argmax(between) + 1])


def scene_water(green, swir, qa):
    """Clear-pixel mask, water mask and threshold for one scene."""
    valid = ((green > NODATA + 1) & (swir > NODATA + 1) & (qa > 0)
             & np.isfinite(green) & np.isfinite(swir) & np.isfinite(qa))
    if valid.sum() < 1000:
        return None
    if np.median(qa[valid]) < 1:
        raise RuntimeError("Pixel QA values look rescaled; cannot read cloud bits.")
    if np.median(green[valid]) > 2:           # raw DNs: apply USGS scale/offset
        green = green * 0.0000275 - 0.2
        swir = swir * 0.0000275 - 0.2
    qa_i = np.where(valid, np.rint(qa), 1).astype(np.int64)
    denom = green + swir
    ok = valid & ((qa_i & QA_BAD) == 0) & (np.abs(denom) > 1e-6)
    if ok.sum() < 1000:
        return None
    mndwi = np.zeros(green.shape, np.float32)
    mndwi[ok] = (green[ok] - swir[ok]) / denom[ok]
    if THRESHOLD == "otsu":
        thr = min(max(otsu(mndwi[ok][::3]), OTSU_MIN), OTSU_MAX)
    else:
        thr = float(THRESHOLD)
    return ok, ok & (mndwi > thr), thr


def build_lake(clear_n, water_n, seed_rc):
    """Lake mask from per-pixel counts; returns (lake, edge_gap %, lake_gap %)."""
    resolved = clear_n > 0
    water = resolved & (water_n.astype(np.float32) >=
                        WATER_FRAC * clear_n.astype(np.float32))
    lab, _ = ndimage.label(water, structure=CROSS)
    ids = set()
    nrow, ncol = lab.shape
    for n, (r, c) in enumerate(seed_rc):
        found = False
        if 0 <= r < nrow and 0 <= c < ncol:
            win = lab[max(r - 5, 0):r + 6, max(c - 5, 0):c + 6]
            nz = win[win > 0]
            if nz.size:
                ids.add(int(np.bincount(nz).argmax()))
                found = True
        if not found:
            log(f"  NOTE: seed {n + 1} {SEEDS[n]} is not in water this year; any part "
                "of the lake reachable only from it is left out.")
    if not ids:
        return None, None, None
    lake = np.isin(lab, sorted(ids))

    added = np.zeros_like(lake)
    if FILL_CLOUD_HOLES:
        holes = ndimage.binary_fill_holes(lake, structure=CROSS) & ~lake
        hl, hn = ndimage.label(holes, structure=CROSS)
        if hn:
            frac = np.asarray(ndimage.mean(~resolved, labels=hl,
                                           index=np.arange(1, hn + 1)))
            added = np.isin(hl, np.arange(1, hn + 1)[frac > 0.5])
    lake = lake | added

    # shoreline pixels, and those that border unobserved pixels
    unres_out = ~resolved & ~lake
    lp = np.pad(lake, 1, constant_values=False)
    up = np.pad(unres_out, 1, constant_values=False)
    perim = np.zeros_like(lake)
    gap = np.zeros_like(lake)
    for dr, dc in ((-1, 0), (1, 0), (0, -1), (0, 1)):
        sl = (slice(1 + dr, 1 + dr + nrow), slice(1 + dc, 1 + dc + ncol))
        perim |= ~lp[sl]
        gap |= up[sl]
    perim &= lake
    gap &= lake
    edge_gap = 100.0 * gap.sum() / max(1, perim.sum())
    lake_gap = 100.0 * added.sum() / max(1, lake.sum())
    return lake, edge_gap, lake_gap


# ------------------------------- output ------------------------------------
FIELDS = [("label", "TEXT", 60), ("start_dt", "TEXT", 10),
          ("end_dt", "TEXT", 10), ("year", "SHORT", None),
          ("n_scenes", "SHORT", None), ("sensors", "TEXT", 30),
          ("thr_mean", "DOUBLE", None), ("edge_gap", "DOUBLE", None),
          ("lake_gap", "DOUBLE", None), ("quality", "TEXT", 8),
          ("area_km2", "DOUBLE", None), ("area_mi2", "DOUBLE", None)]


def save_array(arr, path, grid, nodata):
    xmin0, ymax0, _, nrow = grid
    ras = arcpy.NumPyArrayToRaster(arr, arcpy.Point(xmin0, ymax0 - nrow * CELL),
                                   CELL, CELL, nodata)
    ras.save(path)
    arcpy.management.DefineProjection(path, SR_OUT)
    return path


def to_feature_class(lake_path, out_fc, attrs):
    tmp_p, tmp_d = r"memory\lp", r"memory\ld"
    arcpy.conversion.RasterToPolygon(lake_path, tmp_p, "SIMPLIFY", "Value")
    arcpy.management.Dissolve(tmp_p, tmp_d)
    if MIN_HOLE_KM2 > 0:
        arcpy.management.EliminatePolygonPart(
            tmp_d, out_fc, "AREA", f"{MIN_HOLE_KM2} SquareKilometers",
            part_option="CONTAINED_ONLY")
    else:
        arcpy.management.CopyFeatures(tmp_d, out_fc)
    for name, ftype, length in FIELDS:
        arcpy.management.AddField(out_fc, name, ftype, field_length=length)
    names = [f[0] for f in FIELDS]
    with arcpy.da.UpdateCursor(out_fc, ["SHAPE@"] + names) as cur:
        for row in cur:
            a = row[0].getArea("GEODESIC", "SQUAREKILOMETERS")
            attrs.update(area_km2=round(a, 2), area_mi2=round(a / 2.589988, 2))
            cur.updateRow([row[0]] + [attrs[n] for n in names])
    for fc in (tmp_p, tmp_d):
        arcpy.management.Delete(fc)
    return attrs["area_km2"]


def windows():
    out = [(f"Fall {y}", f"{y}-{WINDOW_START}", f"{y}-{WINDOW_END}")
           for y in range(FIRST_YEAR, LAST_YEAR + 1)]
    return out + list(EXTRA_WINDOWS)


# -------------------------------- main -------------------------------------
def gdb_writable(gdb):
    """Create and delete a scratch table to prove the geodatabase accepts writes."""
    t = os.path.join(gdb, "zz_write_test")
    try:
        if arcpy.Exists(t):
            arcpy.management.Delete(t)
        arcpy.management.CreateTable(gdb, "zz_write_test")
        arcpy.management.Delete(t)
        return True, ""
    except Exception as e:
        return False, str(e).strip().splitlines()[0]


def main():
    arcpy.env.overwriteOutput = True
    os.makedirs(OUT_FOLDER, exist_ok=True)
    workdir = os.path.join(OUT_FOLDER, "rasters")
    os.makedirs(workdir, exist_ok=True)
    gdb = os.path.join(OUT_FOLDER, "GSL_landsat.gdb")
    if not arcpy.Exists(gdb):
        arcpy.management.CreateFileGDB(OUT_FOLDER, "GSL_landsat.gdb")
    if "onedrive" in OUT_FOLDER.lower():
        log("NOTE: OUT_FOLDER is inside OneDrive. Pause OneDrive syncing for the "
            "whole run, or ArcGIS may fail with ERROR 000210.")
    ok, why = gdb_writable(gdb)
    if not ok:
        log(f"Cannot write to {gdb}:\n  {why}\nClose any ArcGIS Pro map or "
            "catalog view using it, pause OneDrive, or set OUT_FOLDER to a local folder.")
        return

    grid = make_grid()
    xmin0, ymax0, ncol, nrow = grid
    log(f"Grid: {ncol} x {nrow} cells at {CELL} m")
    seed_rc = []
    for lon, lat in SEEDS:
        p = arcpy.PointGeometry(arcpy.Point(lon, lat),
                                arcpy.SpatialReference(4326)).projectAs(SR_OUT)
        seed_rc.append((int((ymax0 - p.firstPoint.Y) / CELL),
                        int((p.firstPoint.X - xmin0) / CELL)))

    today = dt.date.today().isoformat()
    fcs, scene_log, stopped = [], [], False
    for label, start, end in windows():
        if start > today:
            continue
        end = min(end, today)
        year = int(start[:4])
        tag = label.lower().replace(" ", "_").replace("-", "")
        out_fc = os.path.join(gdb, f"w_{tag}")
        log(f"\n=== {label} ({start} to {end}) ===")
        if SKIP_EXISTING and arcpy.Exists(out_fc) and year not in RERUN_YEARS:
            log("  Already in the geodatabase; reusing it.")
            fcs.append(out_fc)
            continue
        ok, why = gdb_writable(gdb)             # check before spending time downloading
        if not ok:
            log(f"  STOPPED: the geodatabase no longer accepts writes:\n  {why}\n"
                "  Close anything using it, pause OneDrive or move OUT_FOLDER to a "
                "local folder, then run again; finished years are kept.")
            stopped = True
            break
        try:
            s = session()                       # fresh token each window
            scenes = query_scenes(s, start, end)
            order = order_scenes(scenes, year)
            if scenes:
                ccs = [sc["cloudcover"] for sc in scenes]
                allowed = sensors_for(year) + (["LE07"] if USE_L7_FILL else [])
                n_sens = sum(sc["landsat_product_id"][:4] in allowed for sc in scenes)
                n_tier = sum(sc["landsat_product_id"][:4] in allowed and
                             (not TIER_1_ONLY or sc["landsat_product_id"].endswith("_T1"))
                             for sc in scenes)
                log(f"  {len(scenes)} scenes found (cloud {min(ccs):.0f}-"
                    f"{max(ccs):.0f}%): {n_sens} right sensor, {n_tier} also "
                    f"Tier 1, {len(order)} also <= {CLOUD_MAX}% cloud")
            if not order:
                log("  Skipped: no usable scenes.")
                continue

            clear_n = np.zeros((nrow, ncol), np.uint8)
            water_n = np.zeros((nrow, ncol), np.uint8)
            used, thrs = [], []
            for sc in order:
                if len(used) >= MAX_SCENES:
                    break
                b = scene_bounds(sc, grid)
                if b is None:
                    continue
                arr = download_scene(s, sc, b, grid, workdir)
                res = scene_water(arr[0], arr[1], arr[2])
                if res is None:
                    continue
                ok, wat, thr = res
                r0, r1, c0, c1 = b
                clear_n[r0:r1, c0:c1] += ok.astype(np.uint8)
                water_n[r0:r1, c0:c1] += wat.astype(np.uint8)
                used.append(sc)
                thrs.append(thr)
                cov = float((clear_n >= MIN_CLEAR_OBS).mean())
                log(f"  + {sc['landsat_product_id']}  cloud {sc['cloudcover']:.0f}%"
                    f"  thr {thr:.3f}  -> {cov:.0%} of area has "
                    f">= {MIN_CLEAR_OBS} clear looks")
                if cov >= STOP_FRAC:
                    break
            if not used:
                log("  Skipped: no scene had clear pixels.")
                continue

            lake, edge_gap, lake_gap = build_lake(clear_n, water_n, seed_rc)
            if lake is None:
                log("  Skipped: no seed point fell in water.")
                continue
            quality = ("good" if edge_gap <= EDGE_GAP_GOOD else
                       "check" if edge_gap <= EDGE_GAP_CHECK else "poor")

            lake_path = save_array(lake.astype(np.uint8),
                                   os.path.join(workdir, f"{tag}_lake.tif"), grid, 0)
            if KEEP_RASTERS:
                frac = np.where(clear_n > 0, water_n / np.maximum(clear_n, 1), -1)
                save_array(frac.astype(np.float32),
                           os.path.join(workdir, f"{tag}_waterfreq.tif"), grid, -1)
                save_array(clear_n, os.path.join(workdir, f"{tag}_clearlooks.tif"),
                           grid, 255)

            ids = [sc["landsat_product_id"] for sc in used]
            attrs = {"label": label, "start_dt": start, "end_dt": end,
                     "year": year, "n_scenes": len(used),
                     "sensors": "+".join(sorted({i[:4] for i in ids})),
                     "thr_mean": round(float(np.mean(thrs)), 4),
                     "edge_gap": round(edge_gap, 2),
                     "lake_gap": round(lake_gap, 2), "quality": quality}
            area = to_feature_class(lake_path, out_fc, attrs)
            log(f"  area {area:,.0f} km2 ({area / 2.589988:,.0f} mi2); "
                f"edge gap {edge_gap:.1f}%, filled {lake_gap:.1f}% -> {quality}")
            fcs.append(out_fc)
            scene_log.append([label, start, end, len(used),
                              ";".join(f"{sc['landsat_product_id']} ({sc['cloudcover']:.0f}%)"
                                       for sc in used)])
            if not KEEP_RASTERS:
                remove(lake_path)
        except Exception as e:                  # keep going; report at end
            log(f"  WARNING: {label} failed: {e}")

    if stopped:
        return
    if not fcs:
        log("No outlines were produced.")
        return

    merged = r"memory\merged"
    polys = os.path.join(gdb, "GSL_outlines_polygons")
    lines = os.path.join(gdb, "GSL_outlines_lines")
    arcpy.management.Merge(fcs, merged)
    arcpy.management.Sort(merged, polys, [["area_km2", "DESCENDING"]])
    arcpy.management.PolygonToLine(polys, lines, "IGNORE_NEIGHBORS")
    arcpy.management.JoinField(lines, "ORIG_FID", polys, "OBJECTID",
                               [f[0] for f in FIELDS])
    arcpy.conversion.ExportFeatures(
        polys, os.path.join(OUT_FOLDER, "GSL_Landsat_outlines.shp"))
    arcpy.conversion.ExportFeatures(
        lines, os.path.join(OUT_FOLDER, "GSL_Landsat_outlines_lines.shp"))
    # scene list: keep rows from earlier runs, replace years redone in this one
    csv_path = os.path.join(OUT_FOLDER, "GSL_Landsat_scenes.csv")
    rows = {}
    if os.path.exists(csv_path):
        with open(csv_path, newline="") as fh:
            for r in list(csv.reader(fh))[1:]:
                if r:
                    rows[r[0]] = r
    for r in scene_log:
        rows[r[0]] = r
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["label", "start", "end", "n_scenes", "scenes (cloud %)"])
        w.writerows(sorted(rows.values(), key=lambda r: r[1]))
    log(f"\nDone: {len(fcs)} outlines in {OUT_FOLDER}")

    if COPY_RESULTS_TO:
        try:
            os.makedirs(COPY_RESULTS_TO, exist_ok=True)
            for name in ("GSL_Landsat_outlines", "GSL_Landsat_outlines_lines"):
                arcpy.conversion.ExportFeatures(os.path.join(OUT_FOLDER, name + ".shp"),
                                                os.path.join(COPY_RESULTS_TO, name + ".shp"))
            shutil.copy2(csv_path, COPY_RESULTS_TO)
            log(f"Copied shapefiles and scene list to {COPY_RESULTS_TO}")
        except Exception as e:
            log(f"Copy to {COPY_RESULTS_TO} failed ({e}); the results are in {OUT_FOLDER}.")


if __name__ == "__main__":
    main()
