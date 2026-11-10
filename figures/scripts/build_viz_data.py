"""Build the data the interactive page draws, from the outline shapefile and prep_outlines.py's
outlines.json, with three fixes and the water needed to reach 4,207 ft (historic, long-term average level).

1. Causeway: the railroad causeway splits the lake into two polygons whose shared edge draws as a
   straight line across the lake. Inside a box around the causeway only, the two arms are joined
   (morphological closing), so the drawing follows the lake's outer shoreline. Areas are unchanged.
2. Gaps: 2012 (Landsat 7 SLC-off stripes) is closed across the stripes with its own pixels;
   2024 has thin scene-seam slivers closed the same way; 2017 is missing the north arm west of a
   scene edge, which is filled from the year whose imaged north-arm shoreline matches 2017's best.
   Shoreline that was filled rather than imaged is kept separately so the page can draw it dotted.
3. Fill: every lake part is triangulated (the old build only filled the first ring, so years with the
   north arm as a separate polygon had no north-arm fill).
Water needed: volume between the fall level and 4,207 ft, from the USGS bathymetric areas at each
year's fall level (repo CSV), integrated over elevation with a piecewise-linear area curve.

Usage: python build_viz_data.py <outlines.shp> <levels.csv> <outlines.json> <viz_data.json> [outlines_fixed.json]
outlines_fixed.json (optional) is outlines.json with the fixed rings and areas, for make_static.py.
"""
import json, sys, csv, math
from shapely import affinity, constrained_delaunay_triangles
from shapely.geometry import Polygon, MultiPolygon, LineString, box
from shapely.ops import unary_union
import geopandas as gpd

shp, levels_csv, page_in, out = sys.argv[1:5]
fixed_out = sys.argv[5] if len(sys.argv) > 5 else None
X0, Y0 = 330000, 4490000
TOL_KM, MIN_PART, MIN_HOLE = 0.06, 5.0, 1.0      # same as prep_outlines.py (m -> km)
SMOOTH_KM, MIN_HOLE_LINE, FILL_SMOOTH, MIN_HOLE_FILL = 0.10, 3.0, 0.4, 20.0   # build_viz_data.py
REF_FT = 4207.0          # historic (long-term) average lake level, NGVD29 ft; must match prep_outlines.py
FT_KM = 0.0003048

def km(g): return affinity.affine_transform(g, [1e-3, 0, 0, 1e-3, -X0 / 1e3, -Y0 / 1e3])
def parts(g):
    if g.is_empty: return []
    return list(g.geoms) if g.geom_type in ("MultiPolygon", "GeometryCollection") else [g]
def polys_only(g): return unary_union([p for p in parts(g) if p.geom_type in ("Polygon", "MultiPolygon")])
def close(g, d): return g.buffer(d, join_style="mitre", mitre_limit=3).buffer(-d, join_style="mitre", mitre_limit=3)
def open_(g, d): return g.buffer(-d, join_style="mitre", mitre_limit=3).buffer(d, join_style="mitre", mitre_limit=3)

gdf = gpd.read_file(shp)
gdf = gdf[gdf.label.str.match(r"Fall \d{4}$")]
raw = {int(r.year): km(r.geometry).buffer(0) for r in gdf.itertuples()}
meta = {int(r.year): r for r in gdf.itertuples()}
page = json.load(open(page_in))
pyears = {d["year"]: d for d in page["years"]}

CAUSEWAY = box(13.0, 73.6, 46.5, 76.4)
NORTH = box(-15.0, 75.6, 46.0, 130.0)          # north arm (Gunnison Bay), above the causeway
R2012 = box(-8.0, 77.0, 25.0, 124.0)           # SLC-off striped part of 2012

fixed, filled_zone, notes = {}, {}, {}
for y, g in raw.items():
    g2, zone = g, None
    if y == 2012:
        g2 = unary_union([g, close(g, 0.8).intersection(R2012)]).buffer(0)
        notes[y] = "Landsat 7 striping (SLC-off) across the north arm, closed across the stripes from this year's own pixels"
    if y == 2024:
        g2 = unary_union([g, close(g, 0.45).intersection(NORTH)]).buffer(0)
        notes[y] = "thin scene-seam gaps in the north arm closed from this year's own pixels"
    fixed[y] = g2

# Donor fill: the year (not poor, fall level within 2 ft) whose north arm best matches the target's
# reliably imaged north arm (symmetric difference inside `match`, normalised) supplies the water
# that is missing from the target. Only wide missing pieces are added (opening removes the thin
# fringe a slightly different level would add along an imaged shore).
elev = {int(r["year"]): float(r["elev_ngvd29"]) for r in csv.DictReader(open(levels_csv))}
DONORS = {}
def donor_fill(yt, match, where):
    t = fixed[yt].intersection(match)
    imaged = t.buffer(1.0).intersection(match)
    best = None
    for y, g in raw.items():
        if y in (2012, 2017, yt) or meta[y].quality == "poor" or abs(elev[y] - elev[yt]) > 2: continue
        score = g.intersection(imaged).symmetric_difference(t).area / t.area
        if best is None or score < best[0]: best = (score, y)
    donor = raw[best[1]].intersection(where)
    miss = open_(donor.difference(fixed[yt]), 0.6).intersection(donor)
    miss = unary_union([p for p in parts(miss) if p.area > 5])
    fixed[yt] = unary_union([fixed[yt], miss.buffer(0.05).intersection(donor)]).buffer(0)
    DONORS[yt] = best[1]
    print(yt, "donor", best[1], f"mismatch {best[0]:.3f}", f"added {miss.area:.0f} km2")
    return best

# 2012: Landsat 7 SLC-off; the north arm west of x~22 is striped (gaps reach the shoreline).
b = donor_fill(2012, NORTH.difference(R2012), R2012)
notes[2012] = (f"Landsat 7 striping (SLC-off) across the north arm; stripes closed from this year's own pixels "
               f"and the unimaged shoreline filled with the {DONORS[2012]} shoreline, the year within 2 ft whose "
               f"imaged north arm matches 2012's best")
# 2017: north arm west of a Landsat scene edge (and its north tip) was never imaged.
b = donor_fill(2017, NORTH, NORTH)
notes[2017] = (f"north arm west of a Landsat scene edge was not imaged; filled with the {DONORS[2017]} shoreline, "
               f"the year within 2 ft whose imaged north arm matches 2017's best")
DONOR17 = DONORS[2017]

ratio = sorted(meta[y].area_km2 / raw[y].area for y in raw)[len(raw) // 2]   # geodesic / planar

# ---- water needed: integrate USGS area (km2) over level (ft) from the fall level to 4,207 ft ----
lev = sorted(((float(r["elev_ngvd29"]), float(r["usgs_total"]), int(r["year"])) for r in csv.DictReader(open(levels_csv))))
curve = sorted({e: a for e, a, _ in lev}.items())
def area_at(e):
    for (e0, a0), (e1, a1) in zip(curve, curve[1:]):
        if e0 <= e <= e1: return a0 + (a1 - a0) * (e - e0) / (e1 - e0)
    raise ValueError(e)
def vol_km3(e0, e1, n=400):   # km3 between levels e0 < e1
    h = (e1 - e0) / n
    return sum((area_at(e0 + i * h) + area_at(e0 + (i + 1) * h)) / 2 * h * FT_KM for i in range(n))
M3_AF, M3_GAL, POOL_M3 = 1233.48184, 264.172052, 2500.0   # acre-foot, US gallons, Olympic pool (50x25x2 m)
def need(e):
    if e >= REF_FT: return None
    m3 = vol_km3(e, REF_FT) * 1e9
    return {"m3": round(m3), "af": round(m3 / M3_AF), "gal": round(m3 * M3_GAL), "pools": round(m3 / POOL_M3)}
levels = [{"year": yy, "elev_ft": round(e, 2), "dev_ft": round(e - REF_FT, 2), "need": need(e)} for e, _, yy in sorted(lev, key=lambda r: r[2])]
print(f"USGS area at {REF_FT:.0f}:", round(area_at(REF_FT)), "km2")

# ---- drawing data ----
def r2(v): return round(v, 2)
def line_xy(geoms):
    xs, ys = [], []
    for g in geoms:
        for ln in ([g] if g.geom_type in ("LineString", "LinearRing") else list(getattr(g, "geoms", []))):
            if ln.is_empty or ln.length < 0.05: continue
            cx, cy = ln.coords.xy
            xs += [r2(v) for v in cx] + [None]; ys += [r2(v) for v in cy] + [None]
    return xs, ys

out_years, fixed_years = [], []
for y in sorted(fixed):
    d = {k: v for k, v in pyears[y].items() if k != "rings"}
    g_true = fixed[y]                                       # data area (no causeway merge)
    g = unary_union([g_true, close(g_true, 0.6).intersection(CAUSEWAY)]).buffer(0)   # drawing
    # simplify like prep_outlines.py
    polys = []
    for p in parts(g):
        if p.area < MIN_PART: continue
        p = Polygon(p.exterior, [h for h in p.interiors if Polygon(h).area >= MIN_HOLE])
        polys.append(p.simplify(TOL_KM, preserve_topology=True))
    # filled (not imaged) shoreline: new boundary that is not on the imaged boundary
    zone = None
    if y in notes:
        imaged_b = raw[y].boundary.buffer(0.2)
        new_b = unary_union([p.boundary for p in polys])
        z = new_b.difference(imaged_b)
        z = [ln for ln in (z.geoms if hasattr(z, "geoms") else [z]) if ln.length > 0.4]
        zone = unary_union(z).buffer(0.12) if z else None
    fixed_years.append({**d, "rings": [[[r2(v) for v in r.coords.xy[0]], [r2(v) for v in r.coords.xy[1]]]
                                       for p in polys for r in [p.exterior] + list(p.interiors)]})
    rings = []
    for p in polys:
        for j, r in enumerate([p.exterior] + list(p.interiors)):
            if j > 0 and Polygon(r).area < MIN_HOLE_LINE: continue
            rings.append(LineString(r.coords).simplify(SMOOTH_KM))
    if zone is not None:
        main = [r.difference(zone) for r in rings]
        filled = [r.intersection(zone) for r in rings]
    else:
        main, filled = rings, []
    d["x"], d["y"] = line_xy(main)
    px, py = line_xy(filled)
    # fill: every lake part, larger islands cut out
    vx, vy, I, J, K, index = [], [], [], [], [], {}
    for p in polys:
        q = Polygon(p.exterior, [h for h in p.interiors if Polygon(h).area >= MIN_HOLE_FILL])
        q = q.simplify(FILL_SMOOTH, preserve_topology=True).buffer(0)
        for qq in parts(q):
            for tri in constrained_delaunay_triangles(qq).geoms:
                ids = []
                for c in list(tri.exterior.coords)[:3]:
                    k = (r2(c[0]), r2(c[1]))
                    if k not in index: index[k] = len(vx); vx.append(k[0]); vy.append(k[1])
                    ids.append(index[k])
                I.append(ids[0]); J.append(ids[1]); K.append(ids[2])
    d.update(fx=vx, fy=vy, i=I, j=J, k=K)
    if y in notes:
        est = round(g_true.area * ratio, 1)
        d.update(px=px, py=py, filled_note=notes[y], area_measured_km2=d["area_km2"], area_km2=est,
                 area_mi2=round(est / 2.589988, 1))
        a_h = page["reference"]["area_km2"]
        d["vs_ref_km2"] = round(est - a_h, 1)
        d["exposed_km2"] = round(max(0.0, a_h - est), 1); d["exposed_mi2"] = round(d["exposed_km2"] / 2.589988, 1)
        print(y, "measured", d["area_measured_km2"], "-> gap-filled", est, "km2;", len([v for v in px if v is not None]), "filled pts")
    lv = next(l for l in levels if l["year"] == y)
    d["need"] = lv["need"]
    out_years.append(d)
    for k in ("area_km2", "area_mi2", "exposed_km2", "exposed_mi2", "vs_ref_km2", "need"):
        fixed_years[-1][k] = d.get(k)

if fixed_out:
    json.dump({**{k: v for k, v in page.items() if k != "years"}, "years": fixed_years}, open(fixed_out, "w"), separators=(",", ":"))
json.dump({"union_km2": page["union_km2"], "reference": page["reference"], "donors": {str(k): v for k, v in DONORS.items()},
           "levels": levels, "years": out_years}, open(out, "w"), separators=(",", ":"))
print("years", len(out_years), "line pts", sum(len(d["x"]) for d in out_years))
