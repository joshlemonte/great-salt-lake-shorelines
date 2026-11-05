"""Read per-year GSL outlines from the geodatabase and write a compact JSON
for the 3D figure. Usage: python prep_outlines.py <GSL_landsat.gdb> <out.json>"""
import json, sys, re
import pyogrio, geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union

gdb, out = sys.argv[1], sys.argv[2]
TOL = 60          # m, simplification tolerance
MIN_HOLE = 1e6    # m2, islands smaller than this are dropped
MIN_PART = 5e6    # m2, detached lake parts smaller than this are dropped
X0, Y0 = 330000, 4490000   # km origin (UTM 12N)

def year_layers(path):
    """(year, GeoDataFrame) for each fall outline, from the geodatabase
    (one w_fall_YYYY layer per year) or the merged shapefile (label field)."""
    if path.lower().endswith(".shp"):
        allg = gpd.read_file(path)
        for lab, g in allg.groupby("label"):
            m = re.match(r"Fall (\d{4})$", str(lab))
            if m:
                yield int(m.group(1)), g
        return
    for name, _ in pyogrio.list_layers(path):
        m = re.match(r"w_fall_(\d{4})$", name)
        if m:
            yield int(m.group(1)), gpd.read_file(path, layer=name)


years, full = [], []
for yr, g in year_layers(gdb):
    if g.crs is not None and g.crs.to_epsg() != 32612:
        g = g.to_crs(32612)
    a = g.drop(columns="geometry").iloc[0]
    geom = g.geometry.union_all()
    full.append(geom)
    polys = list(geom.geoms) if isinstance(geom, MultiPolygon) else [geom]
    rings = []
    for p in polys:
        if p.area < MIN_PART:
            continue
        p = Polygon(p.exterior, [h for h in p.interiors if Polygon(h).area >= MIN_HOLE])
        p = p.simplify(TOL, preserve_topology=True)
        for r in [p.exterior] + list(p.interiors):
            xs, ys = r.coords.xy
            rings.append([[round((x - X0) / 1000, 2) for x in xs],
                          [round((y - Y0) / 1000, 2) for y in ys]])
    years.append({"year": yr, "area_km2": round(float(a["area_km2"]), 1),
                  "area_mi2": round(float(a["area_mi2"]), 1), "quality": a["quality"],
                  "edge_gap": round(float(a["edge_gap"]), 2), "n_scenes": int(a["n_scenes"]),
                  "sensors": a["sensors"], "rings": rings})
years_unsorted = list(years)
years.sort(key=lambda d: d["year"])
HEALTHY_FT = 4198.0      # minimum healthy lake level (NGVD29 ft)
LEVELS = sys.argv[3] if len(sys.argv) > 3 else None   # fall gage levels CSV
# Exposed lakebed: dry land inside the largest combined extent of the series
# (union of every year's lake). Each year lies inside the union, so
# exposed = union area - lake area. Planar UTM area is scaled to geodesic
# with the median ratio from the per-year areas (differs by <0.1%).
union = unary_union(full)
ratio = sorted(y["area_km2"] / (f.area / 1e6) for y, f in zip(years_unsorted, full))
ratio = ratio[len(ratio) // 2]
union_km2 = union.area / 1e6 * ratio
for y in years:
    y["exposed_km2"] = round(union_km2 - y["area_km2"], 1)
    y["exposed_mi2"] = round(y["exposed_km2"] / 2.589988, 1)
healthy = None
if LEVELS:
    # Fall lake level (Sep 1-Nov 15 mean, Saltair gage, NGVD29) for each year, and
    # the Landsat-equivalent lake area at the healthy minimum: a linear fit of
    # area on level for good/check years within 3 ft of it. Using the outlines'
    # own area-level relation keeps the comparison on one method.
    import csv, numpy as np
    lev = {int(r["year"]): r for r in csv.DictReader(open(LEVELS))}
    pts = [(float(lev[y["year"]]["elev_ngvd29"]), y["area_km2"]) for y in years
           if y["year"] in lev and y["quality"] != "poor"
           and abs(float(lev[y["year"]]["elev_ngvd29"]) - HEALTHY_FT) <= 3]
    c = np.polyfit([p[0] for p in pts], [p[1] for p in pts], 1)
    a_h = float(np.polyval(c, HEALTHY_FT))
    usgs_h = None
    healthy = {"elev_ft": HEALTHY_FT, "area_km2": round(a_h, 1), "fit_n": len(pts),
               "fit_slope_km2_per_ft": round(float(c[0]), 1)}
    for y in years:
        r = lev.get(y["year"])
        if r is None:
            continue
        y["elev_ft"] = round(float(r["elev_ngvd29"]), 2)
        y["dev_ft"] = round(y["elev_ft"] - HEALTHY_FT, 2)
        y["vs_healthy_km2"] = round(y["area_km2"] - a_h, 1)
        # lakebed that is under water at the healthy minimum but dry this fall
        y["exposed_km2"] = round(max(0.0, a_h - y["area_km2"]), 1)
        y["exposed_mi2"] = round(y["exposed_km2"] / 2.589988, 1)
    print(f"healthy minimum {HEALTHY_FT} ft -> {a_h:.0f} km2 "
          f"(n={len(pts)}, {c[0]:.0f} km2/ft)")
json.dump({"origin_utm12n": [X0, Y0], "union_km2": round(union_km2, 1), "healthy": healthy,
           "years": years}, open(out, "w"), separators=(",", ":"))
print("max combined extent", round(union_km2), "km2; geodesic/planar", round(ratio, 5))
print(len(years), "years;", sum(len(r[0]) for d in years for r in d["rings"]), "points")
