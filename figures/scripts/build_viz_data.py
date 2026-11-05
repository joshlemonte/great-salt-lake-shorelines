"""Turn outlines.json into the compact data the interactive page draws:
smoothed rings (for lines) and a triangulated fill per year.
Usage: python build_viz_data.py outlines.json viz_data.json"""
import json, sys
from shapely.geometry import LinearRing, Polygon
from shapely import constrained_delaunay_triangles

SMOOTH_KM, MIN_HOLE_KM2, FILL_SMOOTH_KM = 0.10, 3.0, 0.4
src, out = sys.argv[1], sys.argv[2]
d = json.load(open(src))
years = []
for y in d["years"]:
    rings = []
    for j, r in enumerate(y["rings"]):
        if len(r[0]) < 4:
            continue
        lr = LinearRing(list(zip(r[0], r[1])))
        if j > 0 and Polygon(lr).area < MIN_HOLE_KM2:
            continue
        rings.append(lr.simplify(SMOOTH_KM))
    # one line trace per year: rings joined with null breaks
    xs, ys = [], []
    for lr in rings:
        cx, cy = lr.coords.xy
        xs += [round(v, 2) for v in cx] + [None]
        ys += [round(v, 2) for v in cy] + [None]
    # fill: exterior with the larger islands cut out, triangulated
    poly = Polygon(rings[0], [h for h in rings[1:]
                              if Polygon(h).area >= 20 and Polygon(rings[0]).contains(Polygon(h))])
    poly = poly.simplify(FILL_SMOOTH_KM, preserve_topology=True).buffer(0)
    polys = list(poly.geoms) if poly.geom_type == "MultiPolygon" else [poly]
    vx, vy, I, J, K, index = [], [], [], [], [], {}
    for p in polys:
        tris = constrained_delaunay_triangles(p)
        for t in tris.geoms:
            ids = []
            for c in list(t.exterior.coords)[:3]:
                key = (round(c[0], 2), round(c[1], 2))
                if key not in index:
                    index[key] = len(vx); vx.append(key[0]); vy.append(key[1])
                ids.append(index[key])
            I.append(ids[0]); J.append(ids[1]); K.append(ids[2])
    years.append({k: y.get(k) for k in ("year", "area_km2", "area_mi2", "exposed_km2", "exposed_mi2",
                                         "elev_ft", "dev_ft", "vs_healthy_km2",
                                         "quality", "edge_gap", "n_scenes", "sensors")} |
                 {"x": xs, "y": ys, "fx": vx, "fy": vy, "i": I, "j": J, "k": K})
json.dump({"union_km2": d["union_km2"], "healthy": d.get("healthy"), "years": years}, open(out, "w"), separators=(",", ":"))
print(len(years), "years,",
      sum(len(v["x"]) for v in years), "line pts,",
      sum(len(v["i"]) for v in years), "triangles")
