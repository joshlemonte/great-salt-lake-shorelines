"""
gsl_lake_levels.py

Bathtub-fill outlines of Great Salt Lake at a series of water-surface
elevations (record high to record low), built from a topobathymetric DEM.

For each level the script:
  1. flags every cell at or below the water surface,
  2. keeps only the cells hydraulically connected to "seed" points in the
     lake (so diked ponds, playas and closed depressions don't count as lake),
  3. converts that region to a polygon, removes tiny holes, and records area.
It then merges all levels into one polygon feature class and one outline
(line) feature class, sorted highest to lowest.

Run in the ArcGIS Pro Python window, a notebook, or with propy.bat.
Requires the Spatial Analyst extension.

Recommended DEMs
  * USGS TBDEM (Root, 2023; doi:10.5066/P9DGG75W): 0.5 m, NAVD88 meters,
    2005-06 sonar + 2016 lidar. Best surface; very large (~32 GB), so let the
    script aggregate it (WORK_CELLSIZE). Set DATUM_OFFSET_FT (see below).
  * Tarboton GSLDEM.tif (HydroShare 582060f00f6b443bb26e896426d9f62a):
    NGVD29 meters, extended into the West Desert. Use DATUM_OFFSET_FT = 0.
    Use this one if the USGS DEM's footprint clips the high-stand shoreline.
"""

import os
import arcpy
from arcpy.sa import Con, RegionGroup, Raster, Aggregate

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------
DEM = r"C:\GSL\GSLDEM.tif"                 # topobathymetric DEM
OUT_GDB = r"C:\GSL\GSL_lake_levels.gdb"    # created if missing

DEM_Z_UNITS = "m"        # vertical units of the DEM: "m" or "ft"

# Lake levels below are NGVD29 feet (the datum the USGS gage records use).
# DATUM_OFFSET_FT = NAVD88 minus NGVD29 at the lake, in feet.
#   0.0   if the DEM is NGVD29 (Tarboton GSLDEM)
#   3.48  if the DEM is NAVD88 (USGS TBDEM): the USGS 2023 elevation-area-volume
#         tables list NAVD88 = NGVD29 + 3.48 ft for Great Salt Lake.
DATUM_OFFSET_FT = 0.0

# Working resolution in DEM horizontal units (m). Aggregating the 0.5 m
# USGS DEM to 5-10 m makes this run in minutes instead of hours.
# None = use native resolution.
WORK_CELLSIZE = 10
# MINIMUM keeps the deepest cell (lake connects more easily);
# MAXIMUM preserves thin barriers such as dikes and the causeway;
# MEAN is a neutral compromise.
AGGREGATE_TYPE = "MEAN"

# Mode 1: named levels (NGVD29 ft, label). Edit freely.
LEVELS = [
    (4211.6, "1986-87 record high"),
    (4208.0, ""),
    (4205.0, "Upper bound, healthy range"),
    (4200.0, "Historical average"),
    (4198.0, "Lower bound, healthy range"),
    (4195.0, ""),
    (4191.3, "1963 low"),
    (4188.5, "Nov 2022 record low"),
]

# Mode 2: regular interval. Set USE_STEP = True to ignore LEVELS above.
USE_STEP = False
STEP_TOP_FT, STEP_BOTTOM_FT, STEP_FT = 4212.0, 4188.0, 2.0

# Seed points (lon, lat, WGS84) inside the lake that are wet even at the
# lowest level. Any connected water body touching a seed is kept.
# Add seeds for Farmington Bay / Bear River Bay if you want them included
# when they disconnect from the main lake at low levels.
SEEDS = [
    (-112.55, 41.00),   # south arm (Gilbert Bay), west of Antelope Island
    (-112.70, 41.45),   # north arm (Gunnison Bay), west of Promontory
]

CONNECTIVITY = "FOUR"     # FOUR is stricter about diagonal leaks through dikes
MIN_HOLE_KM2 = 0.05       # fill interior holes smaller than this (0 = keep all)
SIMPLIFY = "SIMPLIFY"     # or "NO_SIMPLIFY" for stair-step raster edges

# ---------------------------------------------------------------------------
FT_TO_M = 1200.0 / 3937.0   # US survey foot (NGVD29 gage elevations)


def level_list():
    if not USE_STEP:
        return sorted(LEVELS, key=lambda t: t[0], reverse=True)
    out, z = [], STEP_TOP_FT
    while z >= STEP_BOTTOM_FT - 1e-9:
        out.append((round(z, 2), ""))
        z -= STEP_FT
    return out


def to_dem_z(level_ngvd29_ft):
    z_ft = level_ngvd29_ft + DATUM_OFFSET_FT
    return z_ft * FT_TO_M if DEM_Z_UNITS.lower().startswith("m") else z_ft


def main():
    arcpy.CheckOutExtension("Spatial")
    arcpy.env.overwriteOutput = True

    folder, gdb_name = os.path.split(OUT_GDB)
    if not arcpy.Exists(OUT_GDB):
        arcpy.management.CreateFileGDB(folder, gdb_name)

    dem = Raster(DEM)
    sr = dem.spatialReference
    arcpy.env.outputCoordinateSystem = sr
    arcpy.env.extent = dem.extent
    arcpy.env.snapRaster = DEM

    # Optional aggregation to a coarser working grid
    if WORK_CELLSIZE:
        factor = int(round(WORK_CELLSIZE / dem.meanCellWidth))
        if factor > 1:
            arcpy.AddMessage(f"Aggregating DEM by {factor}x ({AGGREGATE_TYPE})...")
            agg_path = os.path.join(arcpy.env.scratchGDB, "dem_work")
            Aggregate(dem, factor, AGGREGATE_TYPE, "TRUNCATE", "DATA").save(agg_path)
            dem = Raster(agg_path)
            arcpy.env.snapRaster = agg_path

    # Seeds in DEM coordinates
    wgs84 = arcpy.SpatialReference(4326)
    seed_xy = []
    for lon, lat in SEEDS:
        p = arcpy.PointGeometry(arcpy.Point(lon, lat), wgs84).projectAs(sr)
        seed_xy.append((p.firstPoint.X, p.firstPoint.Y))

    poly_fcs = []
    rg_path = os.path.join(arcpy.env.scratchGDB, "rg_tmp")

    for level_ft, label in level_list():
        z = to_dem_z(level_ft)
        arcpy.AddMessage(f"Level {level_ft:.2f} ft NGVD29  ->  {z:.3f} {DEM_Z_UNITS} in DEM")

        wet = Con(dem <= z, 1)
        RegionGroup(wet, CONNECTIVITY, "WITHIN", "NO_LINK").save(rg_path)

        ids = set()
        for x, y in seed_xy:
            v = arcpy.management.GetCellValue(rg_path, f"{x} {y}").getOutput(0)
            if v.strip().lower() in ("nodata", ""):
                arcpy.AddWarning(f"  Seed ({x:.0f}, {y:.0f}) is dry at this level.")
            else:
                ids.add(int(float(v)))
        if not ids:
            arcpy.AddWarning("  No wet seeds; level skipped.")
            continue

        where = "VALUE IN ({})".format(",".join(str(i) for i in sorted(ids)))
        lake = Con(Raster(rg_path), 1, "", where)

        tag = f"{level_ft:.2f}".replace(".", "_")
        tmp_poly = rf"memory\poly_{tag}"
        tmp_diss = rf"memory\diss_{tag}"
        out_fc = os.path.join(OUT_GDB, f"lvl_{tag}")

        arcpy.conversion.RasterToPolygon(lake, tmp_poly, SIMPLIFY, "Value")
        arcpy.management.Dissolve(tmp_poly, tmp_diss, multi_part="MULTI_PART")
        if MIN_HOLE_KM2 and MIN_HOLE_KM2 > 0:
            arcpy.management.EliminatePolygonPart(
                tmp_diss, out_fc, "AREA", f"{MIN_HOLE_KM2} SquareKilometers",
                part_option="CONTAINED_ONLY")
        else:
            arcpy.management.CopyFeatures(tmp_diss, out_fc)

        for name, ftype, length in [("level_ft_ngvd29", "DOUBLE", None),
                                    ("label", "TEXT", 60),
                                    ("area_km2", "DOUBLE", None),
                                    ("area_mi2", "DOUBLE", None)]:
            arcpy.management.AddField(out_fc, name, ftype, field_length=length)

        with arcpy.da.UpdateCursor(
                out_fc, ["SHAPE@", "level_ft_ngvd29", "label", "area_km2", "area_mi2"]) as cur:
            for row in cur:
                a_km2 = row[0].getArea("GEODESIC", "SQUAREKILOMETERS")
                row[1:] = [level_ft, label, round(a_km2, 2), round(a_km2 / 2.589988, 2)]
                cur.updateRow(row)
                arcpy.AddMessage(f"  area = {a_km2:,.0f} km2 ({a_km2 / 2.589988:,.0f} mi2)")

        poly_fcs.append(out_fc)
        for fc in (tmp_poly, tmp_diss):
            arcpy.management.Delete(fc)

    if not poly_fcs:
        arcpy.AddError("No levels produced output.")
        return

    all_poly = os.path.join(OUT_GDB, "GSL_levels_polygons")
    all_line = os.path.join(OUT_GDB, "GSL_levels_outlines")
    arcpy.management.Merge(poly_fcs, all_poly)
    # Levels are nested and overlap, so treat each polygon independently
    arcpy.management.PolygonToLine(all_poly, all_line, "IGNORE_NEIGHBORS")
    # PolygonToLine keeps only ORIG_FID; join the attributes back
    arcpy.management.JoinField(all_line, "ORIG_FID", all_poly, "OBJECTID",
                               ["level_ft_ngvd29", "label", "area_km2", "area_mi2"])

    arcpy.AddMessage(f"Done. Polygons: {all_poly}\n      Outlines: {all_line}")


if __name__ == "__main__":
    main()
