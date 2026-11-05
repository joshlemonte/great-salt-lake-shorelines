# Great Salt Lake shorelines, 1985–2026

Fall shorelines of Great Salt Lake mapped from Landsat imagery, stacked into a
3D figure and compared with the lake's 4,198 ft minimum healthy level.

![Stacked fall shorelines of Great Salt Lake, 1985–2026](figures/output/GSL_shoreline_stack_1985-2026.png)

**[Open the interactive, rotatable version](https://joshlemonte.github.io/great-salt-lake-shorelines/)** in your browser
(no download needed). The source file is
[`figures/output/gsl_shorelines_interactive.html`](figures/output/gsl_shorelines_interactive.html).

## What is here

| Path | Contents |
|---|---|
| `arcgis/gsl_landsat_outlines_arcpy.py` | ArcGIS Pro script that maps one fall shoreline per year from Landsat |
| `arcgis/gsl_dem_bathtub_levels.py` | Earlier alternative: "bathtub" outlines at set lake levels from a lakebed DEM |
| `data/outlines/GSL_Landsat_outlines.shp` | Shoreline polygons, one per year (UTM 12N, WGS 84) |
| `data/outlines/GSL_Landsat_scenes.csv` | Landsat scenes used for each year, with cloud cover |
| `data/fall_levels_and_usgs_areas.csv` | Fall lake level at the Saltair gage and the USGS bathymetric area at that level |
| `figures/scripts/` | Python that builds the static PNG and the interactive page from the shapefile |
| `figures/output/` | The finished figures |

## Method

**Imagery.** Landsat Collection 2 Level-2 surface reflectance (Landsat 4, 5, 7, 8, 9)
from the Esri Living Atlas "Landsat Level-2" imagery layer. It needs an ArcGIS
Online organizational account and does not use credits.

**One outline per fall (September 1 to November 15).**
1. Scenes over the lake are ordered from least to most cloudy. Landsat 7 is used only
   in 2012, or to fill gaps, because of its scan-line stripes.
2. Cloud, cloud shadow, cirrus and fill pixels are masked from the Pixel QA band.
3. MNDWI = (green − SWIR1) / (green + SWIR1) is computed, and water is classified
   with an Otsu threshold set separately for each scene.
4. A pixel counts as water if it was water in at least half of its clear looks.
5. Only water connected to seed points in the south arm (Gilbert Bay) and north arm
   (Gunnison Bay) is kept. Interior holes that were never observed are filled, so
   real islands stay as land.

**Quality fields** on each outline:
- `edge_gap`: percent of the shoreline that borders unobserved (cloud) pixels.
- `quality`: good (≤ 2%), check (≤ 10%) or poor (> 10%).

**Healthy minimum.** The 4,198 ft (NGVD29) level is about 3,836 km² of lake on this
Landsat basis. That number comes from a straight-line fit of outline area against
fall lake level for the 13 years within 3 ft of 4,198 ft (about 218 km² per ft).
USGS bathymetry gives 3,684 km² at 4,198 ft. "Lakebed exposed below the healthy
level" is 3,836 km² minus that fall's lake area, or zero when the lake is larger.

**Check against USGS.** Each fall's Landsat area is within 0.93–1.14 of the USGS
bathymetric area at that fall's mean Saltair gage level (r = 0.986). The Landsat
areas run about 5% higher on average, consistent with wet margins and connected
ponds being counted as water.

## Known limitations

- **2025** has no outline because processing failed in the last run.
- **Low-confidence years:** 1993, 2009, 2012 and 2018 have cloud gaps along more
  than 10% of the shoreline.
- **Bays:** Farmington and Bear River bays are included only while connected to
  the main lake.
- **1987:** the outline is cut off where flooded ground west of the lake ran past
  the western edge of the area searched (longitude −113.35°).
- **Typical shoreline:** each outline is the typical shoreline for the fall
  window, not a single date.

## Rebuilding the figures

Requires Python with geopandas, pyogrio, shapely ≥ 2.1, numpy and matplotlib.

```bash
cd figures/scripts
python prep_outlines.py ../../data/outlines/GSL_Landsat_outlines.shp outlines.json ../../data/fall_levels_and_usgs_areas.csv
python build_viz_data.py outlines.json viz_data.json
python build_page.py template.html viz_data.json gsl_shorelines.html
python make_static.py outlines.json GSL_shoreline_stack.png
```

## Data sources

- Landsat Collection 2 Level-2 imagery courtesy of the U.S. Geological Survey, via the
  [Esri Living Atlas Landsat Level-2 layer](https://www.arcgis.com/home/item.html?id=bd6b545b95654d91a0b7faf7b5e010f5).
- Lake level: USGS gage 10010000, Great Salt Lake at Saltair Boat Harbor (NGVD29).
- Elevation–area tables: [USGS Great Salt Lake elevation-area-volume tables, 2023](https://www.sciencebase.gov/catalog/item/6467b42fd34ec11ae4a8afb1).
