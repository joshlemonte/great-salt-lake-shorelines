"""Static figure: stacked Great Salt Lake shorelines (3D) + lake area by year.
Usage: python make_static.py outlines.json out.png [azim] [elev]"""
import json, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

src, out = sys.argv[1], sys.argv[2]
AZIM = float(sys.argv[3]) if len(sys.argv) > 3 else -62
ELEV = float(sys.argv[4]) if len(sys.argv) > 4 else 34
Z_RATIO = 0.55        # stack height relative to the lake's east-west width
SMOOTH_KM = 0.15      # drawing simplification (the data are 30 m)
FILL_ALPHA = 0.07     # translucent water fill per layer
MIN_RING_KM2 = 5      # drop islands / holes smaller than this in the drawing

from shapely.geometry import LinearRing, Polygon
d = json.load(open(src))["years"]
for y in d:
    keep = []
    for j, r in enumerate(y["rings"]):
        if len(r[0]) < 4:
            continue
        ring = LinearRing(list(zip(r[0], r[1])))
        if j > 0 and Polygon(ring).area < MIN_RING_KM2:
            continue
        xs, ys = ring.simplify(SMOOTH_KM).coords.xy
        keep.append([list(xs), list(ys)])
    y["rings"] = keep
years = np.array([y["year"] for y in d])
areas = np.array([y["area_km2"] for y in d])
beds = np.array([y["exposed_km2"] for y in d])
H = json.load(open(src))["healthy"]
A_H, H_FT = H["area_km2"], H["elev_ft"]
devs = np.array([y["dev_ft"] for y in d])
DEV_SPAN = 12          # color saturates at ±12 ft from the healthy minimum
LAKE_C, BED_C = "#2a78d6", "#eb6834"

plt.rcParams.update({
    "font.family": ["Liberation Sans", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 10, "axes.edgecolor": "#c9cfcd", "axes.labelcolor": "#3d4547",
    "xtick.color": "#5b6466", "ytick.color": "#5b6466", "text.color": "#15191a",
})
INK, MUTED, RULE, BG = "#15191a", "#5b6466", "#d9dddb", "#ffffff"
# diverging: fall level far below the healthy minimum (orange) -> at it (gray) -> far above (blue)
cmap = LinearSegmentedColormap.from_list(
    "level", ["#9c3a12", "#d9581f", "#f0a07a", "#b3b0a9", "#7fb0ee", "#2a78d6", "#0d366b"])
norm = Normalize(vmin=-DEV_SPAN, vmax=DEV_SPAN)
REF_C = "#5b6466"

imax, imin = int(np.argmax(areas)), int(np.argmin(areas))
label_years = sorted({years[0], years[imax], years[imin], years[-1],
                      *[y for y in (1992, 1999) if y in years]})

fig = plt.figure(figsize=(13, 7.2), dpi=200, facecolor=BG)
ax = fig.add_axes([-0.04, 0.0, 0.70, 0.98], projection="3d")
ax.set_facecolor(BG)
z0 = years[0] - 1
for i, y in enumerate(d):
    segs = [np.column_stack([r[0], r[1], np.full(len(r[0]), y["year"])]) for r in y["rings"]]
    poor = y["quality"] == "poor"
    lw = 1.7 if i in (imax, imin) else 0.75
    lc = Line3DCollection(segs, colors=[cmap(norm(y["dev_ft"]))], linewidths=lw,
                          linestyles=(0, (3, 2)) if poor else "solid", alpha=0.95)
    ax.add_collection3d(lc)
    ext = y["rings"][0]
    fill = Poly3DCollection([np.column_stack([ext[0], ext[1], np.full(len(ext[0]), y["year"])])],
                            facecolors=[cmap(norm(y["dev_ft"]))], edgecolors="none", alpha=FILL_ALPHA)
    ax.add_collection3d(fill)

allx = np.concatenate([r[0] for y in d for r in y["rings"]])
ally = np.concatenate([r[1] for y in d for r in y["rings"]])
ax.set_xlim(allx.min(), allx.max()); ax.set_ylim(ally.min(), ally.max())
ax.set_zlim(z0, years[-1] + 1)
ax.set_box_aspect((np.ptp(allx), np.ptp(ally), Z_RATIO * np.ptp(allx)))
ax.view_init(elev=ELEV, azim=AZIM)
ax.set_xticks([]); ax.set_yticks([])
ymax_, ymin_ = int(years[imax]), int(years[imin])
have = set(years.tolist())
zt = sorted({ymax_, ymin_} | {z for z in range(1990, int(years[-1]) + 1, 5)
                             if z in have and min(abs(z - ymax_), abs(z - ymin_)) >= 4})
if int(years[-1]) - max(zt) >= 3:
    zt.append(int(years[-1]))
zl = [f"{z}  high stand" if z == ymax_ else f"{z}  low" if z == ymin_ else str(z) for z in zt]
ax.set_zticks(zt); ax.set_zticklabels(zl, fontsize=9, color=MUTED)
ax.zaxis.set_tick_params(pad=14)
for a in (ax.xaxis, ax.yaxis, ax.zaxis):
    a.pane.set_facecolor((1, 1, 1, 0)); a.pane.set_edgecolor(RULE)
    a._axinfo["grid"].update(color=RULE, linewidth=0.5)
ax.xaxis.line.set_color((1, 1, 1, 0)); ax.yaxis.line.set_color((1, 1, 1, 0))
ax.zaxis.line.set_color(RULE)

# label selected rings at their easternmost point
ax.text2D(0.05, 0.93, "Great Salt Lake shoreline, each fall", transform=ax.transAxes,
          fontsize=15, fontweight="bold", color=INK)
ax.text2D(0.05, 0.885, f"{years[0]}–{years[-1]}  ·  one ring per year, colored by fall lake level "
          f"relative to the {H_FT:,.0f} ft healthy minimum  ·  Landsat 30 m",
          transform=ax.transAxes, fontsize=9.5, color=MUTED)

# right: area by year
ax2 = fig.add_axes([0.70, 0.39, 0.24, 0.40])
def gapped(v):
    xs, ys = [], []
    for i, y in enumerate(years):
        if i and y - years[i - 1] > 1:
            xs.append(y - 0.5); ys.append(np.nan)
        xs.append(y); ys.append(v[i])
    return xs, ys
MISSING = sorted(set(range(years[0], years[-1] + 1)) - set(years.tolist()))
ax2.plot(*gapped(areas), color=LAKE_C, lw=1.6, zorder=1)
ax2.plot(*gapped(beds), color=BED_C, lw=1.6, zorder=1)
ax2.axhline(A_H, color=REF_C, lw=1, ls=(0, (4, 3)), zorder=0)
ax2.text(years[-1] + 1, A_H + 60, f"Healthy minimum ({H_FT:,.0f} ft)", fontsize=8, color=MUTED,
         ha="right", va="bottom")
ax2.scatter(years, beds, color=BED_C, s=12, zorder=2, edgecolors="none")
ax2.text(years[-1] + 0.6, areas[-1], "Lake", fontsize=8.5, color=INK, va="center")
ax2.text(years[-1] + 0.6, beds[-1], "Exposed\nbelow\nhealthy", fontsize=8.5, color=INK, va="center")
ax2.scatter(years, areas, c=[cmap(norm(v)) for v in devs], s=18, zorder=2,
            edgecolors="none")
for i in (imax, imin):
    hi = i == imax
    ax2.annotate(f"{years[i]}\n{areas[i]:,.0f} km²", (years[i], areas[i]),
                 textcoords="offset points", xytext=(7, -2) if hi else (0, -10),
                 ha="left" if hi else "center", va="top",
                 fontsize=8, color=INK)
poor = [i for i, y in enumerate(d) if y["quality"] == "poor"]
if poor:
    ax2.scatter(years[poor], areas[poor], s=26, facecolors=BG,
                edgecolors=[cmap(norm(devs[i])) for i in poor], linewidths=1.2, zorder=3)
ax2.set_title("Area (km²)", loc="left", fontsize=10.5, fontweight="bold", color=INK)
ax2.spines[["top", "right"]].set_visible(False)
ax2.grid(axis="y", color=RULE, lw=0.5); ax2.set_axisbelow(True)
ax2.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter("{x:,.0f}"))
ax2.set_xlim(years[0] - 1, years[-1] + 1)
ax2.set_ylim(0, np.ceil((areas.max() + 300) / 1000) * 1000)
ax2.tick_params(labelsize=8.5, length=0)

# ring-color key
cax = fig.add_axes([0.70, 0.282, 0.24, 0.014])
cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
cb.outline.set_visible(False)
cb.set_ticks([-DEV_SPAN, -6, 0, 6, DEV_SPAN])
cb.set_ticklabels([f"−{DEV_SPAN} ft", "−6", f"{H_FT:,.0f} ft", "+6", f"+{DEV_SPAN} ft"])
cax.tick_params(labelsize=8, length=0, colors=MUTED)
cax.set_title("Ring color: fall lake level vs. healthy minimum", loc="left", fontsize=8.5,
              color=INK, pad=4)

notes = ("Each outline is the typical shoreline over Sep 1–Nov 15, from all clear Landsat looks.\n"
         "Open circle / dashed ring: cloud gaps along part of the shore (lower confidence).\n"

         f"Healthy minimum {H_FT:,.0f} ft = {A_H:,.0f} km² of lake, from a fit of these Landsat areas\n"
         f"to fall lake level (Saltair gage, Sep 1–Nov 15 mean, NGVD29).\n"
         f"USGS bathymetry gives 3,684 km² at {H_FT:,.0f} ft.\n"
         f"Exposed below healthy level = {A_H:,.0f} km² minus that fall's lake area (0 if larger).\n"
         "Includes Farmington and Bear River bays only while connected to the main lake.\n"
         "1987 lobe at lower left: flooded low ground west of the lake, cut off at\n"
         "the western edge of the area searched."
         + (f"\nNo outline for {', '.join(map(str, MISSING))} (processing failed)." if MISSING else ""))
fig.text(0.70, 0.025, notes, fontsize=7.0, color=MUTED, va="bottom", linespacing=1.45)
fig.savefig(out, dpi=200, facecolor=BG)
print("saved", out)
