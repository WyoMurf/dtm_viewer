# DTM Terrain Viewer

A raylib/C tool for exploring real DTM (digital terrain model) elevation
data: `terrain_viewer` stands a first-person camera directly on the
terrain and lets you walk around, or plots an elevation profile between
two lat/lon points.

(Originally built alongside github.com/WyoMurf/kdtree, whose `earth_viewer`
provides a matching whole-globe orbit-camera view with a real Earth
texture and city dots -- this project is standalone and doesn't depend on
it, but the two pair well if you want both a globe view and a ground-level
view of the same place.)

DATA:

This was built against three wyolidar (Wyoming lidar, https://lidar.wygisc.org/wyolidar)
DTM GeoTIFFs covering Cody, Meeteetse, and the Bighorn Mountains -- Cloud
Optimized GeoTIFFs (BigTIFF), single-band 16-bit unsigned int elevation in
whole meters, EPSG:4326 (geographic lat/lon, ~1m/pixel), LZW-compressed,
with an 11-level overview pyramid. No ingestion tool is needed -- `dtm.c`
reads these files directly.

By default `terrain_viewer` looks for these three specific files at fixed
paths on the machine it was developed on (100.74.88.66):

    /home/murf/wyodem/lidar/W109N044_Deg_Cog.tif  (Cody)
    /home/murf/wyodem/lidar/W108N044_Deg_Cog.tif  (Meeteetse)
    /home/murf/wyodem/lidar/W107N044_Deg_Cog.tif  (Bighorn Mountains / Cloud Peak)

To point it at different DTM files, edit `kDefaultDtmPaths` in
`terrain_viewer.c`.

BUILDING:

Needs a C11 compiler, `libtiff-devel` (headers), and for `terrain_viewer`
(not `dtm_test`) a locally-built raylib -- see raylib.com/ for build
instructions; `Makefile`'s `RAYLIB_DIR` points at where it expects to find
`libraylib.a`/`.so` and `raylib.h`/`rlgl.h`/`raymath.h`.

    sudo dnf install libtiff-devel   # (or your distro's equivalent)
    make dtm_test        # standalone sanity check, no raylib needed
    make terrain_viewer   # the raylib walkthrough + profile tool

`dtm_test` is worth running first against your own DTM files -- it prints
elevation at a handful of hardcoded probe points plus a coarse per-file
min/max scan, both useful for confirming a new file's tags parsed correctly
before trusting it in the viewer.

RUNNING -- walkthrough mode:

    ./terrain_viewer [cody|meeteetse|bighorn]

Defaults to `cody`. WASD to move, Shift to sprint, mouse to look around
(pitch clamped to +/-89 degrees), Esc to quit. The terrain is a single
~1.8km-square heightfield window centered on you, rebuilt (not every frame)
once you've walked far enough from its center -- see KNOWN LIMITATIONS
below. It's shaded by an elevation color ramp (dark green low, tan mid,
white near the top of whatever range the current window spans) since no
orthophoto imagery is loaded.

`TV_LAT`/`TV_LON` override the starting point; `TV_SCREENSHOT` (a filename)
plus `TV_SCREENSHOT_FRAME` (a frame count, default 30) take a screenshot
and exit -- a headless-testing convenience for automated/remote runs.
`TV_YAW`/`TV_PITCH` set the starting look direction (mouse-look is skipped
entirely when `TV_SCREENSHOT` is set, since a headless run has no real
pointer and would otherwise drift non-deterministically).

RUNNING -- elevation profile mode:

    ./terrain_viewer --profile lat1 lon1 lat2 lon2 [screenshot.png]

Plots elevation vs. distance between the two points (great-circle path, via
`geo_utils.h`'s `haversine_distance`/`initial_bearing`/`destination_point`),
with a curved "sightline" reference (thin red) representing the actual
straight 3D line of sight between the two endpoints once Earth's curvature
and standard 4/3-radius atmospheric refraction are accounted for --
textbook ITU-R P.526-style: a straight line between two points d1 and d2
away from each end (d1+d2 = total distance) clears the curved Earth's
surface by `d1*d2 / (2*k*R)` less than flat-plane geometry suggests. Any
stretch of raw terrain (thick brown) that pokes above that reference line
is highlighted in bright red -- it blocks direct line-of-sight between the
two points. Points outside all three loaded DTM tiles' coverage show as a
gap in the terrain line (NODATA), which is expected for a long profile that
leaves Wyoming.

KNOWN LIMITATIONS (deliberate v1 scope, not oversights):

- The walkthrough is one regenerable heightfield window, not a proper
  multi-ring clipmap -- there's a small pop when it rebuilds, and no
  distant low-res terrain silhouette beyond the window's edge (you see sky
  past ~900m). Good enough to "stand on the surface and look around"; a
  clipmap would be the natural next step if the pop-in becomes annoying.
- No orthophoto texture (elevation-ramp shading only), except that
  `hillshade_W109N44_Deg_Cog.tif` (already present alongside the Cody DTM
  on the development machine) is a candidate texture for that one tile
  specifically, not yet wired up.
- Elevation values are whole meters in the source data -- a short, nearly
  flat profile can show a visible staircase pattern at zoomed-in scale.
  That's the source data's own precision, not a bug in the sampling code.
- `DtmSampleMeters`'s "first file in the set that covers this point wins"
  rule assumes the DTM tiles you add don't overlap, which the wyolidar
  1-degree tiles don't.
