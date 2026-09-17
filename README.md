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
(pitch clamped to +/-89 degrees), scroll wheel to zoom, Esc to quit. It's
shaded by an elevation color ramp (dark green low, tan mid, white near the
top of whatever range the current window spans) since no orthophoto
imagery is loaded. The top-right HUD shows your current compass heading
(degrees, 0 = north, clockwise -- plus the 16-point label, N/NNE/NE/.../
NNW) and the current field of view.

The default 45 degree vertical FOV (`FOV_DEFAULT_DEG`) is already narrower
than raylib's usual 70 degree default on purpose: a wide FOV is much wider
than the angle your monitor actually subtends in your real field of view
from a normal desk-viewing distance, so distant features (a mountain miles
off) render proportionally smaller on screen than they actually look
standing there in person, even though every distance/geometry is correct.
Scroll to zoom from `FOV_MIN_DEG` (4, a tight telephoto-like look at
something distant) up to `FOV_MAX_DEG` (90, wide-angle); `TV_FOV` sets the
starting value (mainly for `TV_SCREENSHOT` testing, see below).

The terrain is one regenerable heightfield window, rebuilt (not every
frame) once you've walked far enough from its center or turned more than
~20 degrees from the facing it was last built for. Its shape and extent
are visibility-driven, not a fixed size: `ComputeVisibleDistanceKm` in
`terrain_viewer.c` ray-marches out from the camera along several bearings
across your current facing (plus two more to each side), using the same
Earth-curvature-and-refraction model as `--profile` mode's sightline, to
find how far anything could actually be *seen* in that direction -- capped
by the horizon curvature itself, or by a nearer ridge poking up and
blocking the rest, whichever comes first (hard-capped at 20km regardless,
past which the fixed vertex budget is too coarse to be worth extending
further). The window is built oriented to that facing, wide/deep in front
of you and only a small fixed margin behind, with vertex density and DTM
pyramid level both graded by distance from the camera (dense/full-res
near, sparse/coarse far) -- so how much DTM data actually gets touched
tracks what's genuinely visible: standing in the open Bighorn Basin
looking down its length reaches on the order of 10-15km, while facing
straight into a nearby slope on the Bighorn tile drops to a few hundred
meters. The bottom HUD line reports the current window's computed ahead/
side extents.

DTM pixel data is memory-mapped, not read via ordinary file I/O -- but not
by anything in this project's own code: libtiff's default `TIFFOpen(path,
"r")` already `mmap()`s the whole file at open time (confirmed with
`strace`: one `mmap()` of the entire multi-GB file, zero further `read`/
`pread` calls afterward), so tile decompression already pulls compressed
bytes straight out of that mapping, lazily faulting in pages from disk
only as tiles are actually touched -- the same "mmap once, let the OS
handle residency" idea as `kdtree`'s own shard files, just arrived at via
libtiff's own defaults rather than anything explicit in `dtm.c`.

`TV_LAT`/`TV_LON` override the starting point; `TV_SCREENSHOT` (a filename)
plus `TV_SCREENSHOT_FRAME` (a frame count, default 30) take a screenshot
and exit -- a headless-testing convenience for automated/remote runs.
`TV_YAW`/`TV_PITCH`/`TV_FOV` set the starting look direction and zoom
(mouse-look and scroll-to-zoom are both skipped entirely when
`TV_SCREENSHOT` is set, since a headless run has no real pointer and would
otherwise drift non-deterministically).

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
  multi-ring clipmap -- there's a small pop when it rebuilds. Good enough
  to "stand on the surface and look around"; a clipmap (or caching more
  than one window so a rebuild reuses recently-seen terrain instead of
  resampling it) would be the natural next step if the pop-in becomes
  annoying.
- A distant feature's rendered silhouette can still shift slightly between
  rebuilds -- every rebuild resamples the whole window (including far
  vertices) from scratch, at a coarse DTM pyramid level, through a sparse,
  yaw-rotated vertex grid, so the exact points sampled along, say, a
  faraway ridgeline aren't identical from one rebuild to the next. The
  side-axis warp (symmetric now, matching the forward axis -- see
  `BuildTerrainWindow`) and a larger yaw-triggered rebuild threshold
  (`TERRAIN_REBUILD_YAW_DEG`, 60 rather than 20) cut this down a lot --
  scanning across a peak now shifts it smoothly across frame rather than
  visibly reshaping it -- but a fully stable distant skyline would need a
  separate, far-terrain representation that persists across near-window
  rebuilds instead of being resampled by them (e.g. a cached horizon strip
  built from the same per-bearing ray march `ComputeVisibleDistanceKm`
  already does), which is a bigger change than this pass made.
- The rebuild trigger compares the camera's position against the
  *current* window's own forward/side/behind extents (a fixed fraction of
  each, `TERRAIN_REBUILD_MARGIN`), not a flat distance -- important because
  those extents are themselves visibility-driven and can be much smaller
  than a fixed threshold would assume (e.g. ~500m facing straight into a
  slope). A flat-distance trigger let the camera walk past a small
  window's real edge before ever rebuilding, which looked like "the
  ground disappears and I can see under the map."
- The visibility ray march decides how *far* to extend geometry in a
  direction, not which interior points to skip -- a distant peak poking
  over a closer ridge correctly extends that direction's render distance,
  but the hidden valley in between still gets drawn (it's one continuous
  heightfield mesh, not a masked one). An actual occlusion mask would need
  irregular mesh topology or a ray-marched/voxel renderer.
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
