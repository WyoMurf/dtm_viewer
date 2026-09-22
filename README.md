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

DATA SOURCES:

`terrain_viewer` locates and downloads whatever DTM tile(s) cover wherever
you point it, anywhere in the US -- no manual data prep needed. Three real,
public, unauthenticated data sources, tried in order, best resolution
first (`dtm_fetch.c`):

1. **Wyoming wyolidar** (https://lidar.wygisc.org/wyolidar, bucket
   `wyolidar.s3.arcc.uwyo.edu`) -- ~1m lidar, statewide (40 tiles,
   `W104N041`..`W111N045`). Cloud Optimized GeoTIFFs (BigTIFF), single-band
   16-bit unsigned int elevation in whole meters, EPSG:4326 (geographic
   lat/lon), LZW-compressed, 11-level overview pyramid. This project's
   original three manually-downloaded tiles (Cody/Meeteetse/Bighorn) are
   this same format and are reused as-is, never re-fetched. Tile name is a
   pure formula from lat/lon -- no network round trip needed just to find
   the right one.
2. **USGS 3DEP 1-meter DEM** (also `prd-tnm.s3.amazonaws.com`, under
   `StagedProducts/Elevation/1m/`) -- matches Wyoming's own resolution
   wherever it's been flown, which by now is most of the country. Unlike
   the other two sources, tiles here are organized by named lidar-
   acquisition "Projects" on a UTM grid rather than a clean lat/lon
   formula, so `dtm_fetch.c` queries USGS's TNM Access product-search API
   (`tnmaccess.nationalmap.gov/api/v1/products` -- the same API USGS's own
   download tools use) for the specific tile/URL covering a point, rather
   than computing a filename. Classic TIFF, single-band 32-bit float
   elevation in meters, **UTM projected coordinates** (not geographic
   lat/lon -- confirmed directly by parsing a real tile's GeoKeyDirectory:
   EPSG:26911, NAD83 / UTM zone 11N), NODATA -999999.
3. **USGS 3DEP 1/3 arc-second** (also `prd-tnm.s3.amazonaws.com`) --
   nationwide (~10m) seamless DEM, the final fallback wherever neither of
   the above has coverage. Classic TIFF, single-band 32-bit float
   elevation in meters, clean degree-aligned tiles (`n{lat}w{lon}`
   naming), NODATA -999999. NAD83 vs. Wyoming's WGS84 datum -- a ~1-2m
   horizontal discrepancy in CONUS, well below this tool's precision (also
   true of the 1-meter source's NAD83 UTM tiles).

`dtm.c` reads all three layouts directly (uint16 or float32 samples,
geographic or UTM-projected coordinates, all detected per-file from the
TIFF's own GeoTIFF tags) -- no ingestion tool needed for any of them. UTM
support means parsing the file's GeoKeyDirectoryTag (34735) for
GTModelTypeGeoKey/ProjectedCSTypeGeoKey to detect a projected CRS and
decode its EPSG code into a UTM zone/hemisphere, then a standard
non-iterative Snyder transverse Mercator transform (forward for sampling
a query point against the file's native UTM bounds, inverse only for
the friendly lon/lat range printed when a file opens) -- verified
end-to-end against a real downloaded 1-meter tile: sampling returned real
terrain (2007.4m in the Nevada high desert) instead of silently reporting
NODATA for a point the tile genuinely covers.

At startup, `terrain_viewer` scans two directories for `*.tif` and loads
whatever's already there: `/home/murf/wyodem/lidar` (the original
hand-downloaded Wyoming tiles) and `~/dtm_cache` (where auto-downloaded
tiles, from any of the three sources, land and get reused on future
runs). To add
your own pre-downloaded tiles, just drop them in either directory --
edit `kDtmSearchDirs` in `terrain_viewer.c` to add more.

Once you point the tool at an actual location (walk mode's start point;
`--profile`'s two endpoints plus points every ~15% along the path;
`--viewshed`'s center plus points along 8 compass bearings at
`maxDistanceKm`), it checks whether a loaded file already covers each
point and, if not, looks for the right tile at each source in turn: a
HEAD request confirms it exists and reports its size, then -- unless
`TV_AUTO_DOWNLOAD=1` is set -- prompts `y/N` on stdin before actually
downloading (`curl`) into `~/dtm_cache`. Coverage is checked once at each
mode's natural anchor points, not continuously, so walking far enough out
of a walk-mode session's initial coverage during a long walk can still run
off the edge of loaded data (matching this project's existing "known
limitations, not oversights" pattern -- restart at the new location rather
than expecting it to fetch mid-walk). A point with no coverage at either
source (ocean, outside the US) just prints a warning and shows as NODATA.

`TV_AUTO_DOWNLOAD=1` skips the confirmation prompt -- needed for
`TV_SCREENSHOT` headless/scripted runs, or anyone who wants this
non-interactive.

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

`make dtm_fetch_test` builds a standalone smoke test for the auto-download
logic (no raylib needed): `./dtm_fetch_test cacheDir lat lon` ensures
coverage at that point (downloading if needed, prompting unless
`TV_AUTO_DOWNLOAD=1`) and prints the sampled elevation there.

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

    ./terrain_viewer --profile lat1 lon1 [+height1] lat2 lon2 [+height2] [screenshot.png]

`+height` is an optional mast/tower height above ground at that
endpoint -- e.g. a cell tower's antenna height, not a person standing
there -- defaulting to 0 (ground level) when omitted. It must start with
`+` (a bare number there is parsed as the next lat/lon instead), which
also means it can't disambiguate a `+height` from a positive
(eastern-hemisphere) longitude typed with an explicit `+` -- not a
concern for anything in Wyoming, where longitude is always negative.
Add an `m` suffix for meters (`+30m`) or `f` for feet (`+100f`, case-
insensitive either way); no suffix defaults to feet (`+100` == `+100f`),
matching this project's original convention from before unit suffixes
existed.
Endpoint markers get a short vertical mast line up from the ground when a
height is given, and the "A"/"B" labels show ground/tower/total elevation
separately (e.g. `1524m ground + 100ft = 1554m`) -- the tower height gets
added on top of the ground elevation at that point for every part of the
curvature/refraction sightline math and the endpoint markers, not just
display; raising a tower's height on one or both ends is exactly how you'd
check whether it clears a hill that blocks line-of-sight at ground level.

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
two points. Points with no DTM coverage at all (ocean, outside the US --
see DATA SOURCES above) show as a gap in the terrain line (NODATA).

RUNNING -- viewshed mode:

    ./terrain_viewer --viewshed lat lon [+height] maxDistanceKm [options] [screenshot.png]

    Options (any order, after maxDistanceKm):
      --erp watts    transmitter ERP in watts -- ALSO enables the RF
                     signal-strength model below (default: off, pure
                     geometric LOS map)
      --freq MHz     carrier frequency in MHz (default: 869, cellular Band A)
      --sens dBm     receiver sensitivity threshold (default: -100)

Given a tower's position and antenna height (`+height`, same `+`/unit-
suffix convention as `--profile`, 0 if omitted), rasterizes a
top-down map (north up, tower at center) of coverage within
`maxDistanceKm`. Distance rings mark 25/50/75/100% of `maxDistanceKm`.

Without `--erp` (the default), this is a **pure geometric line-of-sight
map**: green if a straight line from the tower clears the terrain
(curvature and refraction accounted for), red if a nearer ridge blocks it,
gray where the loaded DTM tiles don't reach. No path loss, antenna
radiation pattern, Fresnel-zone clearance, or signal strength is modeled --
green means "a straight line to the tower isn't blocked by terrain", not
"you'll get a strong signal there".

With `--erp`, a real (if simplified) **RF signal-strength model** replaces
the binary map with a continuous one: free-space path loss (`20*log10(d_km)
+ 20*log10(f_MHz) + 32.44`) plus single-knife-edge diffraction loss
(ITU-R P.526 -- same family as `--profile`'s curvature formula) computed
from the single worst obstruction found so far along each ray, relative to
the signal's own Fresnel-zone geometry at that point. Color is a red
(weak, right at the sensitivity threshold) -> yellow -> green (strong)
ramp; gray means the point doesn't reach the sensitivity threshold at all.
The green end of the ramp auto-scales to whatever the strongest margin
actually computed for that run turns out to be (`ComputeViewshedRaster`
does the full sweep into a margin buffer first, then colors it against
that discovered ceiling) rather than a fixed dB span -- a flat 0-40dB
range washes out to solid green the moment a scenario has more than 40dB
of headroom anywhere, which is most of the map for a strong transmitter
at short range (e.g. 140W at 1.5km leaves ~59dB of margin on a clear
path). The color bar's own printed endpoints always show the true dBm
range in play for that run. `--erp` is watts (matching how FCC license records
usually report it, as ERP referenced to a dipole); it's converted to EIRP
internally (`+2.15 dB`) since that's what the path-loss math needs. This
tends to show meaningfully more coverage than the plain LOS map, because
real signals bend around obstacles they can't see past in a straight
line -- particularly noticeable at low cellular bands like 850MHz
(`--freq`'s default), which is exactly why carriers favor low-band
spectrum for rural/terrain-heavy coverage.

What this does **not** model, even with `--erp`: multiple obstructions
along one path (only the single dominant one, unlike full multi-edge
models such as Longley-Rice/ITM -- the FCC's own standard, and a much
bigger undertaking were it ever worth porting), antenna radiation pattern
(isotropic assumed), ground conductivity/reflection, foliage or building
clutter, or troposcatter. Treat this as a solid, physically-grounded
estimate for scouting -- not a substitute for a real RF site survey.

Each map point is tested at `VIEWSHED_RECEIVER_HEIGHT_M` (1.5m, a typical
handset/vehicle height) above the DTM's bare-ground elevation there, while
the terrain itself (not that same receiver height) is what determines
whether a closer point blocks a farther one along each ray -- standard
"observer offset / target offset" viewshed practice, and the reason a
point's own ground elevation can't disqualify itself while still
correctly blocking everything behind it. `ComputeViewshedRaster` in
`terrain_viewer.c` reuses the exact curvature/refraction model `--profile`
does, generalized from "how far can I see in one direction" (a single
running-max-angle horizon march) to a full radial sweep across ~360-3600
bearings (adaptive to `maxDistanceKm` and the fixed 900px raster
resolution) that rasterizes every sample as it's computed, rather than
just the farthest visible point. Runs in under 2 seconds even at a 20km
radius against the full DTM set.

Real streets are overlaid on top of the coverage raster for whatever
community happens to fall inside the mapped circle -- `osm_roads.c` fetches
every `highway=*` way in a padded bounding box around the map from
OpenStreetMap's public API (`api.openstreetmap.org/api/0.6/map`, no auth,
no key), hand-parses the XML (two linear scans: collect every referenced
node, then resolve each way's `<nd ref>` list against them -- not a general
XML parser, matching how `dtm.c` already hand-parses GeoTIFF tags instead
of depending on more of libtiff than it needs), and draws each way as a
thin line in the same tower-centered, north-up local-meters projection
`ComputeViewshedRaster` used for the coverage colors. Best-effort and never
fatal: a fetch failure or zero mapped roads just mean no streets get drawn
that run, not a crash. `TV_NO_ROADS=1` skips the fetch entirely (headless/
offline testing).

OSM's live API is the production editing backend, not a CDN built for
heavy polling -- it can hand back a brief failure (rate limiting, a
momentary hiccup) that moments later just works again with no change at
all (observed directly: a real run failed outright, and the exact same
URL succeeded seconds later). `FetchUrl` retries a failed fetch up to
three times with a short backoff before actually giving up, so a
transient blip costs a couple extra seconds rather than losing the whole
overlay for that run.

OSM's `/map` API refuses anything over 0.25 square degrees (confirmed
directly: a real request that size gets a real HTTP 400, "The maximum
bbox size is 0.250000") -- roughly a 20-22km viewshed radius at Wyoming's
latitude. Past that, the request is split into an N x N grid of sub-boxes
that each stay comfortably under the limit, fetched and merged into one
combined result; a way whose nodes straddle a tile boundary comes back in
full from every tile it touches, so `osm_roads.c` tracks already-added way
ids across tiles to avoid drawing (and storing) the same road several
times over. Capped at a 6x6 grid (~80km radius) so a huge radius can't
fire off dozens of sequential requests -- a few failed tiles among many
just means partial coverage near the edges, not a lost overlay.

SCOUTING TOOLS (`make probe_points los_check`):

Two small headless helpers (no raylib, no window) for batch terrain
scouting -- e.g. searching many candidate sites at once for something like
a passive RF reflector location, rather than checking points one at a time
in `--profile`. Both read lines from stdin and print one result line per
input line, so they're meant to be driven by a generated grid of points
(a short script looping over `geo_utils.h`'s `destination_point` math, or
any other source of lat/lon pairs) rather than typed by hand.

    ./probe_points file.tif [file.tif ...]

Reads `lat,lon` per line, prints `lat,lon,elevM` (or `DTM_NODATA`'s
sentinel value if outside coverage). Just a bare elevation lookup -- useful
for finding local high points or low saddles in a region before checking
which ones actually have useful sightlines.

    ./los_check file.tif [file.tif ...]

Reads `lat1,lon1,h1ft,lat2,lon2,h2ft` per line (h in feet above ground at
each point -- always feet here, no `m`/`f` suffix like `--profile`'s
`+height` args take, since this is a plain CSV batch format), prints
`lat1,lon1,lat2,lon2,CLEAR,marginM` or `...,BLOCKED,marginM,atKm` --
reuses the exact same curvature/refraction sightline math `--profile`
plots (worst point along the path relative to the direct line, not just a
single-obstruction check), just without opening a window per pair, so
hundreds of candidate pairs can be tested in well under a second. This is
how the two-tower/five-point Meeteetse reflector-site search was actually
done: generate a grid of candidate points, batch-test all of them against
the towers and the target addresses with `los_check`, then confirm the
winner visually with a real `--profile` screenshot before trusting it.

CASE STUDY -- siting a passive reflector for Meeteetse (KNKN312/KNKN244):

Worked example of using this tool end to end for a real passive-reflector
siting problem, kept here because the method (not just the specific
answer) is reusable for the next one.

**The problem.** Two real cell sites, KNKN312 (44.14206,-108.824333) and
KNKN244 (44.143611,-108.822222), sit just east of Meeteetse; seven
specific residential addresses on the west side of town -- identifying
which houses currently lack coverage, so deliberately not named here even
though this repo is private -- don't have a clear line from either tower,
per `--profile`/`los_check` against the loaded Wyoming lidar tiles. (An
eighth candidate address was dropped once confirmed to already have
reception.) The addresses are referred to below as Target 1-7.

**Site search.** A grid of candidate points (`destination_point` fanned
out at various bearings/distances from the address cluster, piped through
`los_check`) found plenty of points with geometric line-of-sight to both
the towers and the addresses -- but LOS alone isn't enough. The first
"perfect LOS" candidate sat almost exactly on the straight line between a
tower and two of the addresses (included angle ~175-179 deg), which is
useless for a flat mirror: at that near-zero deviation a reflector would
need to operate at grazing incidence (see the sizing formula below --
effective area scales with cos(theta), theta = deviation/2, and cos(90
deg) = 0). The fix was re-scanning for sites *off* the direct tower-
address line, filtering for a genuine bend angle (roughly 60-130 deg is
comfortable), not just raw visibility. Testing actual public-facility
addresses (fire district, sheriff's office -- named, known lease-free
land along the through-town highway corridor) turned out to beat every
private-parcel candidate found by blind grid search, both for siting ease
and for reflector size.

**Sizing a flat-plate reflector.** For two path legs of length d1, d2 (km)
via a flat plate of physical area A (m^2) at incidence angle theta from
the plate's normal (effective area Ae = A*cos(theta)), the isolated-
antenna path loss through the reflector is:

    L_dB = 142.0 + 20*log10(d1) + 20*log10(d2) - 20*log10(Ae)

(derived from treating the plate as an aperture that both intercepts and
re-radiates -- notably frequency-independent, which is the standard,
slightly counterintuitive result for flat passive reflectors: loss
depends on plate area, path lengths, and incidence angle, not on
wavelength.) Solving for the minimum Ae that closes the link against a
transmitter's EIRP and a receiver's sensitivity threshold, then A =
Ae/cos(theta), gives the required physical plate size -- this is exactly
what turned "the same site is either <1m^2 or basically infinite depending
on where you put it" from a vague warning into a concrete go/no-go check
per candidate.

**Result.** Best site found: the Park County Sheriff's Office lot
(44.157895,-108.870449), on a 40ft mast. Confirmed clear (fresh
`los_check` run, not just distance/elevation estimates) to both towers
(4.08-4.16km, 6+m margin) and all 7 remaining addresses (318-832m,
2.7-3.0m margin), with bend angles from 75-129 deg -- comfortably clear
of the grazing problem. Required plate size peaks at 1.23 m^2 (Target 3,
the longest/worst-angled leg); every other address needs less.

A backup was checked the same way ~250m down the same street: the
Meeteetse Fire District lot (44.157664,-108.870856) also clears every
target, with bend angles 69-130 deg and a slightly *smaller* worst-case
plate (1.16 m^2, still Target 3), but needs a taller 50ft mast to get
there (three of the seven addresses are marginally blocked at 40ft, same
as the Sheriff's site's own initial 30ft check). RF-wise the two sites
are essentially interchangeable -- pick whichever turns out easier on the
ground (pole placement, existing infrastructure to tie into, county
sign-off).

**Construction.** At 869 MHz (cellular Band A, lambda ~345mm), flatness
only needs to hold to about lambda/16-lambda/20 (~17-22mm) -- the
Rayleigh criterion used for reflector-antenna surfaces generally, and
comfortably loose for a ~1m panel built with ordinary sheet-metal
fabrication. No optical-grade polish needed; "smooth" here means smooth
relative to a 345mm wavelength, not visible light. A mesh (not a solid
plate) works identically as long as the openings stay under roughly
lambda/10 (~35mm); this cuts wind loading substantially for a panel this
size on a 40-80ft mast. Standard stainless steel window screen (~1.1-
1.4mm openings) is 15-30x finer than that threshold and, despite being a
noticeably worse conductor than aluminum (~1.45e6 vs ~3.5e7 S/m), still
reflects with negligible extra loss at this frequency -- the same reason
stainless mesh is standard material for RF shielding enclosures generally.

**Orientation.** By the law of reflection, the plate's face (its outward
normal) should bisect the bearing to the transmitter and the bearing to
the target, for each address individually; since a panel this size (a
few wavelengths across at 869MHz) has a fairly wide natural beam, one
compromise orientation covering the whole spread of target bearings is
sufficient rather than needing per-address aiming. For the Sheriff's
Office site: tower bearing ~114 deg, target bearings 188-242 deg, giving
individual ideal facings of 151-178 deg and an overall compromise of
**160 deg true**. Bearings from this tool (`geo_utils.h`) are relative to
true north; converting to a compass heading needs the local magnetic
declination (9.92 deg E at this exact site as of September 2026, per
NOAA/NCEI's WMM-2025 calculator -- https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml,
drifting about -0.09 deg/year here), giving **~150 deg on a hand
compass**. The Fire District backup works out nearly identically --
tower bearing ~114 deg, target bearings 148-178 deg, compromise **159 deg
true**, and the same 9.92 deg E declination (only 250m away) gives **~149
deg magnetic** -- unsurprising given how close the two sites are to each
other. Verify with a real signal meter at a few target addresses once
mounted, whichever site is used -- this estimate is a starting point for
the physical install, not a substitute for field tuning.

**Not covered by this tool**: structural/lightning-protection grounding
for the mast and panel (NEC Article 810 territory, especially relevant
here since the site is a public-safety facility) -- that belongs with a
licensed electrician or tower installer, not a terrain/RF estimate.

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
