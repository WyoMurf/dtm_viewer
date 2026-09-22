#define _POSIX_C_SOURCE 200809L /* glob, mkdir */
#define _DEFAULT_SOURCE 1 /* M_PI and friends in math.h, which _POSIX_C_SOURCE alone hides */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <glob.h>
#include <sys/stat.h>

#include "raylib.h"
#include "rlgl.h"
#include "raymath.h"
#include "dtm.h"
#include "dtm_fetch.h"
#include "osm_roads.h"
#include "geo_utils.h"
#include "longley_rice.h"

static void EnsureCoverageOrWarn(double lat, double lon);

/*
 * Ground-level companion to earth_viewer.c: instead of an orbit camera
 * looking at the whole globe, this stands a first-person camera directly
 * on one of the wyolidar DTM tiles (see README-terrain.md) and lets you
 * walk around on the actual terrain, or (--profile mode) plots an
 * elevation-vs-distance graph between two lat/lon points with Earth's
 * curvature accounted for.
 *
 * Coordinates: world space is a local flat-earth tangent plane centered on
 * wherever the walk started (DtmLonLatToLocalMeters/DtmLocalMetersToLonLat,
 * see dtm.h) -- +X is east, +Z is south, so -Z (raylib's default forward)
 * points north, matching earth_viewer.c's LonLatToCartesian convention of
 * "north is the direction you'd expect on a map". World Y is elevation in
 * meters minus the walk's starting elevation (g_originElevation), just to
 * keep the numbers small/centered near the camera instead of ~1500-4000.
 */

typedef struct { const char *name; double lat, lon; } NamedPlace;
static const NamedPlace kPlaces[] = {
    { "cody",      44.52634, -109.05653 },
    { "meeteetse", 44.14670, -108.87100 },
    { "bighorn",   44.38530, -107.16830 }, /* Cloud Peak area */
};
#define PLACE_COUNT (int)(sizeof(kPlaces) / sizeof(kPlaces[0]))

/* Where DTM tiles are looked for at startup (whatever's already there
 * loads immediately, no network needed) and where EnsureDtmCoverage (see
 * dtm_fetch.h) downloads new ones on demand once a target location is
 * known. The first entry is this project's original manually-downloaded
 * Wyoming tiles; the second (computed at startup into g_dtmCacheDir) is
 * where auto-downloaded tiles -- Wyoming or USGS 3DEP -- land and get
 * reused on future runs. See README.md's DATA SOURCES section. */
static const char *kDtmSearchDirs[] = {
    "/home/murf/wyodem/lidar",
};
#define DTM_SEARCH_DIR_COUNT (int)(sizeof(kDtmSearchDirs) / sizeof(kDtmSearchDirs[0]))

static char g_dtmCacheDir[600];
static DtmSet g_dtm = { 0 };

#define EYE_HEIGHT_M 1.7f
#define WALK_SPEED_MPS 4.0f
#define SPRINT_MULT 3.0f
#define MOUSE_SENSITIVITY 0.12f

/* Default vertical FOV, and the scroll-wheel zoom range around it. 70 deg
 * (raylib/many games' usual default) is much wider than the angle a
 * monitor actually subtends in your real field of view from a normal
 * viewing distance -- rendering at that FOV makes anything far away (a
 * mountain miles off) look proportionally smaller on screen than it
 * really looks standing there in person, even though the geometry/
 * distances are all correct. 45 deg is a closer match for a typical
 * "sit at a desk" monitor distance; scrolling zooms further in for a
 * telephoto-like close look at something distant (e.g. Heart Mountain). */
#define FOV_DEFAULT_DEG 45.0f
#define FOV_MIN_DEG 4.0f
#define FOV_MAX_DEG 90.0f
#define FOV_ZOOM_SPEED_DEG 3.0f

/* Smallest signed difference b-a, in degrees, wrapped to (-180, 180] --
 * for comparing two yaw angles that both accumulate unbounded (mouse-look
 * just adds/subtracts degrees with no wraparound) without a 359-vs-1
 * comparison reading as a huge turn. */
static float AngleDiffDeg(float a, float b) {
    float d = fmodf(b - a + 180.0f, 360.0f);
    if (d < 0) d += 360.0f;
    return d - 180.0f;
}

/* wc.yaw (see WalkForward) already IS a compass bearing -- 0 at north,
 * increasing clockwise through east -- it just accumulates unbounded
 * (mouse-look adds/subtracts degrees with no wraparound), so this only
 * needs to fold it into [0, 360). */
static float CompassHeadingDeg(float yawDeg) {
    float h = fmodf(yawDeg, 360.0f);
    if (h < 0) h += 360.0f;
    return h;
}

/* Standard 16-point compass rose abbreviation for a heading already
 * normalized to [0, 360) by CompassHeadingDeg. */
static const char *CompassLabel(float headingDeg) {
    static const char *labels[16] = {
        "N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
        "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"
    };
    int idx = (int)((headingDeg + 11.25f) / 22.5f) % 16;
    return labels[idx];
}

/* --- Gribb-Hartmann frustum extraction, verbatim from earth_viewer.c/viewer.c --- */
typedef struct { float a, b, c, d; } Plane;
static void ExtractFrustumPlanes(Matrix m, Plane out[6]) {
    out[0] = (Plane){ m.m3+m.m0, m.m7+m.m4, m.m11+m.m8,  m.m15+m.m12 };
    out[1] = (Plane){ m.m3-m.m0, m.m7-m.m4, m.m11-m.m8,  m.m15-m.m12 };
    out[2] = (Plane){ m.m3+m.m1, m.m7+m.m5, m.m11+m.m9,  m.m15+m.m13 };
    out[3] = (Plane){ m.m3-m.m1, m.m7-m.m5, m.m11-m.m9,  m.m15-m.m13 };
    out[4] = (Plane){ m.m3+m.m2, m.m7+m.m6, m.m11+m.m10, m.m15+m.m14 };
    out[5] = (Plane){ m.m3-m.m2, m.m7-m.m6, m.m11-m.m10, m.m15-m.m14 };
}

/* --- Terrain window: one regenerable heightfield mesh centered on wherever
 * the camera last was when it was (re)built. Not a multi-ring clipmap --
 * see README-terrain.md's known-limitations note -- rebuilt whenever the
 * camera has wandered far enough from its center, OR turned far enough
 * from the facing direction it was last built for (see below: the
 * window's shape now depends on which way you're looking, not just where
 * you are). --- */
#define TERRAIN_RESOLUTION 254       /* 255x255 vertices = 65,025, just under raylib's 65536-per-mesh index limit */
#define TERRAIN_REBUILD_MARGIN 0.55f /* rebuild once the camera has covered this fraction of the window's own extent in some direction -- see the rebuild-trigger comment in RunWalkMode for why a fixed distance isn't safe here */
#define TERRAIN_REBUILD_YAW_DEG 60.0f /* ...or turned this many degrees from the window's build-time facing. Larger than you might expect: a full rebuild resamples the whole window (see the side-axis warp comment below), and rebuilding on every modest turn was itself the main cause of a distant peak's rendered shape visibly shifting while just panning across it. */
#define TERRAIN_BEHIND_M 250.0f      /* fixed, small, NOT ray-marched -- see ComputeVisibleDistanceKm's comment */
#define TERRAIN_MIN_EXTENT_M 200.0f  /* floor on forward/side extent so the window is never degenerately small */

/* Earth-curvature/refraction horizon-visibility model, shared with
 * ComputeVisibleDistanceKm below and RunProfileMode's sightline plot --
 * see the latter for the derivation (textbook ITU-R P.526-style bulge). */
#define VISIBILITY_REFRACTION_K (4.0 / 3.0)
#define VISIBILITY_MAX_ELEV_M 4500.0    /* headroom over the Bighorn tile's known ~4013m max -- see README.md */
#define VISIBILITY_MAX_RANGE_KM 20.0    /* hard ceiling: matches rlSetClipPlanes' far plane, and beyond this TERRAIN_RESOLUTION vertices is too coarse to bother */

/* Horizon-silhouette ray march along one compass bearing from
 * (camLat,camLon): how far out, in this direction, could anything loaded
 * actually be SEEN from an eye at eyeElevM (absolute elevation, meters),
 * given Earth's curvature, standard atmospheric refraction, and whatever
 * terrain relief actually exists along the way? This is exactly
 * "lazily depending on view direction, curvature, and intervening
 * obstacles" -- the caller uses the returned distance to decide how far
 * to extend the terrain mesh (and therefore how much DTM data actually
 * gets sampled/touched) in that direction, instead of always building a
 * fixed-size window regardless of what's actually visible.
 *
 * Uses the same curvature+refraction model as RunProfileMode's sightline
 * (destination_point, the d1*d2/(2*k*R) bulge term) but as a
 * single-observer running-max-angle horizon check rather than a
 * two-endpoint chord: a sample at distance d is "visible" only if its
 * curvature-corrected elevation angle beats every closer sample's --
 * otherwise a nearer ridge (or Earth's own curvature) hides it. */
static double ComputeVisibleDistanceKm(double camLat, double camLon, double eyeElevM, double bearingDeg) {
    double effectiveRadiusM = EARTH_RADIUS_KM * VISIBILITY_REFRACTION_K * 1000.0;
    double maxAngle = -1e18;
    double lastVisibleKm = 0.05; /* always at least a little ground right around the camera */

    for (double d = 0.05; d < VISIBILITY_MAX_RANGE_KM; d *= 1.4) {
        double dM = d * 1000.0;

        /* Bound check: even the tallest elevation anywhere in the loaded
         * DTMs can no longer beat maxAngle at this distance -- the
         * curvature penalty term only grows with distance, so once the
         * best possible case falls behind, nothing farther can ever catch
         * back up. Stop the ray here rather than sampling (and thus
         * touching/paging in) DTM data that provably couldn't be seen. */
        double boundAngle = (VISIBILITY_MAX_ELEV_M - eyeElevM) / dM - dM / (2.0 * effectiveRadiusM);
        if (boundAngle < maxAngle) break;

        double lat, lon;
        destination_point(camLat, camLon, bearingDeg, d, EARTH_RADIUS_KM, &lat, &lon);
        int level = d < 2.0 ? 0 : d < 8.0 ? 2 : 4; /* coarser DTM pyramid level for farther, cheaper samples -- same idea RunProfileMode uses */
        double m = DtmSampleMeters(&g_dtm, lat, lon, level);
        if (m == DTM_NODATA) break; /* no more loaded data exists this direction */

        double angle = (m - eyeElevM) / dM - dM / (2.0 * effectiveRadiusM);
        if (angle > maxAngle) {
            maxAngle = angle;
            lastVisibleKm = d;
        }
    }
    return lastVisibleKm;
}

/* Max of ComputeVisibleDistanceKm over several bearings, in meters --
 * used to size the terrain window's forward extent (a fan of bearings
 * across the camera's facing) and its side extent (bearings at +/-90). */
static float MaxVisibleExtentM(double camLat, double camLon, double eyeElevM, const double *bearingsDeg, int count) {
    double maxKm = 0.0;
    for (int i = 0; i < count; i++) {
        double km = ComputeVisibleDistanceKm(camLat, camLon, eyeElevM, bearingsDeg[i]);
        if (km > maxKm) maxKm = km;
    }
    float m = (float)(maxKm * 1000.0);
    return m < TERRAIN_MIN_EXTENT_M ? TERRAIN_MIN_EXTENT_M : m;
}

typedef struct {
    Mesh mesh;
    Material material;
    int built;
    double centerLat, centerLon; /* the (lat,lon) this window's local (0,0) corresponds to */
    float centerX, centerZ;      /* that same center, in the walk's own local-meters frame */
    float yaw;                   /* the facing this window's shape was built for -- see TERRAIN_REBUILD_YAW_DEG */
    float forwardM, sideM;       /* last computed visibility extents, for the HUD -- see MaxVisibleExtentM */
} TerrainWindow;

/* Whether the camera has strayed close enough to this window's own actual
 * boundary (in whichever direction) to need a rebuild -- NOT simply "moved
 * more than a fixed distance from center". The window's forward/side
 * extents are themselves visibility-driven (see BuildTerrainWindow) and
 * can be much smaller than any fixed threshold -- e.g. ~500m facing
 * straight into a nearby slope -- so a fixed-distance trigger can let the
 * camera walk right off the edge of a small window before ever
 * triggering, which is exactly what "walking forward and seeing under the
 * map" turned out to be: the window's edge, not a hole in it. */
static int TerrainWindowNeedsRebuild(const TerrainWindow *tw, Vector3 camPos, float camYaw) {
    float dx = camPos.x - tw->centerX, dz = camPos.z - tw->centerZ;
    float tyr = tw->yaw * DEG2RAD;
    float tFwdX = sinf(tyr), tFwdZ = -cosf(tyr);
    float tRightX = cosf(tyr), tRightZ = sinf(tyr);
    float alongForward = dx * tFwdX + dz * tFwdZ;
    float alongSide = dx * tRightX + dz * tRightZ;

    if (alongForward > tw->forwardM * TERRAIN_REBUILD_MARGIN) return 1;
    if (alongForward < -TERRAIN_BEHIND_M * TERRAIN_REBUILD_MARGIN) return 1;
    if (fabsf(alongSide) > tw->sideM * TERRAIN_REBUILD_MARGIN) return 1;
    if (fabsf(AngleDiffDeg(tw->yaw, camYaw)) > TERRAIN_REBUILD_YAW_DEG) return 1;
    return 0;
}

/* Elevation color ramp -- no orthophoto imagery is available, so terrain is
 * shaded by height: dark green low, tan/brown mid, white near the top of
 * whatever range this window actually spans (computed per-window, not a
 * fixed global range, so a low-relief window still shows visible shading). */
static Color ElevationColor(float t) { /* t in [0,1] */
    Color low  = (Color){ 46, 82, 42, 255 };
    Color mid  = (Color){ 133, 115, 74, 255 };
    Color high = (Color){ 245, 245, 250, 255 };
    if (t < 0.5f) {
        float u = t / 0.5f;
        return (Color){
            (unsigned char)(low.r + (mid.r - low.r) * u),
            (unsigned char)(low.g + (mid.g - low.g) * u),
            (unsigned char)(low.b + (mid.b - low.b) * u),
            255
        };
    }
    float u = (t - 0.5f) / 0.5f;
    return (Color){
        (unsigned char)(mid.r + (high.r - mid.r) * u),
        (unsigned char)(mid.g + (high.g - mid.g) * u),
        (unsigned char)(mid.b + (high.b - mid.b) * u),
        255
    };
}

/* Per-vertex working data for the two-pass mesh build below -- position
 * (already the mesh's own final local x/z, before the center translate
 * applied at draw time) plus height, so the second pass can compute real
 * normals from actual neighbor positions instead of assuming a uniform,
 * axis-aligned grid spacing (no longer true now that the window is
 * yaw-oriented and non-uniformly warped along its forward axis). */
typedef struct { float x, z, h; } TerrainVert;

static void BuildTerrainWindow(TerrainWindow *tw, double originLat, double originLon,
                                double centerLat, double centerLon, double originElevation, float yawDeg) {
    if (tw->built) UnloadMesh(tw->mesh);

    int res = TERRAIN_RESOLUTION;
    int vertsPerSide = res + 1;
    int vertexCount = vertsPerSide * vertsPerSide;
    int triangleCount = res * res * 2;

    Mesh mesh = { 0 };
    mesh.vertexCount = vertexCount;
    mesh.triangleCount = triangleCount;
    mesh.vertices = RL_MALLOC(sizeof(float) * 3 * (size_t)vertexCount);
    mesh.normals = RL_MALLOC(sizeof(float) * 3 * (size_t)vertexCount);
    mesh.colors = RL_MALLOC(sizeof(unsigned char) * 4 * (size_t)vertexCount);
    mesh.indices = RL_MALLOC(sizeof(unsigned short) * 3 * (size_t)triangleCount);
    if (!mesh.vertices || !mesh.normals || !mesh.colors || !mesh.indices) {
        fprintf(stderr, "BuildTerrainWindow: out of memory\n");
        exit(1);
    }

    double centerEastM, centerNorthM;
    DtmLonLatToLocalMeters(originLat, originLon, centerLat, centerLon, &centerEastM, &centerNorthM);

    /* Eye elevation at the window's center, for the visibility ray march --
     * falls back to the walk's starting elevation if this exact point
     * happens to be NODATA (shouldn't normally happen, since it's wherever
     * the camera currently is). */
    double centerElevM = DtmSampleMeters(&g_dtm, centerLat, centerLon, 0);
    double eyeElevM = (centerElevM == DTM_NODATA ? originElevation : centerElevM) + EYE_HEIGHT_M;

    /* Forward/right unit vectors in world (x,z), rotated by yawDeg -- same
     * convention as WalkForward/UpdateWalkCamera's strafe direction (right
     * = normalize(cross(flatForward, worldUp))), so "ahead" here always
     * matches whichever way the camera is actually facing. */
    float yr = yawDeg * DEG2RAD;
    float fwdX = sinf(yr), fwdZ = -cosf(yr);
    float rightX = cosf(yr), rightZ = sinf(yr);

    double forwardBearings[5] = { yawDeg - 40.0, yawDeg - 20.0, yawDeg, yawDeg + 20.0, yawDeg + 40.0 };
    double sideBearings[2] = { yawDeg - 90.0, yawDeg + 90.0 };
    float forwardM = MaxVisibleExtentM(centerLat, centerLon, eyeElevM, forwardBearings, 5);
    float sideM = MaxVisibleExtentM(centerLat, centerLon, eyeElevM, sideBearings, 2);
    float behindM = TERRAIN_BEHIND_M; /* fixed, not ray-marched -- see ComputeVisibleDistanceKm's comment header */

    /* Forward axis (index i) is split into a small uniform "behind" run and
     * a quadratically-warped "ahead" run -- dense near the camera (i just
     * past the behind section), sparse out toward forwardM -- so the fixed
     * TERRAIN_RESOLUTION vertex budget still gives good close-up detail
     * even when forwardM is many kilometers on a clear, unobstructed view.
     * Side axis (index j) is warped the same way, symmetrically about
     * sideOffset=0 (dense near the camera's own line of travel, sparse out
     * toward +/-sideM) -- it used to be uniformly spaced, which gave a
     * distant feature a noticeably different vertex density (and therefore
     * a different simplified silhouette) depending on whether it currently
     * fell under the forward axis's warp or the side axis's uniform
     * spacing -- e.g. a mountain visibly changing shape as you panned it
     * from "ahead" to "off to the side". Matching both axes' density
     * falloff makes a given distant point's approximation far more
     * consistent regardless of which direction it's in. */
    int behindVerts = vertsPerSide / 8;
    if (behindVerts < 2) behindVerts = 2;
    int aheadVerts = vertsPerSide - behindVerts;

    TerrainVert *verts = malloc(sizeof(TerrainVert) * (size_t)vertexCount);
    float minH = 1e9f, maxH = -1e9f;

    for (int j = 0; j < vertsPerSide; j++) {
        float sj = (float)j / res * 2.0f - 1.0f; /* -1..+1 */
        float sideOffset = (sj < 0 ? -1.0f : 1.0f) * sideM * sj * sj;
        for (int i = 0; i < vertsPerSide; i++) {
            float alongForward;
            if (i < behindVerts) {
                alongForward = -behindM + (float)i * (behindM / (behindVerts - 1));
            } else {
                float k = (float)(i - behindVerts) / (float)(aheadVerts - 1);
                alongForward = forwardM * k * k;
            }

            float x = fwdX * alongForward + rightX * sideOffset;
            float z = fwdZ * alongForward + rightZ * sideOffset;

            /* z is a world-Z offset, which per the file header convention
             * means south -- so it maps to DECREASING northM, the same
             * negation UpdateWalkCamera/the rebuild trigger apply via
             * `-wc->position.z`. See the bug this fixed: a from-scratch
             * reader might reasonably add rather than subtract here, which
             * silently samples each vertex's height from the north-south
             * MIRROR of where it actually sits in world space. */
            double localEast = centerEastM + x;
            double localNorth = centerNorthM - z;
            double lat, lon;
            DtmLocalMetersToLonLat(originLat, originLon, localEast, localNorth, &lat, &lon);

            /* Coarser DTM pyramid level for vertices far from the camera --
             * consistent with ComputeVisibleDistanceKm's own level choice,
             * and another concrete way the amount of DTM data actually
             * touched now tracks what's genuinely resolvable at that
             * distance rather than always reading full resolution. */
            float radialDist = sqrtf(x * x + z * z);
            int level = radialDist < 200.0f ? 0 : radialDist < 1000.0f ? 1 : radialDist < 4000.0f ? 2 : radialDist < 10000.0f ? 3 : 4;

            double m = DtmSampleMeters(&g_dtm, lat, lon, level);
            float h = (m == DTM_NODATA) ? 0.0f : (float)(m - originElevation);

            TerrainVert *tv = &verts[j * vertsPerSide + i];
            tv->x = x; tv->z = z; tv->h = h;
            if (h < minH) minH = h;
            if (h > maxH) maxH = h;
        }
    }
    if (maxH <= minH) maxH = minH + 1.0f; /* avoid divide-by-zero on a perfectly flat window */

    int v = 0;
    for (int j = 0; j < vertsPerSide; j++) {
        for (int i = 0; i < vertsPerSide; i++) {
            TerrainVert *tv = &verts[j * vertsPerSide + i];
            mesh.vertices[v * 3 + 0] = tv->x;
            mesh.vertices[v * 3 + 1] = tv->h;
            mesh.vertices[v * 3 + 2] = tv->z;

            /* Normal from actual neighbor positions (central differences,
             * one-sided at the window's edges) -- the grid is no longer
             * uniformly spaced or world-axis-aligned, so the old shortcut
             * (assume a constant "step" in both x and z) no longer applies;
             * this version works for any grid shape. Cross order verified
             * to give a +Y (up) normal for flat terrain, matching
             * rightDir x forwardDir = +Y (since forwardDir x rightDir =
             * forwardDir x (forwardDir x up) = -up by the vector triple
             * product identity, as forwardDir.up = 0). */
            int iL = i > 0 ? i - 1 : i, iR = i < res ? i + 1 : i;
            int jU = j > 0 ? j - 1 : j, jD = j < res ? j + 1 : j;
            TerrainVert *tL = &verts[j * vertsPerSide + iL], *tR = &verts[j * vertsPerSide + iR];
            TerrainVert *tU = &verts[jU * vertsPerSide + i], *tD = &verts[jD * vertsPerSide + i];
            Vector3 tangentI = { tR->x - tL->x, tR->h - tL->h, tR->z - tL->z };
            Vector3 tangentJ = { tD->x - tU->x, tD->h - tU->h, tD->z - tU->z };
            Vector3 n = Vector3Normalize(Vector3CrossProduct(tangentJ, tangentI));
            mesh.normals[v * 3 + 0] = n.x;
            mesh.normals[v * 3 + 1] = n.y;
            mesh.normals[v * 3 + 2] = n.z;

            Color c = ElevationColor((tv->h - minH) / (maxH - minH));
            mesh.colors[v * 4 + 0] = c.r;
            mesh.colors[v * 4 + 1] = c.g;
            mesh.colors[v * 4 + 2] = c.b;
            mesh.colors[v * 4 + 3] = 255;
            v++;
        }
    }
    free(verts);

    int idx = 0;
    for (int j = 0; j < res; j++) {
        for (int i = 0; i < res; i++) {
            unsigned short a = (unsigned short)(j * vertsPerSide + i);
            unsigned short b = (unsigned short)(a + vertsPerSide);
            unsigned short c = (unsigned short)(a + 1);
            unsigned short d = (unsigned short)(b + 1);
            mesh.indices[idx++] = a; mesh.indices[idx++] = b; mesh.indices[idx++] = c;
            mesh.indices[idx++] = c; mesh.indices[idx++] = b; mesh.indices[idx++] = d;
        }
    }

    UploadMesh(&mesh, false);
    tw->mesh = mesh;
    if (!tw->built) tw->material = LoadMaterialDefault();
    tw->built = 1;
    tw->centerLat = centerLat;
    tw->centerLon = centerLon;
    tw->centerX = (float)centerEastM;
    tw->centerZ = (float)centerNorthM;
    tw->yaw = yawDeg;
    tw->forwardM = forwardM;
    tw->sideM = sideM;
}

/* --- Walk mode --- */

typedef struct {
    Vector3 position; /* world meters, relative to walk origin */
    float yaw, pitch;  /* degrees */
    float fovy;        /* degrees; scroll wheel zooms, see FOV_* constants */
} WalkCamera;

static Vector3 WalkForward(const WalkCamera *wc) {
    float yr = wc->yaw * DEG2RAD, pr = wc->pitch * DEG2RAD;
    return (Vector3){ sinf(yr) * cosf(pr), sinf(pr), -cosf(yr) * cosf(pr) };
}

static void UpdateWalkCamera(WalkCamera *wc, double originLat, double originLon, double originElevation, float deltaTime) {
    /* Headless screenshot runs (see TV_SCREENSHOT in RunWalkMode) have no
     * real pointer -- GLFW/X11 can report a spurious nonzero delta on such
     * a virtual/unfocused window, which would make TV_YAW/TV_PITCH
     * non-deterministic frame to frame. Skip mouse-look entirely in that
     * mode; interactive runs are unaffected. */
    if (!getenv("TV_SCREENSHOT")) {
        Vector2 md = GetMouseDelta();
        wc->yaw += md.x * MOUSE_SENSITIVITY;
        wc->pitch -= md.y * MOUSE_SENSITIVITY;

        wc->fovy -= GetMouseWheelMove() * FOV_ZOOM_SPEED_DEG;
        if (wc->fovy < FOV_MIN_DEG) wc->fovy = FOV_MIN_DEG;
        if (wc->fovy > FOV_MAX_DEG) wc->fovy = FOV_MAX_DEG;
    }
    if (wc->pitch > 89.0f) wc->pitch = 89.0f;
    if (wc->pitch < -89.0f) wc->pitch = -89.0f;

    Vector3 forward = WalkForward(wc);
    Vector3 flatForward = Vector3Normalize((Vector3){ forward.x, 0, forward.z });
    Vector3 right = Vector3Normalize(Vector3CrossProduct(flatForward, (Vector3){ 0, 1, 0 }));

    float speed = WALK_SPEED_MPS * (IsKeyDown(KEY_LEFT_SHIFT) ? SPRINT_MULT : 1.0f);
    Vector3 move = { 0, 0, 0 };
    if (IsKeyDown(KEY_W) || IsKeyDown(KEY_UP)) move = Vector3Add(move, flatForward);
    if (IsKeyDown(KEY_S) || IsKeyDown(KEY_DOWN)) move = Vector3Subtract(move, flatForward);
    if (IsKeyDown(KEY_D)) move = Vector3Add(move, right);
    if (IsKeyDown(KEY_A)) move = Vector3Subtract(move, right);
    if (Vector3Length(move) > 0.0001f) {
        move = Vector3Scale(Vector3Normalize(move), speed * deltaTime);
        wc->position.x += move.x;
        wc->position.z += move.z;
    }

    double lat, lon;
    DtmLocalMetersToLonLat(originLat, originLon, wc->position.x, -wc->position.z, &lat, &lon);
    double groundM = DtmSampleMeters(&g_dtm, lat, lon, 0);
    float groundY = (groundM == DTM_NODATA) ? wc->position.y - EYE_HEIGHT_M : (float)(groundM - originElevation);
    wc->position.y = groundY + EYE_HEIGHT_M;
}

static int RunWalkMode(double startLat, double startLon) {
    double originLat = startLat, originLon = startLon;
    EnsureCoverageOrWarn(originLat, originLon);
    double originElevation = DtmSampleMeters(&g_dtm, originLat, originLon, 0);
    if (originElevation == DTM_NODATA) {
        fprintf(stderr, "terrain_viewer: (%.5f, %.5f) isn't covered by any loaded DTM file\n", startLat, startLon);
        return 1;
    }
    printf("Starting elevation: %.1f m\n", originElevation);

    SetConfigFlags(FLAG_WINDOW_RESIZABLE | FLAG_VSYNC_HINT);
    InitWindow(1280, 720, "DTM Terrain Walkthrough");
    if (!getenv("TV_SCREENSHOT")) DisableCursor(); /* headless screenshot runs don't want an input-captured window */

    WalkCamera wc = { 0 };
    wc.position = (Vector3){ 0, (float)EYE_HEIGHT_M, 0 };
    wc.yaw = getenv("TV_YAW") ? (float)atof(getenv("TV_YAW")) : 0.0f;
    wc.pitch = getenv("TV_PITCH") ? (float)atof(getenv("TV_PITCH")) : -5.0f;
    wc.fovy = getenv("TV_FOV") ? (float)atof(getenv("TV_FOV")) : FOV_DEFAULT_DEG;

    TerrainWindow terrain = { 0 };
    BuildTerrainWindow(&terrain, originLat, originLon, originLat, originLon, originElevation, wc.yaw);

    Camera3D camera = { 0 };
    camera.projection = CAMERA_PERSPECTIVE;
    camera.up = (Vector3){ 0, 1, 0 };

    const char *screenshotPath = getenv("TV_SCREENSHOT");
    int screenshotFrame = getenv("TV_SCREENSHOT_FRAME") ? atoi(getenv("TV_SCREENSHOT_FRAME")) : 30;
    int frameCount = 0;

    SetTargetFPS(60);
    rlSetClipPlanes(0.1, 20000.0);

    while (!WindowShouldClose()) {
        frameCount++;
        float dt = GetFrameTime();
        UpdateWalkCamera(&wc, originLat, originLon, originElevation, dt);

        if (TerrainWindowNeedsRebuild(&terrain, wc.position, wc.yaw)) {
            double lat, lon;
            DtmLocalMetersToLonLat(originLat, originLon, wc.position.x, -wc.position.z, &lat, &lon);
            BuildTerrainWindow(&terrain, originLat, originLon, lat, lon, originElevation, wc.yaw);
        }

        camera.position = wc.position;
        Vector3 forward = WalkForward(&wc);
        camera.target = Vector3Add(wc.position, forward);
        camera.fovy = wc.fovy;

        BeginDrawing();
        ClearBackground((Color){ 137, 187, 227, 255 });

        BeginMode3D(camera);
            rlEnableBackfaceCulling();
            /* Terrain vertices are already in the walk's local-meters frame
             * centered on terrain.centerX/centerZ -- translate the mesh
             * (not the vertices) so it renders at its real offset from
             * wherever the camera currently is. */
            Matrix xform = MatrixTranslate(terrain.centerX, 0, terrain.centerZ);
            DrawMesh(terrain.mesh, terrain.material, xform);
        EndMode3D();

        DrawFPS(10, 10);
        {
            double lat, lon;
            DtmLocalMetersToLonLat(originLat, originLon, wc.position.x, -wc.position.z, &lat, &lon);
            double elev = wc.position.y - EYE_HEIGHT_M + originElevation;
            DrawText(TextFormat("lat=%.5f  lon=%.5f  elev=%.1fm  yaw=%.0f pitch=%.0f", lat, lon, elev, wc.yaw, wc.pitch), 10, 35, 18, RAYWHITE);
            DrawText(TextFormat("visible range: %.0fm ahead, %.0fm to the sides (curvature/obstacle-limited)", terrain.forwardM, terrain.sideM), 10, 58, 16, RAYWHITE);
        }
        {
            float heading = CompassHeadingDeg(wc.yaw);
            const char *compassText = TextFormat("%03.0f deg %s", heading, CompassLabel(heading));
            int fontSize = 26;
            int tw = MeasureText(compassText, fontSize);
            DrawText(compassText, GetScreenWidth() - tw - 12, 10, fontSize, RAYWHITE);

            const char *zoomText = TextFormat("FOV %.0f deg (scroll to zoom)", wc.fovy);
            int tw2 = MeasureText(zoomText, 16);
            DrawText(zoomText, GetScreenWidth() - tw2 - 12, 40, 16, SKYBLUE);
        }
        DrawText("WASD move, Shift sprint, mouse look, scroll to zoom, Esc to quit", 10, GetScreenHeight() - 30, 14, SKYBLUE);
        EndDrawing();

        if (screenshotPath && frameCount == screenshotFrame) {
            TakeScreenshot(screenshotPath);
            break;
        }
    }

    if (terrain.built) { UnloadMesh(terrain.mesh); UnloadMaterial(terrain.material); }
    CloseWindow();
    return 0;
}

/* --- Profile mode: elevation vs. distance between two points, with Earth's
 * curvature shown relative to the straight sightline chord between the two
 * endpoints (the standard "does the terrain block line-of-sight" plot). --- */

#define PROFILE_SAMPLES 500

#define FEET_TO_METERS 0.3048

/* Parses a +<number>[m|f] height argument (e.g. "+195", "+50m", "+30f") into
 * feet above ground -- 'm'/'M' means the number given is meters, 'f'/'F' or
 * no suffix means feet, matching this project's original convention (kept
 * as the default since it predates the unit suffix). Internally everything
 * still runs in feet, same as before, so this is purely an input-parsing
 * change -- display/math elsewhere in the file is untouched. */
static double ParseHeightFeetArg(const char *arg) {
    char *end;
    double value = strtod(arg + 1, &end); /* +1 skips the leading '+' */
    if (*end == 'm' || *end == 'M') return value / FEET_TO_METERS;
    return value;
}

static int RunProfileMode(double lat1, double lon1, double heightFt1, double lat2, double lon2, double heightFt2, const char *screenshotPath) {
    double totalKm = haversine_distance(lat1, lon1, lat2, lon2, EARTH_RADIUS_KM);
    double bearing = initial_bearing(lat1, lon1, lat2, lon2);
    printf("Profile: (%.5f,%.5f)+%.0fft -> (%.5f,%.5f)+%.0fft, %.3f km, initial bearing %.1f deg\n",
           lat1, lon1, heightFt1, lat2, lon2, heightFt2, totalKm, bearing);

    EnsureCoverageOrWarn(lat1, lon1);
    EnsureCoverageOrWarn(lat2, lon2);
    for (double frac = 0.15; frac < 1.0; frac += 0.15) {
        double lat, lon;
        destination_point(lat1, lon1, bearing, totalKm * frac, EARTH_RADIUS_KM, &lat, &lon);
        EnsureCoverageOrWarn(lat, lon);
    }

    double distKm[PROFILE_SAMPLES], elevM[PROFILE_SAMPLES];
    int haveData[PROFILE_SAMPLES];
    double elev1 = DtmSampleMeters(&g_dtm, lat1, lon1, 2);
    double elev2 = DtmSampleMeters(&g_dtm, lat2, lon2, 2);

    for (int i = 0; i < PROFILE_SAMPLES; i++) {
        double d = totalKm * i / (PROFILE_SAMPLES - 1);
        double lat, lon;
        if (i == PROFILE_SAMPLES - 1) { lat = lat2; lon = lon2; }
        else destination_point(lat1, lon1, bearing, d, EARTH_RADIUS_KM, &lat, &lon);
        double m = DtmSampleMeters(&g_dtm, lat, lon, 2); /* level 2 (~4m/px) plenty for a multi-km profile */
        distKm[i] = d;
        haveData[i] = (m != DTM_NODATA);
        elevM[i] = haveData[i] ? m : 0.0;
    }

    /* Curvature + standard 4/3-Earth-radius refraction correction. Textbook
     * (ITU-R P.526-style) formula: a straight line of sight between two
     * points at distances d1/d2 from each end (d1+d2 = totalKm) clears the
     * curved Earth's surface by less than flat-plane geometry suggests, by
     * bulge(d) = d1*d2 / (2*k*R) -- e.g. two observers at sea level far
     * enough apart can't see each other at all because the ocean itself
     * "bulges up" between them by more than their eye height. So a terrain
     * point at distance d effectively needs to clear chord(d) + bulge(d),
     * not just the flat chord(d), to be visible -- sightlineM below is that
     * raised reference curve, and apparentM is terrain height above it
     * (positive = blocks line-of-sight). */
    double effectiveRadiusKm = EARTH_RADIUS_KM * (4.0 / 3.0);
    /* groundValid* is the bare terrain elevation at each endpoint (for
     * drawing the mast up from the ground below); elevValid* adds the
     * +heightFt mast/tower height on top of it -- e.g. a cell tower's
     * antenna height above ground, not a person standing there -- and is
     * what actually anchors the sightline chord/bulge and the endpoint
     * markers below. */
    double groundValid1 = (elev1 == DTM_NODATA) ? elevM[0] : elev1;
    double groundValid2 = (elev2 == DTM_NODATA) ? elevM[PROFILE_SAMPLES - 1] : elev2;
    double elevValid1 = groundValid1 + heightFt1 * FEET_TO_METERS;
    double elevValid2 = groundValid2 + heightFt2 * FEET_TO_METERS;
    double sightlineM[PROFILE_SAMPLES], apparentM[PROFILE_SAMPLES];
    for (int i = 0; i < PROFILE_SAMPLES; i++) {
        double d = distKm[i];
        double chordM = elevValid1 + (elevValid2 - elevValid1) * (d / totalKm);
        double bulgeM = 1000.0 * (d * (totalKm - d)) / (2.0 * effectiveRadiusKm);
        sightlineM[i] = chordM + bulgeM;
        apparentM[i] = elevM[i] - sightlineM[i];
    }

    SetConfigFlags(FLAG_WINDOW_RESIZABLE | FLAG_VSYNC_HINT);
    InitWindow(1280, 640, "DTM Elevation Profile");
    SetTargetFPS(30);

    int frameCount = 0;
    int screenshotFrame = getenv("TV_SCREENSHOT_FRAME") ? atoi(getenv("TV_SCREENSHOT_FRAME")) : 5;

    while (!WindowShouldClose()) {
        frameCount++;
        int w = GetScreenWidth(), h = GetScreenHeight();
        int marginL = 70, marginR = 30, marginT = 50, marginB = 60;
        int plotW = w - marginL - marginR, plotH = h - marginT - marginB;

        double minE = 1e18, maxE = -1e18, minS = 1e18, maxS = -1e18;
        for (int i = 0; i < PROFILE_SAMPLES; i++) {
            if (sightlineM[i] < minS) minS = sightlineM[i];
            if (sightlineM[i] > maxS) maxS = sightlineM[i];
            if (!haveData[i]) continue;
            if (elevM[i] < minE) minE = elevM[i];
            if (elevM[i] > maxE) maxE = elevM[i];
        }
        double loE = minE < minS ? minE : minS, hiE = maxE > maxS ? maxE : maxS;
        if (hiE <= loE) hiE = loE + 1.0;
        double pad = (hiE - loE) * 0.08;
        loE -= pad; hiE += pad;

        #define SX(d) (marginL + (int)((d) / totalKm * plotW))
        #define SY(e) (marginT + plotH - (int)(((e) - loE) / (hiE - loE) * plotH))

        BeginDrawing();
        ClearBackground(RAYWHITE);

        DrawRectangleLines(marginL, marginT, plotW, plotH, GRAY);
        for (int g = 0; g <= 4; g++) {
            double e = loE + (hiE - loE) * g / 4.0;
            int y = SY(e);
            DrawLine(marginL, y, marginL + plotW, y, (Color){ 220, 220, 220, 255 });
            DrawText(TextFormat("%.0fm", e), 5, y - 6, 14, DARKGRAY);
        }
        for (int g = 0; g <= 5; g++) {
            double d = totalKm * g / 5.0;
            int x = SX(d);
            DrawText(TextFormat("%.1fkm", d), x - 15, marginT + plotH + 6, 14, DARKGRAY);
        }

        /* Raw terrain (brown), the curved sightline reference (red,
         * chord+bulge -- see the comment above where it's computed), and a
         * thicker red overlay on any stretch where terrain pokes above that
         * sightline, i.e. actually blocks line-of-sight. */
        for (int i = 0; i < PROFILE_SAMPLES - 1; i++) {
            if (!haveData[i] || !haveData[i + 1]) continue;
            DrawLine(SX(distKm[i]), SY(elevM[i]), SX(distKm[i + 1]), SY(elevM[i + 1]), (Color){ 90, 60, 20, 255 });
        }
        for (int i = 0; i < PROFILE_SAMPLES - 1; i++) {
            DrawLine(SX(distKm[i]), SY(sightlineM[i]), SX(distKm[i + 1]), SY(sightlineM[i + 1]), (Color){ 200, 30, 30, 160 });
        }
        for (int i = 0; i < PROFILE_SAMPLES - 1; i++) {
            if (!haveData[i] || !haveData[i + 1]) continue;
            if (apparentM[i] > 0 || apparentM[i + 1] > 0) {
                DrawLineEx((Vector2){ (float)SX(distKm[i]), (float)SY(elevM[i]) },
                           (Vector2){ (float)SX(distKm[i + 1]), (float)SY(elevM[i + 1]) }, 3.0f, (Color){ 220, 40, 40, 255 });
            }
        }

        /* Mast/tower line from the ground up to the actual sightline
         * height, when a +heightFt was given -- makes clear the marker is
         * an antenna height above ground, not a person standing there. */
        if (heightFt1 > 0.0) DrawLine(SX(0), SY(groundValid1), SX(0), SY(elevValid1), (Color){ 90, 90, 90, 255 });
        if (heightFt2 > 0.0) DrawLine(SX(totalKm), SY(groundValid2), SX(totalKm), SY(elevValid2), (Color){ 90, 90, 90, 255 });
        DrawCircle(SX(0), SY(elevValid1), 5, BLUE);
        DrawCircle(SX(totalKm), SY(elevValid2), 5, BLUE);
        DrawText(TextFormat("Total distance: %.3f km   (red = terrain blocking line-of-sight; thin red line = direct sightline incl. curvature+refraction)", totalKm), marginL, 8, 15, DARKGRAY);
        if (heightFt1 > 0.0) {
            DrawText(TextFormat("A (%.4f,%.4f) %.0fm ground + %.0fft = %.0fm", lat1, lon1, groundValid1, heightFt1, elevValid1), marginL, marginT - 22, 16, BLUE);
        } else {
            DrawText(TextFormat("A (%.4f,%.4f) %.0fm", lat1, lon1, elevValid1), marginL, marginT - 22, 16, BLUE);
        }
        {
            const char *bText = heightFt2 > 0.0
                ? TextFormat("B (%.4f,%.4f) %.0fm ground + %.0fft = %.0fm", lat2, lon2, groundValid2, heightFt2, elevValid2)
                : TextFormat("B (%.4f,%.4f) %.0fm", lat2, lon2, elevValid2);
            DrawText(bText, marginL + plotW - MeasureText(bText, 16), marginT - 22, 16, BLUE);
        }

        EndDrawing();

        if (screenshotPath && frameCount == screenshotFrame) {
            TakeScreenshot(screenshotPath);
            break;
        }
    }

    CloseWindow();
    return 0;
}

/* --- Viewshed mode: given a tower (lat, lon, antenna height), rasterize a
 * top-down map of which ground points within maxDistanceKm have a clear,
 * curvature-and-terrain-corrected line of sight to it -- e.g. for scouting
 * where cellular coverage from a given tower actually reaches versus where
 * a ridge blocks it. This is a pure geometric line-of-sight map: no path
 * loss, antenna radiation pattern, Fresnel-zone clearance, or signal
 * strength is modeled -- green means "a straight line to the tower clears
 * the terrain", not "you'll have usable signal there". See README.md. */

#define VIEWSHED_MAP_PX 900              /* square raster size, in pixels */
#define VIEWSHED_RECEIVER_HEIGHT_M 1.5   /* typical handset/vehicle height above ground, added at every point being TESTED for visibility -- but NOT to the terrain profile itself, which is what actually blocks the line (see ComputeViewshedRaster's comment) */

/* Optional RF signal-strength model, layered on top of the pure-geometry
 * viewshed when the caller opts in (rf->enabled) by passing --freq/--erp
 * on the command line. Two pieces, both standard textbook formulas:
 *
 * - Free-space path loss: FSPL(dB) = 20*log10(d_km) + 20*log10(f_MHz) +
 *   32.44 -- how much a signal weakens with distance alone, clear sky.
 *
 * - Single knife-edge diffraction loss (ITU-R P.526, same family as the
 *   curvature/refraction formulas already used elsewhere in this file):
 *   given the ONE most-obstructing terrain point found so far along a ray
 *   (tracked during the same running-horizon march the pure-geometry
 *   viewshed already does -- see ComputeViewshedRaster), compute how far
 *   that obstruction poked above the direct tower-to-target line (h,
 *   accounting for Earth's curvature bulge at the obstruction's position)
 *   relative to the signal's own Fresnel-zone size at that geometry (the
 *   diffraction parameter nu), and from that a smooth, continuous loss in
 *   dB -- ~0 for a well-clear path, rising through "marginal Fresnel
 *   clearance" even before full geometric blockage, further rising once
 *   truly blocked. This is why a spot just barely hidden behind a ridge
 *   shows as weak-but-present signal rather than a hard cutoff, which is
 *   physically realistic and matches how real cell edges behave.
 *
 * What this deliberately does NOT model: multiple obstructions along one
 * path (this tracks only the single dominant one, not full multi-edge
 * diffraction like Longley-Rice/ITM), antenna radiation pattern
 * (isotropic assumed), ground conductivity/reflection, foliage/building
 * clutter, or troposcatter. It's a solid, physically-grounded estimate --
 * not the same fidelity as a professional RF planning tool. See README.md.
 */
typedef struct {
    int enabled;
    double freqMHz;
    double eirpDbm;
    double sensDbm; /* receiver sensitivity threshold -- below this, marked as no usable coverage regardless of the continuous color ramp */

    /* --itm: swap the single-knife-edge model above for a real port of the
     * NTIA/ITS Irregular Terrain Model (Longley-Rice), point-to-point mode
     * -- see longley_rice.c/.h. Needs a handful of regional constants the
     * simple model doesn't: ground electrical properties and atmosphere/
     * climate, all defaulted to standard "average" values (see main()'s
     * argument parsing) rather than requiring the user to look anything
     * up, matching how real ITM-based studies work unless someone has a
     * specific reason to override them. */
    int useItm;
    int climate;   /* CLIMATE__* from longley_rice.h */
    int pol;       /* POLARIZATION__* from longley_rice.h */
    double epsilon; /* relative permittivity (dielectric constant) of the ground */
    double sigma;   /* ground conductivity, S/m */
    double N0;      /* surface refractivity, N-units */
} RfParams;

#define SPEED_OF_LIGHT_M_S 299792458.0

static double FreeSpacePathLossDb(double distKm, double freqMHz) {
    if (distKm < 0.001) distKm = 0.001; /* avoid log(0) right at the tower */
    return 20.0 * log10(distKm) + 20.0 * log10(freqMHz) + 32.44;
}

static double KnifeEdgeDiffractionLossDb(double hM, double d1Km, double d2Km, double freqMHz) {
    if (d1Km <= 0.0 || d2Km <= 0.0) return 0.0; /* no controlling obstruction registered yet */
    double wavelengthM = SPEED_OF_LIGHT_M_S / (freqMHz * 1.0e6);
    double d1M = d1Km * 1000.0, d2M = d2Km * 1000.0;
    double nu = hM * sqrt(2.0 * (d1M + d2M) / (wavelengthM * d1M * d2M));
    if (nu <= -0.7) return 0.0; /* well clear of the obstruction */
    double L = 6.9 + 20.0 * log10(sqrt((nu - 0.1) * (nu - 0.1) + 1.0) + nu - 0.1);
    return L < 0.0 ? 0.0 : L;
}

/* Color for a point marginDb above (or below) the sensitivity threshold --
 * gray if it doesn't reach threshold at all (no usable coverage), else a
 * red(weak, right at threshold)-yellow-green(strong, 40dB+ of margin)
 * ramp, the same shape a real carrier coverage map uses. */
static Color SignalColor(double marginDb, double ceilingDb) {
    if (marginDb < 0.0) return (Color){ 140, 140, 140, 255 };
    double t = marginDb / ceilingDb;
    if (t > 1.0) t = 1.0;
    Color weak = (Color){ 205, 60, 40, 255 }, mid = (Color){ 230, 200, 40, 255 }, strong = (Color){ 40, 140, 60, 255 };
    Color a, b; double u;
    if (t < 0.5) { a = weak; b = mid; u = t / 0.5; } else { a = mid; b = strong; u = (t - 0.5) / 0.5; }
    return (Color){
        (unsigned char)(a.r + (b.r - a.r) * u), (unsigned char)(a.g + (b.g - a.g) * u),
        (unsigned char)(a.b + (b.b - a.b) * u), 255
    };
}

/* Same running-max-angle horizon algorithm as ComputeVisibleDistanceKm,
 * generalized to mark every sample along every bearing (not just the
 * farthest visible one) and rasterize the result directly into `pixels`
 * (VIEWSHED_MAP_PX x VIEWSHED_MAP_PX, tower at the center, north up).
 *
 * Two different elevations matter at each step and must NOT be conflated:
 * groundAngle (bare terrain) is what updates the running horizon, because
 * real terrain blocks everything farther along the ray regardless of who's
 * asking; receiverAngle (terrain + VIEWSHED_RECEIVER_HEIGHT_M) is what
 * gets compared against that horizon to decide THIS point's own
 * visibility, because a receiver a couple meters off the ground can see
 * past an obstruction a person standing at ground level couldn't -- but
 * that same couple meters shouldn't let a receiver claim it can see over
 * its OWN terrain and thereby suppress everything behind it. Standard
 * "observer offset / target offset" GIS viewshed practice, just applied
 * along one ray at a time. */
/* Returns the ceiling (dB above sensitivity) the RF color ramp actually
 * used -- see the big comment below for why this can't just be a fixed
 * 40dB. Meaningless (0.0) when rf->enabled is false. */
static double ComputeViewshedRaster(Color *pixels, double towerLat, double towerLon, double towerElevM, double maxDistanceKm, const RfParams *rf) {
    double effectiveRadiusKm = EARTH_RADIUS_KM * VISIBILITY_REFRACTION_K;
    double effectiveRadiusM = effectiveRadiusKm * 1000.0;
    double metersPerPixel = (maxDistanceKm * 2000.0) / VIEWSHED_MAP_PX;
    double stepKm = metersPerPixel / 1000.0;
    if (stepKm < 0.005) stepKm = 0.005;

    double angleStepDeg = (metersPerPixel / (maxDistanceKm * 1000.0)) * (180.0 / M_PI);
    if (angleStepDeg < 0.05) angleStepDeg = 0.05;
    if (angleStepDeg > 2.0) angleStepDeg = 2.0;
    int bearingCount = (int)(360.0 / angleStepDeg);
    if (bearingCount < 360) bearingCount = 360;
    if (bearingCount > 3600) bearingCount = 3600;

    Color visibleColor = (Color){ 70, 160, 90, 255 };
    Color blockedColor = (Color){ 195, 60, 50, 255 };

    /* A fixed 0-40dB ramp (the old behavior) saturates solid green the
     * moment a scenario has more than 40dB of link margin anywhere --
     * which is most of the map for a high-power transmitter at short
     * range (e.g. 140W at 1.5km leaves ~59dB of margin over a clear
     * path, versus the 40dB ceiling). So: run the whole sweep storing
     * each point's raw margin into `margin[]` instead of a color, track
     * the strongest margin actually seen, then color every point in a
     * second pass against THAT ceiling -- the full red-to-green range is
     * always in play, at whatever scale this particular tower/power/
     * distance combination actually produces. NAN marks a pixel no
     * bearing ever reached (stays whatever `pixels[]` was pre-filled
     * with, e.g. the "not computed" background). */
    float *margin = NULL;
    double maxMarginSeen = -1e18;
    if (rf->enabled) {
        margin = malloc(sizeof(float) * VIEWSHED_MAP_PX * VIEWSHED_MAP_PX);
        for (int i = 0; i < VIEWSHED_MAP_PX * VIEWSHED_MAP_PX; i++) margin[i] = NAN;
    }

    if (rf->enabled && rf->useItm) {
        /* ITM needs the actual digitized terrain profile from tower to
         * target (pfl[]), not just a running horizon angle, and each
         * evaluation costs O(profile length) instead of the knife-edge
         * model's O(1) -- calling it at the same fine stepKm resolution
         * the geometry sweep uses (often thousands of steps per bearing)
         * would be far too slow. So this uses a much coarser radial
         * spacing of its own (itmStepKm, capped at ~150 points per ray),
         * building pfl[] incrementally as it walks outward -- exactly the
         * evenly-spaced/elevations-only format the ported reference
         * implementation expects -- and splats each result over a bigger
         * box than the geometry loop's 3x3 to cover the wider radial
         * gaps this coarser spacing leaves. Full angular resolution
         * (bearingCount, same as the geometry sweep) is kept so there are
         * no angular gaps/streaks. */
        double towerGroundM = DtmSampleMeters(&g_dtm, towerLat, towerLon, 0);
        double hTxM = towerElevM - towerGroundM;
        double itmStepKm = maxDistanceKm / 150.0;
        if (itmStepKm < stepKm) itmStepKm = stepKm;
        int maxPts = (int)(maxDistanceKm / itmStepKm) + 4;
        double *pfl = malloc(sizeof(double) * (size_t)(maxPts + 2));
        int splatR = (int)ceil((itmStepKm * 1000.0 / metersPerPixel) / 2.0);
        if (splatR < 1) splatR = 1;
        if (splatR > 6) splatR = 6;

        for (int b = 0; b < bearingCount; b++) {
            double bearing = 360.0 * b / bearingCount;
            double bearingRad = bearing * DEG2RAD;
            double sinB = sin(bearingRad), cosB = cos(bearingRad);
            pfl[2] = towerGroundM; /* profile point 0: the tower's own bare-ground elevation */
            int np = 0;

            for (double d = itmStepKm; d <= maxDistanceKm; d += itmStepKm) {
                double lat, lon;
                destination_point(towerLat, towerLon, bearing, d, EARTH_RADIUS_KM, &lat, &lon);
                int level = d < 2.0 ? 0 : d < 8.0 ? 2 : 4;
                double ground = DtmSampleMeters(&g_dtm, lat, lon, level);
                if (ground == DTM_NODATA) break; /* no more loaded data this direction */
                np++;
                pfl[2 + np] = ground;
                pfl[0] = (double)np;
                pfl[1] = itmStepKm * 1000.0;

                double A_db;
                long warnings;
                int rtn = ITM_P2P_TLS(hTxM, VIEWSHED_RECEIVER_HEIGHT_M, pfl, rf->climate, rf->N0, rf->freqMHz,
                                       rf->pol, rf->epsilon, rf->sigma, MDVAR__BROADCAST_MODE, 50.0, 50.0, 50.0,
                                       &A_db, &warnings);
                if (rtn >= 1000) continue; /* invalid input at this profile length (e.g. path still too short) -- skip, try farther out */

                double thisMargin = rf->eirpDbm - A_db - rf->sensDbm;
                if (thisMargin > maxMarginSeen) maxMarginSeen = thisMargin;

                double dM = d * 1000.0;
                double eastM = dM * sinB, northM = dM * cosB;
                int px = VIEWSHED_MAP_PX / 2 + (int)lround(eastM / metersPerPixel);
                int py = VIEWSHED_MAP_PX / 2 - (int)lround(northM / metersPerPixel);
                for (int oy = -splatR; oy <= splatR; oy++) {
                    int yy = py + oy;
                    if (yy < 0 || yy >= VIEWSHED_MAP_PX) continue;
                    for (int ox = -splatR; ox <= splatR; ox++) {
                        int xx = px + ox;
                        if (xx < 0 || xx >= VIEWSHED_MAP_PX) continue;
                        margin[yy * VIEWSHED_MAP_PX + xx] = (float)thisMargin;
                    }
                }
            }
        }
        free(pfl);
    } else {
        for (int b = 0; b < bearingCount; b++) {
            double bearing = 360.0 * b / bearingCount;
            double bearingRad = bearing * DEG2RAD;
            double sinB = sin(bearingRad), cosB = cos(bearingRad);
            double maxAngle = -1e18;
            double obstDistKm = -1.0, obstElevM = 0.0; /* the single controlling (worst-so-far) obstruction along this ray, for the RF model's diffraction term -- see the RfParams comment */

            for (double d = stepKm; d <= maxDistanceKm; d += stepKm) {
                double lat, lon;
                destination_point(towerLat, towerLon, bearing, d, EARTH_RADIUS_KM, &lat, &lon);
                double dM = d * 1000.0;
                int level = d < 2.0 ? 0 : d < 8.0 ? 2 : 4;
                double ground = DtmSampleMeters(&g_dtm, lat, lon, level);
                if (ground == DTM_NODATA) break; /* no more loaded data this direction */

                double groundAngle = (ground - towerElevM) / dM - dM / (2.0 * effectiveRadiusM);
                double rxElevM = ground + VIEWSHED_RECEIVER_HEIGHT_M;
                double receiverAngle = (rxElevM - towerElevM) / dM - dM / (2.0 * effectiveRadiusM);
                int visible = receiverAngle > maxAngle;

                Color c = blockedColor;
                double thisMargin = 0.0;
                if (rf->enabled) {
                    double diffLossDb = 0.0;
                    if (obstDistKm > 0.0) {
                        double d1 = obstDistKm, d2 = d - obstDistKm;
                        double bulgeAtObstM = 1000.0 * (d1 * d2) / (2.0 * effectiveRadiusKm);
                        double lineElevAtObstM = towerElevM + (rxElevM - towerElevM) * (d1 / d) + bulgeAtObstM;
                        diffLossDb = KnifeEdgeDiffractionLossDb(obstElevM - lineElevAtObstM, d1, d2, rf->freqMHz);
                    }
                    double prDbm = rf->eirpDbm - FreeSpacePathLossDb(d, rf->freqMHz) - diffLossDb;
                    thisMargin = prDbm - rf->sensDbm;
                    if (thisMargin > maxMarginSeen) maxMarginSeen = thisMargin;
                } else {
                    c = visible ? visibleColor : blockedColor;
                }

                if (groundAngle > maxAngle) { maxAngle = groundAngle; obstDistKm = d; obstElevM = ground; }

                double eastM = dM * sinB, northM = dM * cosB; /* bearing convention matches WalkForward/destination_point: 0=north, clockwise */
                int px = VIEWSHED_MAP_PX / 2 + (int)lround(eastM / metersPerPixel);
                int py = VIEWSHED_MAP_PX / 2 - (int)lround(northM / metersPerPixel); /* image Y grows downward; north is up */

                for (int oy = -1; oy <= 1; oy++) {
                    int yy = py + oy;
                    if (yy < 0 || yy >= VIEWSHED_MAP_PX) continue;
                    for (int ox = -1; ox <= 1; ox++) {
                        int xx = px + ox;
                        if (xx < 0 || xx >= VIEWSHED_MAP_PX) continue;
                        if (rf->enabled) margin[yy * VIEWSHED_MAP_PX + xx] = (float)thisMargin;
                        else pixels[yy * VIEWSHED_MAP_PX + xx] = c;
                    }
                }
            }
        }
    }

    if (!rf->enabled) return 0.0;

    double ceilingDb = maxMarginSeen > 3.0 ? maxMarginSeen : 40.0; /* degenerate "nothing has signal" case -- scale doesn't matter, the map's all gray anyway */
    for (int i = 0; i < VIEWSHED_MAP_PX * VIEWSHED_MAP_PX; i++) {
        if (isnan(margin[i])) continue; /* never reached by any bearing -- leave the pre-filled background */
        pixels[i] = SignalColor(margin[i], ceilingDb);
    }
    free(margin);
    return ceilingDb;
}

static void DrawPixelLine(Color *pixels, int w, int h, int x0, int y0, int x1, int y1, Color c) {
    int dx = abs(x1 - x0), dy = abs(y1 - y0);
    int steps = dx > dy ? dx : dy;
    if (steps == 0) {
        if (x0 >= 0 && x0 < w && y0 >= 0 && y0 < h) pixels[y0 * w + x0] = c;
        return;
    }
    for (int i = 0; i <= steps; i++) {
        int x = x0 + (x1 - x0) * i / steps;
        int y = y0 + (y1 - y0) * i / steps;
        if (x >= 0 && x < w && y >= 0 && y < h) pixels[y * w + x] = c;
    }
}

/* Best-effort overlay of real streets (see osm_roads.h) onto the already-
 * computed coverage raster, reusing the exact same tower-centered,
 * north-up local-meters projection ComputeViewshedRaster used to place
 * its own samples. Silently draws nothing if the fetch fails or the area
 * has no mapped roads -- this is a visual aid, never a reason to fail
 * the whole viewshed. */
static void DrawRoadOverlay(Color *pixels, double towerLat, double towerLon, double maxDistanceKm) {
    if (getenv("TV_NO_ROADS")) return;
    OsmRoadSet roads = { 0 };
    if (OsmRoadsFetch(towerLat, towerLon, maxDistanceKm, &roads) != 0) return;

    double metersPerPixel = (maxDistanceKm * 2000.0) / VIEWSHED_MAP_PX;
    Color roadColor = (Color){ 30, 30, 30, 255 };
    for (int i = 0; i < roads.way_count; i++) {
        OsmWay *w = &roads.ways[i];
        int prevPx = 0, prevPy = 0, havePrev = 0;
        for (int j = 0; j < w->count; j++) {
            double eastM, northM;
            DtmLonLatToLocalMeters(towerLat, towerLon, w->lat[j], w->lon[j], &eastM, &northM);
            int px = VIEWSHED_MAP_PX / 2 + (int)lround(eastM / metersPerPixel);
            int py = VIEWSHED_MAP_PX / 2 - (int)lround(northM / metersPerPixel);
            if (havePrev) DrawPixelLine(pixels, VIEWSHED_MAP_PX, VIEWSHED_MAP_PX, prevPx, prevPy, px, py, roadColor);
            prevPx = px; prevPy = py; havePrev = 1;
        }
    }
    OsmRoadsFree(&roads);
}

static int RunViewshedMode(double towerLat, double towerLon, double heightFt, double maxDistanceKm, const RfParams *rf, const char *screenshotPath) {
    EnsureCoverageOrWarn(towerLat, towerLon);
    for (int b = 0; b < 8; b++) {
        double lat, lon;
        destination_point(towerLat, towerLon, b * 45.0, maxDistanceKm, EARTH_RADIUS_KM, &lat, &lon);
        EnsureCoverageOrWarn(lat, lon);
    }

    double towerGround = DtmSampleMeters(&g_dtm, towerLat, towerLon, 0);
    if (towerGround == DTM_NODATA) {
        fprintf(stderr, "terrain_viewer: tower position (%.5f, %.5f) isn't covered by any loaded DTM file\n", towerLat, towerLon);
        return 1;
    }
    double towerElevM = towerGround + heightFt * FEET_TO_METERS;
    printf("Viewshed: tower (%.5f,%.5f) %.0fm ground + %.0fft = %.0fm, range %.1f km, receiver height %.1fm\n",
           towerLat, towerLon, towerGround, heightFt, towerElevM, maxDistanceKm, VIEWSHED_RECEIVER_HEIGHT_M);
    if (rf->enabled) {
        printf("RF model: %.1f MHz, EIRP %.1f dBm, sensitivity threshold %.1f dBm%s\n", rf->freqMHz, rf->eirpDbm, rf->sensDbm,
               rf->useItm ? " (Longley-Rice/ITM point-to-point)" : " (single-knife-edge diffraction)");
    }

    Color *pixels = malloc(sizeof(Color) * VIEWSHED_MAP_PX * VIEWSHED_MAP_PX);
    if (!pixels) { fprintf(stderr, "terrain_viewer: out of memory\n"); return 1; }
    for (int i = 0; i < VIEWSHED_MAP_PX * VIEWSHED_MAP_PX; i++) pixels[i] = (Color){ 224, 224, 224, 255 };

    printf("Computing viewshed raster...\n");
    double ceilingDb = ComputeViewshedRaster(pixels, towerLat, towerLon, towerElevM, maxDistanceKm, rf);
    if (rf->enabled) printf("Strongest margin in view: %.1f dB above sensitivity -- color scale set to match\n", ceilingDb);
    printf("Done.\n");

    DrawRoadOverlay(pixels, towerLat, towerLon, maxDistanceKm);

    SetConfigFlags(FLAG_WINDOW_RESIZABLE | FLAG_VSYNC_HINT);
    InitWindow(VIEWSHED_MAP_PX + 60, VIEWSHED_MAP_PX + 130, "DTM Viewshed");

    Image img = { pixels, VIEWSHED_MAP_PX, VIEWSHED_MAP_PX, 1, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8 };
    Texture2D tex = LoadTextureFromImage(img);
    free(pixels);

    int frameCount = 0;
    int screenshotFrame = getenv("TV_SCREENSHOT_FRAME") ? atoi(getenv("TV_SCREENSHOT_FRAME")) : 5;
    int marginL = 30, marginT = 70;

    while (!WindowShouldClose()) {
        frameCount++;
        BeginDrawing();
        ClearBackground(RAYWHITE);

        DrawTexture(tex, marginL, marginT, WHITE);
        DrawRectangleLines(marginL, marginT, VIEWSHED_MAP_PX, VIEWSHED_MAP_PX, GRAY);

        int cx = marginL + VIEWSHED_MAP_PX / 2, cy = marginT + VIEWSHED_MAP_PX / 2;
        for (int ring = 1; ring <= 4; ring++) {
            float r = (float)VIEWSHED_MAP_PX / 2.0f * ring / 4.0f;
            DrawCircleLines(cx, cy, r, (Color){ 120, 120, 120, 160 });
            DrawText(TextFormat("%.1fkm", maxDistanceKm * ring / 4.0), cx + 4, (int)(cy - r) + 2, 13, DARKGRAY);
        }
        DrawLine(cx, marginT, cx, marginT + VIEWSHED_MAP_PX, (Color){ 150, 150, 150, 90 });
        DrawLine(marginL, cy, marginL + VIEWSHED_MAP_PX, cy, (Color){ 150, 150, 150, 90 });
        DrawCircle(cx, cy, 5, BLACK);
        DrawCircle(cx, cy, 3, YELLOW);
        DrawText("N", cx - 4, marginT - 18, 16, DARKGRAY);

        DrawText(TextFormat("Tower (%.5f, %.5f)  %.0fm ground + %.0fft antenna = %.0fm    range %.1f km",
                             towerLat, towerLon, towerGround, heightFt, towerElevM, maxDistanceKm), marginL, 8, 15, DARKGRAY);
        if (rf->enabled) {
            DrawText(TextFormat("%.1f MHz, EIRP %.1f dBm  --  color = signal margin above %.0f dBm sensitivity threshold (gray = below it, no usable coverage)",
                                 rf->freqMHz, rf->eirpDbm, rf->sensDbm), marginL, 28, 14, DARKGRAY);
            if (rf->useItm) {
                DrawText(TextFormat("Longley-Rice/ITM point-to-point (climate/pol/ground constants below); %.1fm receiver height; NOT a substitute for a real RF survey", VIEWSHED_RECEIVER_HEIGHT_M),
                          marginL, 46, 13, GRAY);
            } else {
                DrawText(TextFormat("free-space path loss + single-knife-edge diffraction; %.1fm receiver height; NOT a substitute for a real RF survey", VIEWSHED_RECEIVER_HEIGHT_M),
                          marginL, 46, 13, GRAY);
            }

            int barX = marginL, barY = marginT + VIEWSHED_MAP_PX + 8, barW = 260, barH = 14;
            for (int i = 0; i < barW; i++) {
                double marginDb = ceilingDb * i / (barW - 1);
                DrawLine(barX + i, barY, barX + i, barY + barH, SignalColor(marginDb, ceilingDb));
            }
            DrawRectangleLines(barX, barY, barW, barH, DARKGRAY);
            DrawText(TextFormat("%.0f dBm", rf->sensDbm), barX, barY + barH + 2, 12, DARKGRAY);
            DrawText(TextFormat("%.0f dBm", rf->sensDbm + ceilingDb), barX + barW - 55, barY + barH + 2, 12, DARKGRAY);
            DrawRectangle(barX + barW + 20, barY, 14, 14, (Color){ 140, 140, 140, 255 });
            DrawText("no coverage", barX + barW + 40, barY, 14, DARKGRAY);
            DrawLine(barX + barW + 150, barY + 7, barX + barW + 164, barY + 7, (Color){ 30, 30, 30, 255 });
            DrawText("streets (OpenStreetMap)", barX + barW + 170, barY, 14, DARKGRAY);
        } else {
            DrawText("green = line-of-sight clear to tower   red = terrain-blocked   gray = outside loaded DTM coverage",
                      marginL, 28, 14, DARKGRAY);
            DrawText(TextFormat("assumes a %.1fm receiver height above ground at each point (not a signal-strength/path-loss model)", VIEWSHED_RECEIVER_HEIGHT_M),
                      marginL, 46, 13, GRAY);
            DrawRectangle(marginL, marginT + VIEWSHED_MAP_PX + 8, 14, 14, (Color){ 70, 160, 90, 255 });
            DrawText("visible", marginL + 20, marginT + VIEWSHED_MAP_PX + 8, 14, DARKGRAY);
            DrawRectangle(marginL + 90, marginT + VIEWSHED_MAP_PX + 8, 14, 14, (Color){ 195, 60, 50, 255 });
            DrawText("blocked", marginL + 110, marginT + VIEWSHED_MAP_PX + 8, 14, DARKGRAY);
            DrawLine(marginL + 180, marginT + VIEWSHED_MAP_PX + 15, marginL + 194, marginT + VIEWSHED_MAP_PX + 15, (Color){ 30, 30, 30, 255 });
            DrawText("streets (OpenStreetMap)", marginL + 200, marginT + VIEWSHED_MAP_PX + 8, 14, DARKGRAY);
        }

        EndDrawing();

        if (screenshotPath && frameCount == screenshotFrame) {
            TakeScreenshot(screenshotPath);
            break;
        }
    }

    UnloadTexture(tex);
    CloseWindow();
    return 0;
}

/* Loads whatever *.tif files are already sitting in the search dirs --
 * no network access, just whatever's local. Safe to call even if nothing's
 * there yet (g_dtm just stays empty); EnsureDtmCoverage below is what
 * actually fetches new tiles once a target location is known. */
static void ScanDtmDirectories(void) {
    const char *home = getenv("HOME");
    snprintf(g_dtmCacheDir, sizeof(g_dtmCacheDir), "%s/dtm_cache", home ? home : ".");
    mkdir(g_dtmCacheDir, 0755); /* ignore EEXIST -- fine if it's already there */

    const char *dirs[DTM_SEARCH_DIR_COUNT + 1];
    int dirCount = 0;
    for (int i = 0; i < DTM_SEARCH_DIR_COUNT; i++) dirs[dirCount++] = kDtmSearchDirs[i];
    dirs[dirCount++] = g_dtmCacheDir;

    for (int d = 0; d < dirCount; d++) {
        char pattern[700];
        snprintf(pattern, sizeof(pattern), "%s/*.tif", dirs[d]);
        glob_t g = { 0 };
        if (glob(pattern, 0, NULL, &g) == 0) {
            for (size_t i = 0; i < g.gl_pathc; i++) DtmSetAddFile(&g_dtm, g.gl_pathv[i]);
            globfree(&g);
        }
    }
}

/* Ensures a point is covered before it's sampled -- downloads a new tile
 * (Wyoming lidar first, USGS 3DEP fallback) if nothing loaded already
 * reaches there. TV_AUTO_DOWNLOAD=1 skips the interactive y/N confirmation
 * (used by our own TV_SCREENSHOT headless testing, and anyone scripting
 * this non-interactively). Not fatal on failure (e.g. the point is in the
 * ocean or outside the US) -- just warns and leaves that spot as NODATA. */
static void EnsureCoverageOrWarn(double lat, double lon) {
    int autoConfirm = getenv("TV_AUTO_DOWNLOAD") != NULL;
    if (EnsureDtmCoverage(&g_dtm, lat, lon, g_dtmCacheDir, autoConfirm) != 0) {
        fprintf(stderr, "terrain_viewer: no DTM coverage at (%.5f, %.5f) -- that area will show as NODATA\n", lat, lon);
    }
}

int main(int argc, char **argv) {
    if (argc >= 2 && strcmp(argv[1], "--profile") == 0) {
        static const char *usage = "usage: %s --profile lat1 lon1 [+height1] lat2 lon2 [+height2] [screenshot.png]\n"
                                    "  +height is an optional mast/tower height above ground at that point\n"
                                    "  (e.g. a cell tower's antenna height) -- must start with '+', or it's\n"
                                    "  read as the next lat/lon instead; 0 (ground level) if omitted. Add an\n"
                                    "  'm' suffix for meters (e.g. +30m) or 'f' for feet (e.g. +100f); no\n"
                                    "  suffix defaults to feet (e.g. +100 == +100f).\n";
        if (argc < 6) {
            fprintf(stderr, usage, argv[0]);
            return 1;
        }
        ScanDtmDirectories();

        int argi = 2;
        double lat1 = atof(argv[argi++]);
        double lon1 = atof(argv[argi++]);
        double heightFt1 = 0.0;
        if (argi < argc && argv[argi][0] == '+') heightFt1 = ParseHeightFeetArg(argv[argi++]);

        if (argi + 1 >= argc) {
            fprintf(stderr, usage, argv[0]);
            return 1;
        }
        double lat2 = atof(argv[argi++]);
        double lon2 = atof(argv[argi++]);
        double heightFt2 = 0.0;
        if (argi < argc && argv[argi][0] == '+') heightFt2 = ParseHeightFeetArg(argv[argi++]);

        const char *screenshot = (argi < argc) ? argv[argi] : getenv("TV_SCREENSHOT");
        return RunProfileMode(lat1, lon1, heightFt1, lat2, lon2, heightFt2, screenshot);
    }

    if (argc >= 2 && strcmp(argv[1], "--viewshed") == 0) {
        static const char *usage = "usage: %s --viewshed lat lon [+height] maxDistanceKm [options] [screenshot.png]\n"
                                    "  +height is an optional antenna height above ground at the tower --\n"
                                    "  must start with '+', or it's read as maxDistanceKm instead; 0 (ground\n"
                                    "  level) if omitted. Add an 'm' suffix for meters (e.g. +10m) or 'f' for\n"
                                    "  feet (e.g. +30f); no suffix defaults to feet (e.g. +30 == +30f).\n"
                                    "  maxDistanceKm is the coverage radius to map.\n"
                                    "  Options (all optional, any order, after maxDistanceKm):\n"
                                    "    --erp watts    transmitter ERP in watts -- ALSO enables the RF signal-\n"
                                    "                   strength model (default: off, pure geometric LOS map)\n"
                                    "    --freq MHz     carrier frequency in MHz (default: 869, cellular Band A)\n"
                                    "    --sens dBm     receiver sensitivity threshold (default: -100)\n"
                                    "    --itm          use the full Longley-Rice/ITM point-to-point model\n"
                                    "                   instead of the simpler single-knife-edge estimate\n"
                                    "                   (requires --erp; slower, more physically complete)\n"
                                    "    --climate NAME one of: equatorial, continental-subtropical,\n"
                                    "                   maritime-subtropical, desert, continental-temperate\n"
                                    "                   (default), maritime-temperate-land,\n"
                                    "                   maritime-temperate-sea -- ITM only\n"
                                    "    --pol h|v      polarization, horizontal or vertical (default: v) --\n"
                                    "                   ITM only\n"
                                    "    --permittivity E  relative permittivity of the ground (default: 15) --\n"
                                    "                   ITM only\n"
                                    "    --conductivity S  ground conductivity, S/m (default: 0.005) -- ITM only\n"
                                    "    --refractivity N  surface refractivity, N-units (default: 301) --\n"
                                    "                   ITM only\n";
        if (argc < 5) {
            fprintf(stderr, usage, argv[0]);
            return 1;
        }
        ScanDtmDirectories();

        int argi = 2;
        double lat = atof(argv[argi++]);
        double lon = atof(argv[argi++]);
        double heightFt = 0.0;
        if (argi < argc && argv[argi][0] == '+') heightFt = ParseHeightFeetArg(argv[argi++]);

        if (argi >= argc) {
            fprintf(stderr, usage, argv[0]);
            return 1;
        }
        double maxDistanceKm = atof(argv[argi++]);

        RfParams rf = { 0, 869.0, 0.0, -100.0, 0, CLIMATE__CONTINENTAL_TEMPERATE, POLARIZATION__VERTICAL, 15.0, 0.005, 301.0 };
        const char *screenshot = getenv("TV_SCREENSHOT");
        while (argi < argc) {
            if (strcmp(argv[argi], "--freq") == 0 && argi + 1 < argc) {
                rf.freqMHz = atof(argv[argi + 1]);
                argi += 2;
            } else if (strcmp(argv[argi], "--erp") == 0 && argi + 1 < argc) {
                double erpWatts = atof(argv[argi + 1]);
                rf.eirpDbm = 10.0 * log10(erpWatts * 1000.0) + 2.15; /* ERP (ref. dipole) -> EIRP (ref. isotropic), the FCC's usual "power" convention */
                rf.enabled = 1;
                argi += 2;
            } else if (strcmp(argv[argi], "--sens") == 0 && argi + 1 < argc) {
                rf.sensDbm = atof(argv[argi + 1]);
                argi += 2;
            } else if (strcmp(argv[argi], "--itm") == 0) {
                rf.useItm = 1;
                argi += 1;
            } else if (strcmp(argv[argi], "--climate") == 0 && argi + 1 < argc) {
                const char *name = argv[argi + 1];
                if (strcmp(name, "equatorial") == 0) rf.climate = CLIMATE__EQUATORIAL;
                else if (strcmp(name, "continental-subtropical") == 0) rf.climate = CLIMATE__CONTINENTAL_SUBTROPICAL;
                else if (strcmp(name, "maritime-subtropical") == 0) rf.climate = CLIMATE__MARITIME_SUBTROPICAL;
                else if (strcmp(name, "desert") == 0) rf.climate = CLIMATE__DESERT;
                else if (strcmp(name, "continental-temperate") == 0) rf.climate = CLIMATE__CONTINENTAL_TEMPERATE;
                else if (strcmp(name, "maritime-temperate-land") == 0) rf.climate = CLIMATE__MARITIME_TEMPERATE_OVER_LAND;
                else if (strcmp(name, "maritime-temperate-sea") == 0) rf.climate = CLIMATE__MARITIME_TEMPERATE_OVER_SEA;
                else { fprintf(stderr, "terrain_viewer: unknown --climate '%s'\n", name); fprintf(stderr, usage, argv[0]); return 1; }
                argi += 2;
            } else if (strcmp(argv[argi], "--pol") == 0 && argi + 1 < argc) {
                if (strcmp(argv[argi + 1], "h") == 0) rf.pol = POLARIZATION__HORIZONTAL;
                else if (strcmp(argv[argi + 1], "v") == 0) rf.pol = POLARIZATION__VERTICAL;
                else { fprintf(stderr, "terrain_viewer: --pol must be 'h' or 'v'\n"); fprintf(stderr, usage, argv[0]); return 1; }
                argi += 2;
            } else if (strcmp(argv[argi], "--permittivity") == 0 && argi + 1 < argc) {
                rf.epsilon = atof(argv[argi + 1]);
                argi += 2;
            } else if (strcmp(argv[argi], "--conductivity") == 0 && argi + 1 < argc) {
                rf.sigma = atof(argv[argi + 1]);
                argi += 2;
            } else if (strcmp(argv[argi], "--refractivity") == 0 && argi + 1 < argc) {
                rf.N0 = atof(argv[argi + 1]);
                argi += 2;
            } else {
                screenshot = argv[argi];
                argi++;
            }
        }
        if (rf.useItm && !rf.enabled) {
            fprintf(stderr, "terrain_viewer: --itm requires --erp (the RF model must be enabled)\n");
            fprintf(stderr, usage, argv[0]);
            return 1;
        }
        return RunViewshedMode(lat, lon, heightFt, maxDistanceKm, &rf, screenshot);
    }

    ScanDtmDirectories();

    double lat = kPlaces[0].lat, lon = kPlaces[0].lon;
    if (argc >= 2) {
        int matched = 0;
        for (int i = 0; i < PLACE_COUNT; i++) {
            if (strcmp(argv[1], kPlaces[i].name) == 0) { lat = kPlaces[i].lat; lon = kPlaces[i].lon; matched = 1; break; }
        }
        if (!matched) fprintf(stderr, "terrain_viewer: unknown place '%s', defaulting to cody (known: cody, meeteetse, bighorn)\n", argv[1]);
    }
    if (getenv("TV_LAT")) lat = atof(getenv("TV_LAT"));
    if (getenv("TV_LON")) lon = atof(getenv("TV_LON"));

    int result = RunWalkMode(lat, lon);
    DtmSetClose(&g_dtm);
    return result;
}
