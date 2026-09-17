#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "raylib.h"
#include "rlgl.h"
#include "raymath.h"
#include "dtm.h"
#include "geo_utils.h"

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

/* Default DTM file locations on this specific machine (100.74.88.66) --
 * see README-terrain.md. Overridable via argv if this ever runs elsewhere. */
static const char *kDefaultDtmPaths[] = {
    "/home/murf/wyodem/lidar/W109N044_Deg_Cog.tif", /* Cody */
    "/home/murf/wyodem/lidar/W108N044_Deg_Cog.tif", /* Meeteetse */
    "/home/murf/wyodem/lidar/W107N044_Deg_Cog.tif", /* Bighorn */
};
#define DEFAULT_DTM_COUNT (int)(sizeof(kDefaultDtmPaths) / sizeof(kDefaultDtmPaths[0]))

static DtmSet g_dtm = { 0 };

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
 * see README-terrain.md's known-limitations note -- just one square,
 * rebuilt whenever the camera has wandered far enough from its center. --- */
#define TERRAIN_HALF_SIZE_M 900.0f   /* window covers a ~1.8km square */
#define TERRAIN_RESOLUTION 240       /* 241x241 vertices = 58,081, under raylib's 65536-per-mesh index limit */
#define TERRAIN_REBUILD_THRESHOLD_M 300.0f /* rebuild once the camera is this far from the window's center */

typedef struct {
    Mesh mesh;
    Material material;
    int built;
    double centerLat, centerLon; /* the (lat,lon) this window's local (0,0) corresponds to */
    float centerX, centerZ;      /* that same center, in the walk's own local-meters frame */
} TerrainWindow;

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

static void BuildTerrainWindow(TerrainWindow *tw, double originLat, double originLon,
                                double centerLat, double centerLon, double originElevation) {
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

    float *heights = malloc(sizeof(float) * (size_t)vertexCount);
    float minH = 1e9f, maxH = -1e9f;
    double step = (2.0 * TERRAIN_HALF_SIZE_M) / res;

    for (int j = 0; j < vertsPerSide; j++) {
        double localNorth = centerNorthM + (-TERRAIN_HALF_SIZE_M + j * step);
        for (int i = 0; i < vertsPerSide; i++) {
            double localEast = centerEastM + (-TERRAIN_HALF_SIZE_M + i * step);
            double lat, lon;
            DtmLocalMetersToLonLat(originLat, originLon, localEast, localNorth, &lat, &lon);
            double m = DtmSampleMeters(&g_dtm, lat, lon, 0);
            float h = (m == DTM_NODATA) ? 0.0f : (float)(m - originElevation);
            heights[j * vertsPerSide + i] = h;
            if (h < minH) minH = h;
            if (h > maxH) maxH = h;
        }
    }
    if (maxH <= minH) maxH = minH + 1.0f; /* avoid divide-by-zero on a perfectly flat window */

    int v = 0;
    for (int j = 0; j < vertsPerSide; j++) {
        float z = (float)(-TERRAIN_HALF_SIZE_M + j * step); /* south is +Z, see file header comment */
        for (int i = 0; i < vertsPerSide; i++) {
            float x = (float)(-TERRAIN_HALF_SIZE_M + i * step);
            float h = heights[j * vertsPerSide + i];
            mesh.vertices[v * 3 + 0] = x;
            mesh.vertices[v * 3 + 1] = h;
            mesh.vertices[v * 3 + 2] = z;

            /* Normal from central differences, clamped to the window edge
             * (one-sided difference there) -- good enough for lighting a
             * heightfield, no need for the exact analytic surface normal. */
            int iL = i > 0 ? i - 1 : i, iR = i < res ? i + 1 : i;
            int jU = j > 0 ? j - 1 : j, jD = j < res ? j + 1 : j;
            float hL = heights[j * vertsPerSide + iL], hR = heights[j * vertsPerSide + iR];
            float hU = heights[jU * vertsPerSide + i], hD = heights[jD * vertsPerSide + i];
            Vector3 n = Vector3Normalize((Vector3){ hL - hR, 2.0f * (float)step, hU - hD });
            mesh.normals[v * 3 + 0] = n.x;
            mesh.normals[v * 3 + 1] = n.y;
            mesh.normals[v * 3 + 2] = n.z;

            Color c = ElevationColor((h - minH) / (maxH - minH));
            mesh.colors[v * 4 + 0] = c.r;
            mesh.colors[v * 4 + 1] = c.g;
            mesh.colors[v * 4 + 2] = c.b;
            mesh.colors[v * 4 + 3] = 255;
            v++;
        }
    }
    free(heights);

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
}

/* --- Walk mode --- */

typedef struct {
    Vector3 position; /* world meters, relative to walk origin */
    float yaw, pitch;  /* degrees */
} WalkCamera;

#define EYE_HEIGHT_M 1.7f
#define WALK_SPEED_MPS 4.0f
#define SPRINT_MULT 3.0f
#define MOUSE_SENSITIVITY 0.12f

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

    TerrainWindow terrain = { 0 };
    BuildTerrainWindow(&terrain, originLat, originLon, originLat, originLon, originElevation);

    Camera3D camera = { 0 };
    camera.fovy = 70.0f;
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

        float dx = wc.position.x - terrain.centerX, dz = wc.position.z - terrain.centerZ;
        if (sqrtf(dx * dx + dz * dz) > TERRAIN_REBUILD_THRESHOLD_M) {
            double lat, lon;
            DtmLocalMetersToLonLat(originLat, originLon, wc.position.x, -wc.position.z, &lat, &lon);
            BuildTerrainWindow(&terrain, originLat, originLon, lat, lon, originElevation);
        }

        camera.position = wc.position;
        Vector3 forward = WalkForward(&wc);
        camera.target = Vector3Add(wc.position, forward);

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
        }
        DrawText("WASD move, Shift sprint, mouse look, Esc to quit", 10, GetScreenHeight() - 30, 14, SKYBLUE);
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

static int RunProfileMode(double lat1, double lon1, double lat2, double lon2, const char *screenshotPath) {
    double totalKm = haversine_distance(lat1, lon1, lat2, lon2, EARTH_RADIUS_KM);
    double bearing = initial_bearing(lat1, lon1, lat2, lon2);
    printf("Profile: (%.5f,%.5f) -> (%.5f,%.5f), %.3f km, initial bearing %.1f deg\n",
           lat1, lon1, lat2, lon2, totalKm, bearing);

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
    double elevValid1 = (elev1 == DTM_NODATA) ? elevM[0] : elev1;
    double elevValid2 = (elev2 == DTM_NODATA) ? elevM[PROFILE_SAMPLES - 1] : elev2;
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

        DrawCircle(SX(0), SY(elevValid1), 5, BLUE);
        DrawCircle(SX(totalKm), SY(elevValid2), 5, BLUE);
        DrawText(TextFormat("Total distance: %.3f km   (red = terrain blocking line-of-sight; thin red line = direct sightline incl. curvature+refraction)", totalKm), marginL, 8, 15, DARKGRAY);
        DrawText(TextFormat("A (%.4f,%.4f) %.0fm", lat1, lon1, elevValid1), marginL, marginT - 22, 16, BLUE);
        DrawText(TextFormat("B (%.4f,%.4f) %.0fm", lat2, lon2, elevValid2), marginL + plotW - 260, marginT - 22, 16, BLUE);

        EndDrawing();

        if (screenshotPath && frameCount == screenshotFrame) {
            TakeScreenshot(screenshotPath);
            break;
        }
    }

    CloseWindow();
    return 0;
}

static void LoadDefaultDtmSet(void) {
    for (int i = 0; i < DEFAULT_DTM_COUNT; i++) {
        if (DtmSetAddFile(&g_dtm, kDefaultDtmPaths[i]) != 0) {
            fprintf(stderr, "terrain_viewer: failed to open %s -- see README-terrain.md\n", kDefaultDtmPaths[i]);
        }
    }
    if (g_dtm.file_count == 0) {
        fprintf(stderr, "terrain_viewer: no DTM files loaded, nothing to show\n");
        exit(1);
    }
}

int main(int argc, char **argv) {
    if (argc >= 2 && strcmp(argv[1], "--profile") == 0) {
        if (argc < 6) {
            fprintf(stderr, "usage: %s --profile lat1 lon1 lat2 lon2 [screenshot.png]\n", argv[0]);
            return 1;
        }
        LoadDefaultDtmSet();
        double lat1 = atof(argv[2]), lon1 = atof(argv[3]), lat2 = atof(argv[4]), lon2 = atof(argv[5]);
        const char *screenshot = argc >= 7 ? argv[6] : getenv("TV_SCREENSHOT");
        return RunProfileMode(lat1, lon1, lat2, lon2, screenshot);
    }

    LoadDefaultDtmSet();

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
