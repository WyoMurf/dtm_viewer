#define _POSIX_C_SOURCE 200809L /* popen/pclose */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

#include "dtm_fetch.h"

/* --- Wyoming wyolidar: statewide 1m lidar, W{lon}N{lat}_Deg_Cog.tif,
 * confirmed by listing the real bucket (40 tiles, W104N041..W111N045).
 * Tile W{n} covers longitude [-n-0.98, -n+0.02] (the ~0.02deg buffer this
 * project's dtm.c already documents); N{m} covers latitude [m-0.02,
 * m+0.98] the same way -- both formulas verified against the three tiles
 * already in use (W109N044 = lon[-109.98,-108.98] lat[43.98,44.98]). */
static void WyTileName(double lat, double lon, char *out, size_t outsz) {
    int lonN = (int)floor(-lon + 0.02);
    int latN = (int)floor(lat + 0.02);
    snprintf(out, outsz, "W%03dN%03d", lonN, latN);
}
static void WyRemoteUrl(const char *tile, char *out, size_t outsz) {
    snprintf(out, outsz, "https://wyolidar.s3.arcc.uwyo.edu/COG_Mosaic/dtm/%s_Deg_Cog.tif", tile);
}
static void WyLocalName(const char *tile, char *out, size_t outsz) {
    snprintf(out, outsz, "%s_Deg_Cog.tif", tile);
}

/* --- USGS 3DEP: nationwide 1/3 arc-second (~10m) seamless DEM,
 * n{lat}w{-lon}, confirmed by downloading and inspecting a real tile
 * (classic TIFF, float32, NODATA=-999999 -- see dtm.c's float32 support).
 * Tile n{X}w{Y} covers latitude [X-1,X], longitude [-Y,-Y+1] (NW corner at
 * X degrees N, Y degrees W) -- clean degree-aligned, no buffer. */
static void UsgsTileName(double lat, double lon, char *out, size_t outsz) {
    int latN = (int)ceil(lat);
    int lonN = (int)ceil(-lon);
    snprintf(out, outsz, "n%dw%d", latN, lonN);
}
static void UsgsRemoteUrl(const char *tile, char *out, size_t outsz) {
    snprintf(out, outsz, "https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/13/TIFF/current/%s/USGS_13_%s.tif", tile, tile);
}
static void UsgsLocalName(const char *tile, char *out, size_t outsz) {
    snprintf(out, outsz, "USGS_13_%s.tif", tile);
}

typedef struct {
    const char *label;
    void (*tileName)(double lat, double lon, char *out, size_t outsz);
    void (*remoteUrl)(const char *tile, char *out, size_t outsz);
    void (*localName)(const char *tile, char *out, size_t outsz);
} DtmSource;

static const DtmSource kWySource = { "Wyoming wyolidar (~1m)", WyTileName, WyRemoteUrl, WyLocalName };
static const DtmSource kUsgsSource = { "USGS 3DEP (~10m)", UsgsTileName, UsgsRemoteUrl, UsgsLocalName };

static int FileExists(const char *path) {
    struct stat st;
    return stat(path, &st) == 0;
}

/* -1 if the file doesn't exist/can't be stat'd. */
static long FileSize(const char *path) {
    struct stat st;
    return stat(path, &st) == 0 ? (long)st.st_size : -1;
}

/* HEAD request via curl -- returns the HTTP status code (0 on a curl/
 * network failure), and via *outSize the Content-Length if present (-1 if
 * not). No auth, no request body, just enough to know whether the tile
 * exists and how big it is before committing to a download. */
static long HeadRequest(const char *url, long *outSize) {
    *outSize = -1;
    char cmd[1200];
    snprintf(cmd, sizeof(cmd), "curl -sI --max-time 15 '%s' 2>/dev/null", url);
    FILE *p = popen(cmd, "r");
    if (!p) return 0;

    char line[512];
    long httpCode = 0;
    int first = 1;
    while (fgets(line, sizeof(line), p)) {
        if (first) {
            sscanf(line, "HTTP/%*s %ld", &httpCode);
            first = 0;
            continue;
        }
        long v;
        if (sscanf(line, "Content-Length: %ld", &v) == 1 || sscanf(line, "content-length: %ld", &v) == 1) {
            *outSize = v;
        }
    }
    pclose(p);
    return httpCode;
}

static int DownloadFile(const char *url, const char *localPath) {
    char cmd[1200];
    /* -f: fail (nonzero exit) on HTTP errors rather than saving an error
     * page; the partial/temp file curl leaves behind on failure is fine to
     * ignore -- a retry just overwrites it. */
    snprintf(cmd, sizeof(cmd), "curl -sf --max-time 600 -o '%s' '%s'", localPath, url);
    return system(cmd) == 0 ? 0 : -1;
}

/* Prints what's about to be downloaded and, unless autoConfirm, asks y/N
 * on stdin. Returns 1 to proceed, 0 to decline. */
static int ConfirmDownload(const DtmSource *src, const char *tile, const char *url, long sizeBytes, int autoConfirm) {
    if (sizeBytes >= 0) {
        printf("dtm_fetch: %s tile %s found (%.0f MB)\n  %s\n", src->label, tile, sizeBytes / (1024.0 * 1024.0), url);
    } else {
        printf("dtm_fetch: %s tile %s found (size unknown)\n  %s\n", src->label, tile, url);
    }
    if (autoConfirm) {
        printf("dtm_fetch: TV_AUTO_DOWNLOAD set, downloading without prompting\n");
        return 1;
    }
    printf("Download this tile? [y/N] ");
    fflush(stdout);
    char resp[16];
    if (!fgets(resp, sizeof(resp), stdin)) return 0;
    return (resp[0] == 'y' || resp[0] == 'Y');
}

/* Tries one source end to end: local cache hit, or HEAD-then-maybe-
 * download. Returns 0 on success (a file is now open in `set` covering
 * the point), -1 if this source has no such tile, -2 if the user declined
 * the download prompt. */
static int TrySource(const DtmSource *src, DtmSet *set, double lat, double lon, const char *cacheDir, int autoConfirm) {
    char tile[32];
    src->tileName(lat, lon, tile, sizeof(tile));

    char localName[128], localPath[600];
    src->localName(tile, localName, sizeof(localName));
    snprintf(localPath, sizeof(localPath), "%s/%s", cacheDir, localName);

    if (FileExists(localPath)) {
        if (DtmSetAddFile(set, localPath) == 0) {
            printf("dtm_fetch: %s tile %s already cached at %s\n", src->label, tile, localPath);
            return 0;
        }
        /* Likely a truncated/corrupt download left over from an interrupted
         * run -- delete it and fall through to re-fetch rather than caching
         * the failure forever. */
        fprintf(stderr, "dtm_fetch: cached file %s failed to open (incomplete download?) -- deleting and re-fetching\n", localPath);
        remove(localPath);
    }

    char url[700];
    src->remoteUrl(tile, url, sizeof(url));
    long size;
    long httpCode = HeadRequest(url, &size);
    if (httpCode != 200) return -1; /* no such tile at this source -- not an error, just "try the next one" */

    if (!ConfirmDownload(src, tile, url, size, autoConfirm)) return -2;

    printf("dtm_fetch: downloading %s -> %s ...\n", url, localPath);
    if (DownloadFile(url, localPath) != 0) {
        fprintf(stderr, "dtm_fetch: download failed for %s\n", url);
        remove(localPath);
        return -1;
    }
    /* Belt-and-suspenders against a truncated file that still happens to
     * parse as a (geographically incomplete) valid TIFF -- e.g. one cut
     * short by a killed process rather than a clean network error curl
     * would have caught itself. Only checked when HEAD gave us a size. */
    if (size >= 0 && FileSize(localPath) != size) {
        fprintf(stderr, "dtm_fetch: downloaded file size mismatch for %s (expected %ld bytes, got %ld) -- discarding\n", localPath, size, FileSize(localPath));
        remove(localPath);
        return -1;
    }
    printf("dtm_fetch: done.\n");
    if (DtmSetAddFile(set, localPath) == 0) return 0;
    remove(localPath);
    return -1;
}

int EnsureDtmCoverage(DtmSet *set, double lat, double lon, const char *cacheDir, int autoConfirm) {
    if (DtmSetCovers(set, lat, lon)) return 0;

    int rc = TrySource(&kWySource, set, lat, lon, cacheDir, autoConfirm);
    if (rc == 0 || rc == -2) return rc;

    rc = TrySource(&kUsgsSource, set, lat, lon, cacheDir, autoConfirm);
    if (rc == 0 || rc == -2) return rc;

    fprintf(stderr, "dtm_fetch: no DTM coverage available at (%.5f, %.5f)\n", lat, lon);
    return -1;
}
