#define _POSIX_C_SOURCE 200809L /* popen/pclose */
#define _DEFAULT_SOURCE 1 /* M_PI, which _POSIX_C_SOURCE alone hides */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h> /* sleep() -- see FetchUrl's retry */

#include "osm_roads.h"

/* Bounded substring search -- strstr needs a NUL terminator, and we don't
 * want to plant one in the middle of the shared XML buffer just to bound
 * a search, so this scans an explicit [start,end) range instead. */
static const char *BoundedFind(const char *start, const char *end, const char *needle) {
    size_t nlen = strlen(needle);
    if (nlen == 0 || end < start) return NULL;
    for (const char *p = start; p + nlen <= end; p++) {
        if (memcmp(p, needle, nlen) == 0) return p;
    }
    return NULL;
}

/* atof/atol naturally stop at the closing '"', so no copying is needed --
 * just point them at the right spot in the original buffer. */
static int GetAttrDouble(const char *start, const char *end, const char *name, double *out) {
    char pat[32];
    snprintf(pat, sizeof(pat), "%s=\"", name);
    const char *p = BoundedFind(start, end, pat);
    if (!p) return 0;
    *out = atof(p + strlen(pat));
    return 1;
}
static int GetAttrLong(const char *start, const char *end, const char *name, long *out) {
    char pat[32];
    snprintf(pat, sizeof(pat), "%s=\"", name);
    const char *p = BoundedFind(start, end, pat);
    if (!p) return 0;
    *out = atol(p + strlen(pat));
    return 1;
}

typedef struct { long id; double lat, lon; } OsmNodeRec;

static int CompareNodeRec(const void *a, const void *b) {
    long ia = ((const OsmNodeRec *)a)->id, ib = ((const OsmNodeRec *)b)->id;
    return (ia > ib) - (ia < ib);
}

static int CompareLong(const void *a, const void *b) {
    long ia = *(const long *)a, ib = *(const long *)b;
    return (ia > ib) - (ia < ib);
}

/* Tracks which way ids have already been added to the accumulating
 * OsmRoadSet across multiple tile fetches (see the tiling comment on
 * OsmRoadsFetch) -- a way with nodes on both sides of a tile boundary
 * comes back in full from EVERY tile it touches, so without this a large
 * request would draw (and store) the same road several times over. */
typedef struct { long *ids; size_t count, cap; } LongSet;

static int LongSetContains(const LongSet *s, long id) {
    return bsearch(&id, s->ids, s->count, sizeof(long), CompareLong) != NULL;
}
static void LongSetAdd(LongSet *s, long id) {
    if (s->count == s->cap) { s->cap = s->cap ? s->cap * 2 : 256; s->ids = realloc(s->ids, s->cap * sizeof(long)); }
    s->ids[s->count++] = id;
    qsort(s->ids, s->count, sizeof(long), CompareLong); /* small/infrequent enough (once per way, not per lookup) that re-sorting here beats a real insert-sorted implementation for how little code it takes */
}

/* Hand-rolled, not a general XML parser -- OSM's /map API output is
 * regular enough (bounds, then every referenced node, then every way)
 * that two linear scans (collect nodes, then resolve each way's <nd
 * ref>s against them) is simpler and lighter than pulling in a real XML
 * library for this one call site, matching how dtm.c already hand-parses
 * GeoTIFF tags instead of depending on more of libtiff than TIFFOpen.
 * Appends into `out`/`seenWays` rather than resetting them, so multiple
 * tile fetches can accumulate into one combined result. */
static void ParseOsmXml(const char *buf, size_t len, OsmRoadSet *out, LongSet *seenWays) {
    const char *end = buf + len;

    size_t nodeCap = 4096, nodeCount = 0;
    OsmNodeRec *nodes = malloc(nodeCap * sizeof(OsmNodeRec));
    const char *p = buf;
    while ((p = BoundedFind(p, end, "<node ")) != NULL) {
        const char *tagEnd = BoundedFind(p, end, ">");
        if (!tagEnd) break;
        long id; double lat, lon;
        if (GetAttrLong(p, tagEnd, "id", &id) && GetAttrDouble(p, tagEnd, "lat", &lat) && GetAttrDouble(p, tagEnd, "lon", &lon)) {
            if (nodeCount == nodeCap) { nodeCap *= 2; nodes = realloc(nodes, nodeCap * sizeof(OsmNodeRec)); }
            nodes[nodeCount].id = id; nodes[nodeCount].lat = lat; nodes[nodeCount].lon = lon;
            nodeCount++;
        }
        p = tagEnd + 1;
    }
    qsort(nodes, nodeCount, sizeof(OsmNodeRec), CompareNodeRec);

    size_t wayCount = (size_t)out->way_count;
    size_t wayCap = wayCount > 0 ? wayCount : 256;
    while (wayCap < wayCount) wayCap *= 2;
    OsmWay *ways = out->ways ? realloc(out->ways, wayCap * sizeof(OsmWay)) : malloc(wayCap * sizeof(OsmWay));
    p = buf;
    while ((p = BoundedFind(p, end, "<way ")) != NULL) {
        const char *openEnd = BoundedFind(p, end, ">");
        if (!openEnd) break;
        int selfClosed = (openEnd > p && *(openEnd - 1) == '/');
        const char *blockEnd = openEnd;
        if (!selfClosed) {
            const char *closeTag = BoundedFind(openEnd, end, "</way>");
            if (!closeTag) break;
            blockEnd = closeTag;
        }

        long wayId = 0;
        int alreadySeen = selfClosed || !GetAttrLong(p, openEnd, "id", &wayId) || LongSetContains(seenWays, wayId);

        if (!alreadySeen && BoundedFind(openEnd, blockEnd, "<tag k=\"highway\"")) {
            size_t refCap = 16, refCount = 0;
            double *lats = malloc(refCap * sizeof(double));
            double *lons = malloc(refCap * sizeof(double));
            const char *q = openEnd, *nd;
            while ((nd = BoundedFind(q, blockEnd, "<nd ref=\"")) != NULL) {
                const char *ndTagEnd = BoundedFind(nd, blockEnd, ">");
                if (!ndTagEnd) break;
                long ref;
                if (GetAttrLong(nd, ndTagEnd, "ref", &ref)) {
                    OsmNodeRec key = { ref, 0, 0 };
                    OsmNodeRec *found = bsearch(&key, nodes, nodeCount, sizeof(OsmNodeRec), CompareNodeRec);
                    if (found) {
                        if (refCount == refCap) { refCap *= 2; lats = realloc(lats, refCap * sizeof(double)); lons = realloc(lons, refCap * sizeof(double)); }
                        lats[refCount] = found->lat; lons[refCount] = found->lon; refCount++;
                    }
                }
                q = ndTagEnd + 1;
            }
            if (refCount >= 2) {
                if (wayCount == wayCap) { wayCap *= 2; ways = realloc(ways, wayCap * sizeof(OsmWay)); }
                ways[wayCount].lat = lats; ways[wayCount].lon = lons; ways[wayCount].count = (int)refCount;
                wayCount++;
                LongSetAdd(seenWays, wayId);
            } else {
                free(lats); free(lons);
            }
        }
        p = selfClosed ? (openEnd + 1) : (blockEnd + 6); /* 6 == strlen("</way>") */
    }

    free(nodes);
    out->ways = ways;
    out->way_count = (int)wayCount;
}

static char *FetchUrlOnce(const char *url, size_t *outLen) {
    char cmd[900];
    snprintf(cmd, sizeof(cmd), "curl -sf --max-time 30 '%s'", url);
    FILE *p = popen(cmd, "r");
    if (!p) return NULL;
    size_t cap = 1 << 20, len = 0;
    char *buf = malloc(cap);
    size_t n;
    while ((n = fread(buf + len, 1, cap - len, p)) > 0) {
        len += n;
        if (len == cap) { cap *= 2; buf = realloc(buf, cap); }
    }
    int rc = pclose(p);
    if (rc != 0 || len == 0) { free(buf); return NULL; }
    if (outLen) *outLen = len;
    return buf;
}

/* OSM's live API is the production editing backend, not a CDN built for
 * heavy polling -- it can hand back a brief failure (rate limiting, a
 * momentary hiccup) that a few seconds later just works again. Confirmed
 * in practice: a real run failed outright on the first attempt, and the
 * exact same URL succeeded (HTTP 200) moments later with no code change
 * at all. Three tries with a short backoff turns that kind of blip into
 * "the overlay took an extra couple seconds" instead of "no roads". */
static char *FetchUrl(const char *url, size_t *outLen) {
    for (int attempt = 1; attempt <= 3; attempt++) {
        char *result = FetchUrlOnce(url, outLen);
        if (result) return result;
        if (attempt < 3) {
            fprintf(stderr, "osm_roads: fetch attempt %d failed, retrying...\n", attempt);
            sleep(attempt); /* 1s, then 2s */
        }
    }
    return NULL;
}

/* OSM's /map API refuses any request over 0.25 square degrees (confirmed
 * directly: a real 25km-radius request got a real HTTP 400, "The maximum
 * bbox size is 0.250000") -- roughly a 20-22km viewshed radius at
 * Wyoming's latitude. TILE_TARGET_DEG2 stays comfortably under that so a
 * bigger request just gets split into an N x N grid of sub-boxes, each
 * fetched and parsed (with LongSet catching ways that straddle a tile
 * boundary and would otherwise come back -- in full -- from every tile
 * they touch) into one combined OsmRoadSet, rather than refusing to draw
 * roads at all past ~20km. */
#define TILE_TARGET_DEG2 0.20

static int FetchOneTile(double minLon, double minLat, double maxLon, double maxLat, OsmRoadSet *out, LongSet *seenWays) {
    char url[400];
    snprintf(url, sizeof(url), "https://api.openstreetmap.org/api/0.6/map?bbox=%.6f,%.6f,%.6f,%.6f", minLon, minLat, maxLon, maxLat);
    size_t len = 0;
    char *xml = FetchUrl(url, &len);
    if (!xml) return -1;
    ParseOsmXml(xml, len, out, seenWays);
    free(xml);
    return 0;
}

int OsmRoadsFetch(double centerLat, double centerLon, double radiusKm, OsmRoadSet *out) {
    out->ways = NULL;
    out->way_count = 0;
    LongSet seenWays = { 0 };

    /* +5% pad so the requested circle is fully inside the fetched box. */
    double dLat = (radiusKm / 111.32) * 1.05;
    double dLon = (radiusKm / (111.32 * cos(centerLat * M_PI / 180.0))) * 1.05;
    double minLat = centerLat - dLat, maxLat = centerLat + dLat;
    double minLon = centerLon - dLon, maxLon = centerLon + dLon;
    double area = (maxLat - minLat) * (maxLon - minLon);

    int n = (int)ceil(sqrt(area / TILE_TARGET_DEG2));
    if (n < 1) n = 1;
    if (n > 6) { /* 6x6 already covers an ~80km radius -- past that, cap it rather than firing off dozens of sequential requests */
        fprintf(stderr, "osm_roads: %.1fkm radius is huge -- capping the road overlay to a 6x6 tile grid (partial coverage near the edges)\n", radiusKm);
        n = 6;
    }

    if (n > 1) printf("osm_roads: fetching street data for the overlay (%dx%d tiles, area too big for one request)...\n", n, n);
    else printf("osm_roads: fetching street data for the overlay...\n");

    int okCount = 0;
    for (int iy = 0; iy < n; iy++) {
        double tMinLat = minLat + (maxLat - minLat) * iy / n;
        double tMaxLat = minLat + (maxLat - minLat) * (iy + 1) / n;
        for (int ix = 0; ix < n; ix++) {
            double tMinLon = minLon + (maxLon - minLon) * ix / n;
            double tMaxLon = minLon + (maxLon - minLon) * (ix + 1) / n;
            if (FetchOneTile(tMinLon, tMinLat, tMaxLon, tMaxLat, out, &seenWays) == 0) okCount++;
        }
    }

    free(seenWays.ids);

    if (okCount == 0) {
        fprintf(stderr, "osm_roads: all %d tile fetch(es) failed -- skipping road overlay\n", n * n);
        return -1;
    }
    if (okCount < n * n) {
        fprintf(stderr, "osm_roads: %d of %d tiles failed to fetch -- showing partial road coverage\n", n * n - okCount, n * n);
    }
    printf("osm_roads: %d road segments found\n", out->way_count);
    return 0;
}

void OsmRoadsFree(OsmRoadSet *set) {
    for (int i = 0; i < set->way_count; i++) {
        free(set->ways[i].lat);
        free(set->ways[i].lon);
    }
    free(set->ways);
    set->ways = NULL;
    set->way_count = 0;
}
