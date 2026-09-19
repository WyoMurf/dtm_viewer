#define _POSIX_C_SOURCE 200809L /* popen/pclose */
#define _DEFAULT_SOURCE 1 /* M_PI, which _POSIX_C_SOURCE alone hides */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

/* Hand-rolled, not a general XML parser -- OSM's /map API output is
 * regular enough (bounds, then every referenced node, then every way)
 * that two linear scans (collect nodes, then resolve each way's <nd
 * ref>s against them) is simpler and lighter than pulling in a real XML
 * library for this one call site, matching how dtm.c already hand-parses
 * GeoTIFF tags instead of depending on more of libtiff than TIFFOpen. */
static void ParseOsmXml(const char *buf, size_t len, OsmRoadSet *out) {
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

    size_t wayCap = 256, wayCount = 0;
    OsmWay *ways = malloc(wayCap * sizeof(OsmWay));
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

        if (!selfClosed && BoundedFind(openEnd, blockEnd, "<tag k=\"highway\"")) {
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

static char *FetchUrl(const char *url, size_t *outLen) {
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

int OsmRoadsFetch(double centerLat, double centerLon, double radiusKm, OsmRoadSet *out) {
    out->ways = NULL;
    out->way_count = 0;

    /* +5% pad so the requested circle is fully inside the fetched box. */
    double dLat = (radiusKm / 111.32) * 1.05;
    double dLon = (radiusKm / (111.32 * cos(centerLat * M_PI / 180.0))) * 1.05;
    double minLat = centerLat - dLat, maxLat = centerLat + dLat;
    double minLon = centerLon - dLon, maxLon = centerLon + dLon;

    if ((maxLat - minLat) * (maxLon - minLon) > 0.24) {
        fprintf(stderr, "osm_roads: %.1fkm radius needs a bounding box too big for OSM's /map API (limit ~0.25 sq deg) -- skipping road overlay\n", radiusKm);
        return -1;
    }

    char url[400];
    snprintf(url, sizeof(url), "https://api.openstreetmap.org/api/0.6/map?bbox=%.6f,%.6f,%.6f,%.6f", minLon, minLat, maxLon, maxLat);

    printf("osm_roads: fetching street data for the overlay...\n");
    size_t len = 0;
    char *xml = FetchUrl(url, &len);
    if (!xml) {
        fprintf(stderr, "osm_roads: fetch failed -- skipping road overlay\n");
        return -1;
    }

    ParseOsmXml(xml, len, out);
    free(xml);
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
