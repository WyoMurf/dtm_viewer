#ifndef OSM_ROADS_H
#define OSM_ROADS_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Road overlay data for terrain_viewer's --viewshed mode: pulls real
 * street geometry from OpenStreetMap's public API (no auth, no key --
 * https://api.openstreetmap.org/api/0.6/map?bbox=... , confirmed
 * directly against real data before this was written, same "genuine
 * public source, verified not assumed" approach as dtm_fetch.h) so a
 * viewshed's coverage colors can be drawn under the streets they'd
 * actually affect, for whatever community happens to fall inside the
 * mapped circle.
 */

typedef struct {
    double *lat;
    double *lon;
    int count;
} OsmWay;

typedef struct {
    OsmWay *ways;
    int way_count;
} OsmRoadSet;

/*
 * Fetches every highway=* way inside a padded bounding box around a
 * circle of `radiusKm` centered on (centerLat, centerLon) and fills
 * `out`. Returns 0 on success (out->way_count may legitimately be 0 --
 * open country has no roads) or -1 if the fetch/parse failed or the
 * requested area exceeds the API's bbox size limit (~0.25 square
 * degrees) -- either way this is a best-effort overlay, never fatal to
 * the caller, so -1 just means "skip drawing roads this run".
 */
int OsmRoadsFetch(double centerLat, double centerLon, double radiusKm, OsmRoadSet *out);

/* Frees everything OsmRoadsFetch allocated into `set`. */
void OsmRoadsFree(OsmRoadSet *set);

#ifdef __cplusplus
}
#endif

#endif /* OSM_ROADS_H */
