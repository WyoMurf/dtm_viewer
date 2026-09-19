#ifndef DTM_FETCH_H
#define DTM_FETCH_H

#include "dtm.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Automatic DTM tile acquisition: given a point, figures out which real
 * data source covers it and gets a file loaded into a DtmSet, downloading
 * it first if necessary. Two sources, tried in order:
 *
 *  1. Wyoming's wyolidar project (wyolidar.s3.arcc.uwyo.edu) -- ~1m lidar,
 *     statewide (40 tiles, W104N041..W111N045), same format/naming this
 *     project's own manually-downloaded tiles already use.
 *  2. USGS 3DEP's nationwide 1/3 arc-second seamless DEM
 *     (prd-tnm.s3.amazonaws.com) -- ~10m, float32 samples, covers the
 *     whole US as a fallback wherever Wyoming's higher-res data doesn't
 *     reach.
 *
 * Both are genuine public, unauthenticated S3-compatible buckets with
 * predictable per-tile URLs -- confirmed directly (not assumed) by
 * listing/downloading real tiles from each before this was written; see
 * README.md's DATA SOURCES section for the URL patterns and the naming
 * formulas.
 */

/*
 * Ensures `set` has a file covering (lat, lon). If one's already loaded
 * (checked via DtmSetCovers), returns 0 immediately with no I/O at all.
 * Otherwise: checks `cacheDir` for an already-downloaded matching file: if
 * present, just opens it. Otherwise tries Wyoming's tile, then USGS's,
 * each via a HEAD request first -- if a tile exists there, its name/
 * source/size are printed and, unless `autoConfirm` is set, a y/N prompt
 * on stdin confirms before actually downloading (curl -o) into `cacheDir`
 * and opening the result.
 *
 * Returns 0 on success (already covered, found locally, or freshly
 * downloaded and opened), -1 if no coverage is available there (e.g.
 * ocean, outside the US) or a download/open attempt failed, and -2 if the
 * user declined the download prompt.
 */
int EnsureDtmCoverage(DtmSet *set, double lat, double lon, const char *cacheDir, int autoConfirm);

#ifdef __cplusplus
}
#endif

#endif /* DTM_FETCH_H */
