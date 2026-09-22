#ifndef DTM_FETCH_H
#define DTM_FETCH_H

#include "dtm.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Automatic DTM tile acquisition: given a point, figures out which real
 * data source covers it and gets a file loaded into a DtmSet, downloading
 * it first if necessary. Three sources, tried in order (best resolution
 * first):
 *
 *  1. Wyoming's wyolidar project (wyolidar.s3.arcc.uwyo.edu) -- ~1m lidar,
 *     statewide (40 tiles, W104N041..W111N045), same format/naming this
 *     project's own manually-downloaded tiles already use. Formula-based
 *     tile naming (WyTileName/WyRemoteUrl) -- no network round trip
 *     needed just to locate the right tile.
 *  2. USGS 3DEP's 1-meter DEM (also on prd-tnm.s3.amazonaws.com, under
 *     StagedProducts/Elevation/1m/) -- matches Wyoming's own resolution
 *     wherever it's flown, which by now is most of the country, but tiles
 *     are organized by named lidar-acquisition "Projects" on a UTM grid
 *     rather than a clean lat/lon formula, so locating the right one
 *     needs a real lookup against USGS's TNM Access product-search API
 *     (tnmaccess.nationalmap.gov -- the same API their own download tools
 *     use) rather than a filename formula.
 *  3. USGS 3DEP's nationwide 1/3 arc-second seamless DEM (also on
 *     prd-tnm.s3.amazonaws.com) -- ~10m, float32 samples, nationally
 *     complete with no gaps, the final fallback wherever neither of the
 *     above has coverage.
 *
 * All three are genuine public, unauthenticated sources with either
 * predictable per-tile URLs or a real public lookup API -- confirmed
 * directly (not assumed) by listing/downloading/querying real tiles from
 * each before this was written; see README.md's DATA SOURCES section for
 * the URL patterns and naming formulas.
 */

/*
 * Ensures `set` has a file covering (lat, lon). If one's already loaded
 * (checked via DtmSetCovers), returns 0 immediately with no I/O at all.
 * Otherwise: checks `cacheDir` for an already-downloaded matching file: if
 * present, just opens it. Otherwise tries Wyoming's tile, then USGS's 1m
 * tile (via a TNM API lookup), then USGS's 1/3 arc-second tile, each via
 * a HEAD request (or, for the 1m source, the lookup response itself)
 * first -- if a tile exists there, its name/source/size are printed and,
 * unless `autoConfirm` is set, a y/N prompt on stdin confirms before
 * actually downloading (curl -o) into `cacheDir` and opening the result.
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
