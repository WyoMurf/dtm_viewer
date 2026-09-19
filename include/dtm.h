#ifndef DTM_H
#define DTM_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Reader for the WyoGISC/wyolidar-style DTM GeoTIFF COGs: EPSG:4326
 * (geographic lat/lon, not projected), single-band 16-bit unsigned int
 * elevation in whole meters, tiled, with a pyramid of overview IFDs (each
 * half the width/height of the one before). Tiles are decoded on demand via
 * libtiff and cached -- opening one of these files never reads the whole
 * multi-GB image into memory, only the tiles a given DtmSampleMeters() call
 * actually touches.
 *
 * This is deliberately narrow: it assumes the specific layout these DTM
 * files have (checked by DtmSetAddFile, which rejects anything that doesn't
 * match) rather than being a general-purpose GeoTIFF reader.
 */

#define DTM_NODATA (-32768.0) /* returned when a point is outside every open file's coverage, or lands where the source raster itself is NODATA (65535 in these files) */
#define DTM_MAX_FILES 32 /* was 8 -- too easy to exceed once ScanDtmDirectories is loading a growing wyodem/lidar plus dtm_cache instead of 3 hardcoded tiles */

typedef struct DtmFile DtmFile; /* opaque; see dtm.c */

typedef struct {
    DtmFile *files[DTM_MAX_FILES];
    int file_count;
} DtmSet;

/*
 * Opens `path` and adds it to `set`. Returns 0 on success, -1 on failure
 * (file won't open, isn't a recognized DTM layout, or `set` is already at
 * DTM_MAX_FILES) -- a message is printed to stderr explaining which.
 */
int DtmSetAddFile(DtmSet *set, const char *path);

/* Closes every file in the set (libtiff handles + cached tile buffers). */
void DtmSetClose(DtmSet *set);

/*
 * Bilinearly-interpolated elevation in meters at (lat, lon) in degrees,
 * sampled from whichever file in the set covers that point (first match
 * wins, checked in the order files were added -- fine as long as the DTM
 * tiles you add don't overlap, which the wyolidar 1-degree tiles don't).
 * `level` selects a pyramid level: 0 is full resolution (~1m/pixel),
 * each level up is roughly half that resolution again, clamped to
 * whatever levels the covering file actually has. Use level 0 for
 * close-up work (the ground walkthrough) and a higher level when sampling
 * over tens of kilometers is all that's needed (e.g. the elevation
 * profile tool) -- it's both faster and enough precision at that scale.
 *
 * Returns DTM_NODATA if no open file covers (lat, lon), or if every pixel
 * bilinear interpolation would blend at that point is the source
 * raster's own NODATA value.
 */
double DtmSampleMeters(DtmSet *set, double lat, double lon, int level);

/*
 * Whether some already-open file in `set` covers (lat, lon) geographically
 * (bounding box only -- doesn't check for a NODATA void within it). Meant
 * for callers deciding whether a *new* tile needs to be located and opened
 * before sampling at a point (see dtm_fetch.h's EnsureDtmCoverage) -- much
 * cheaper than a full DtmSampleMeters call when all you need is yes/no.
 */
int DtmSetCovers(const DtmSet *set, double lat, double lon);

/*
 * Flat-earth tangent-plane projection centered at (originLat, originLon):
 * converts to/from meters east/north of that origin. Accurate to a few cm
 * per km at these latitudes -- meant for placing terrain within a
 * several-km walk radius in a 3D engine's world space, not for anything
 * larger (no earth curvature or ellipsoid correction is applied here).
 */
void DtmLonLatToLocalMeters(double originLat, double originLon, double lat, double lon, double *outEastM, double *outNorthM);
void DtmLocalMetersToLonLat(double originLat, double originLon, double eastM, double northM, double *outLat, double *outLon);

#ifdef __cplusplus
}
#endif

#endif /* DTM_H */
