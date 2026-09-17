#ifndef GEO_UTILS_H
#define GEO_UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Small geo/angle math utilities used by terrain_viewer's --profile mode.
 * Spherical-Earth only (no ellipsoid model) -- exact enough for plotting an
 * elevation profile, not meant as a general geodesy library. (This is a
 * trimmed-down, self-contained copy of a couple of functions from a larger
 * geo_utils.h in a separate project of mine, github.com/WyoMurf/kdtree --
 * kept here rather than as a dependency since this project has nothing else
 * to do with that one.)
 *
 * All latitude/longitude angle inputs and outputs are in degrees.
 */

#define EARTH_RADIUS_KM 6371.0

/*
 * Fast, approximate great-circle distance between two lat/lon points on a
 * perfect sphere of the given radius, using the haversine formula. Units of
 * the return value match the units of `radius` (e.g. pass EARTH_RADIUS_KM
 * for kilometers on Earth modeled as a sphere).
 */
double haversine_distance(double lat1, double lon1, double lat2, double lon2, double radius);

/*
 * Initial bearing (forward azimuth, degrees clockwise from true north, in
 * [0, 360)) of the great-circle path from (lat1,lon1) to (lat2,lon2), on a
 * perfect sphere -- the companion to haversine_distance() for callers that
 * need to walk along that path, not just measure it.
 */
double initial_bearing(double lat1, double lon1, double lat2, double lon2);

/*
 * Direct geodesic problem on a sphere: given a start point, an initial
 * bearing in degrees (as returned by initial_bearing()), and a distance in
 * the same units as `radius` (e.g. EARTH_RADIUS_KM for kilometers), writes
 * the (lat, lon) reached to outLat / outLon. Spherical, not ellipsoidal --
 * consistent with haversine_distance()'s assumption, and exact enough for
 * sampling points along a path (e.g. for an elevation profile).
 */
void destination_point(double lat1, double lon1, double bearingDeg, double distance, double radius,
                        double *outLat, double *outLon);

#ifdef __cplusplus
}
#endif

#endif /* GEO_UTILS_H */
