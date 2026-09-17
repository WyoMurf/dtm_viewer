#include <math.h>

#include "geo_utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define TO_RAD (M_PI / 180.0)

double haversine_distance(double lat1, double lon1, double lat2, double lon2, double radius) {
    double phi1 = lat1 * TO_RAD;
    double phi2 = lat2 * TO_RAD;
    double dphi = (lat2 - lat1) * TO_RAD;
    double dlambda = (lon2 - lon1) * TO_RAD;

    double a = sin(dphi / 2.0) * sin(dphi / 2.0) +
                cos(phi1) * cos(phi2) * sin(dlambda / 2.0) * sin(dlambda / 2.0);

    double c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));

    return radius * c;
}

double initial_bearing(double lat1, double lon1, double lat2, double lon2) {
    double phi1 = lat1 * TO_RAD;
    double phi2 = lat2 * TO_RAD;
    double dlambda = (lon2 - lon1) * TO_RAD;

    double y = sin(dlambda) * cos(phi2);
    double x = cos(phi1) * sin(phi2) - sin(phi1) * cos(phi2) * cos(dlambda);
    double theta = atan2(y, x);
    double deg = theta / TO_RAD;
    return fmod(deg + 360.0, 360.0);
}

void destination_point(double lat1, double lon1, double bearingDeg, double distance, double radius,
                        double *outLat, double *outLon) {
    double phi1 = lat1 * TO_RAD;
    double lambda1 = lon1 * TO_RAD;
    double theta = bearingDeg * TO_RAD;
    double delta = distance / radius; /* angular distance */

    double phi2 = asin(sin(phi1) * cos(delta) + cos(phi1) * sin(delta) * cos(theta));
    double lambda2 = lambda1 + atan2(sin(theta) * sin(delta) * cos(phi1),
                                      cos(delta) - sin(phi1) * sin(phi2));

    *outLat = phi2 / TO_RAD;
    *outLon = lambda2 / TO_RAD;
}
