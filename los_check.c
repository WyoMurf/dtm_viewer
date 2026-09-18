/* Scratch tool: batch line-of-sight check, reusing the exact curvature/
 * refraction sightline math RunProfileMode uses, but headless (no raylib,
 * no window) so many pairs can be tested per second. Not part of the
 * regular build -- ad hoc terrain-scouting helper.
 *
 * Reads lines "lat1,lon1,h1ft,lat2,lon2,h2ft" from stdin (h in feet above
 * ground at each point), prints "lat1,lon1,lat2,lon2,CLEAR" or
 * "...,BLOCKED,worstMarginM,atKm" -- worstMarginM is how far the worst
 * obstruction pokes ABOVE the direct line (negative would mean the whole
 * path clears with that much room to spare, but this tool only reports it
 * for genuinely blocked paths). */
#include <stdio.h>
#include <math.h>
#include "dtm.h"
#include "geo_utils.h"

#define FEET_TO_METERS 0.3048
#define SAMPLES 400

static DtmSet g_dtm = { 0 };

int main(int argc, char **argv) {
    for (int i = 1; i < argc; i++) DtmSetAddFile(&g_dtm, argv[i]);

    char line[256];
    while (fgets(line, sizeof(line), stdin)) {
        double lat1, lon1, h1, lat2, lon2, h2;
        if (sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf", &lat1, &lon1, &h1, &lat2, &lon2, &h2) != 6) continue;

        double g1 = DtmSampleMeters(&g_dtm, lat1, lon1, 0);
        double g2 = DtmSampleMeters(&g_dtm, lat2, lon2, 0);
        if (g1 == DTM_NODATA || g2 == DTM_NODATA) {
            printf("%.6f,%.6f,%.6f,%.6f,NODATA\n", lat1, lon1, lat2, lon2);
            continue;
        }
        double elev1 = g1 + h1 * FEET_TO_METERS, elev2 = g2 + h2 * FEET_TO_METERS;

        double totalKm = haversine_distance(lat1, lon1, lat2, lon2, EARTH_RADIUS_KM);
        double bearing = initial_bearing(lat1, lon1, lat2, lon2);
        double effRadiusKm = EARTH_RADIUS_KM * (4.0 / 3.0);

        double worstMargin = -1e18, worstKm = 0;
        for (int i = 0; i <= SAMPLES; i++) {
            double d = totalKm * i / SAMPLES;
            double lat, lon;
            if (i == 0) { lat = lat1; lon = lon1; }
            else if (i == SAMPLES) { lat = lat2; lon = lon2; }
            else destination_point(lat1, lon1, bearing, d, EARTH_RADIUS_KM, &lat, &lon);
            double m = DtmSampleMeters(&g_dtm, lat, lon, d < 2.0 ? 0 : d < 8.0 ? 2 : 4);
            if (m == DTM_NODATA) continue;
            double chordM = elev1 + (elev2 - elev1) * (totalKm > 0 ? d / totalKm : 0);
            double bulgeM = 1000.0 * (d * (totalKm - d)) / (2.0 * effRadiusKm);
            double margin = m - (chordM + bulgeM);
            if (margin > worstMargin) { worstMargin = margin; worstKm = d; }
        }

        if (worstMargin > 0) {
            printf("%.6f,%.6f,%.6f,%.6f,BLOCKED,%.1f,%.2f\n", lat1, lon1, lat2, lon2, worstMargin, worstKm);
        } else {
            printf("%.6f,%.6f,%.6f,%.6f,CLEAR,%.1f\n", lat1, lon1, lat2, lon2, -worstMargin);
        }
    }
    return 0;
}
