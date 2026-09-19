/* Standalone smoke test for dtm_fetch.c: ./dtm_fetch_test cacheDir lat lon
 * Ensures coverage at (lat,lon), then prints the sampled elevation there.
 * TV_AUTO_DOWNLOAD=1 skips the interactive confirmation prompt. */
#include <stdio.h>
#include <stdlib.h>
#include "dtm_fetch.h"

int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "usage: %s cacheDir lat lon\n", argv[0]);
        return 1;
    }
    const char *cacheDir = argv[1];
    double lat = atof(argv[2]), lon = atof(argv[3]);

    DtmSet set = { 0 };
    int rc = EnsureDtmCoverage(&set, lat, lon, cacheDir, getenv("TV_AUTO_DOWNLOAD") != NULL);
    if (rc != 0) {
        fprintf(stderr, "EnsureDtmCoverage failed: rc=%d\n", rc);
        return 1;
    }
    double m = DtmSampleMeters(&set, lat, lon, 0);
    printf("(%.5f, %.5f) -> %.1f m\n", lat, lon, m);
    DtmSetClose(&set);
    return 0;
}
