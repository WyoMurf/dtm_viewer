/* Scratch tool: read "lat,lon" pairs from stdin, print "lat,lon,elevM".
 * Not part of the regular build -- ad hoc terrain-scouting helper. */
#include <stdio.h>
#include "dtm.h"

int main(int argc, char **argv) {
    DtmSet set = { 0 };
    for (int i = 1; i < argc; i++) DtmSetAddFile(&set, argv[i]);

    char line[256];
    while (fgets(line, sizeof(line), stdin)) {
        double lat, lon;
        if (sscanf(line, "%lf,%lf", &lat, &lon) != 2) continue;
        double m = DtmSampleMeters(&set, lat, lon, 0);
        printf("%.6f,%.6f,%.1f\n", lat, lon, m);
    }
    return 0;
}
