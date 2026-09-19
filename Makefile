.PHONY: clean

CC = gcc
CFLAGS = -std=c11 -g -Wall -Wno-unused -fPIC
LDFLAGS = -lm

# Machine-specific path to a locally-built raylib checkout (needed only for
# `terrain_viewer`, not `dtm_test`) -- see README.md.
RAYLIB_DIR = /home/murf/external-sources/raylib/src

geo_utils.o: geo_utils.c include/geo_utils.h
	$(CC) $(CFLAGS) -Iinclude -c geo_utils.c -o $@

dtm.o: dtm.c include/dtm.h
	$(CC) $(CFLAGS) -Iinclude -c dtm.c -o $@

dtm_fetch.o: dtm_fetch.c include/dtm_fetch.h include/dtm.h
	$(CC) $(CFLAGS) -std=gnu11 -Iinclude -c dtm_fetch.c -o $@

# Standalone sanity check for dtm.c/dtm.h, with no raylib dependency --
# point it at one or more DTM GeoTIFFs and it prints elevation at a few
# known lat/lon probes plus a coarse per-file min/max scan. See README.md.
dtm_test: dtm.c include/dtm.h
	$(CC) $(CFLAGS) -std=gnu11 -O2 -DDTM_TEST -Iinclude dtm.c -o $@ -ltiff -lm

# -O3 matters a lot here: at -O0, raymath.h's Vector3* helpers don't get
# inlined despite being declared `static inline`, costing real per-frame CPU
# in the terrain/profile render loops.
terrain_viewer: terrain_viewer.c include/dtm.h include/dtm_fetch.h include/geo_utils.h geo_utils.o dtm.o dtm_fetch.o
	$(CC) $(CFLAGS) -std=gnu11 -O3 terrain_viewer.c geo_utils.o dtm.o dtm_fetch.o -o $@ -Iinclude -I$(RAYLIB_DIR) -L$(RAYLIB_DIR) -Wl,-rpath,$(RAYLIB_DIR) -lraylib -lGL -lm -lpthread -ldl -lrt -lX11 -ltiff

# Headless terrain-scouting helpers, no raylib dependency -- see README.md's
# SCOUTING TOOLS section. Both read points from stdin so they're easy to
# drive from a generated grid (e.g. a small Python destination_point loop)
# without opening a window per point.
probe_points: probe_points.c include/dtm.h dtm.o
	$(CC) $(CFLAGS) -std=gnu11 -O2 -Iinclude probe_points.c dtm.o -o $@ -ltiff -lm

los_check: los_check.c include/dtm.h include/geo_utils.h dtm.o geo_utils.o
	$(CC) $(CFLAGS) -std=gnu11 -O2 -Iinclude los_check.c dtm.o geo_utils.o -o $@ -ltiff -lm

# Standalone smoke test for dtm_fetch.c -- see README.md's DATA SOURCES
# section. Usage: ./dtm_fetch_test cacheDir lat lon
dtm_fetch_test: dtm_fetch_test.c include/dtm_fetch.h include/dtm.h dtm.o dtm_fetch.o
	$(CC) $(CFLAGS) -std=gnu11 -O2 -Iinclude dtm_fetch_test.c dtm.o dtm_fetch.o -o $@ -ltiff -lm

clean:
	rm -f geo_utils.o dtm.o dtm_fetch.o dtm_test terrain_viewer probe_points los_check dtm_fetch_test
