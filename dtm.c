#define _POSIX_C_SOURCE 200809L /* strdup */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include <tiffio.h>

#include "dtm.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * GeoTIFF's key tags (ModelPixelScale 33550, ModelTiepoint 33922, GDAL's
 * private NODATA 42113) aren't known to plain libtiff -- that's normally
 * libgeotiff's job, which isn't installed here. libtiff *does* have a
 * TIFFSetTagExtender/TIFFFieldInfo mechanism meant for exactly this, but as
 * of libtiff 4.6 it can't actually represent a FIELD_CUSTOM variable-count
 * DOUBLE array (TIFFFetchNormalTag logs "set_field_type ... is
 * TIFF_SETGET_UNDEFINED and thus tag is not read from file" and silently
 * drops it) -- confirmed by testing against these exact files. Rather than
 * reach into libtiff's non-public TIFFField internals to work around that,
 * these three tags are read directly with a small raw IFD-0 scan (plain
 * fopen/fread, no libtiff involved) -- everything else (tile layout,
 * pyramid levels, and the actual LZW-compressed pixel data) still goes
 * through libtiff via TIFF* below, which handles those tags natively.
 */
#define DTM_TAG_GEOPIXELSCALE   33550
#define DTM_TAG_GEOTIEPOINTS    33922
#define DTM_TAG_GDAL_NODATA     42113
#define DTM_TAG_GEOKEYDIRECTORY 34735

/* GeoKey IDs read out of GeoKeyDirectoryTag -- see ReadGeoTagsRaw. Only
 * the two needed to tell "geographic lat/lon" (Wyoming, USGS 1/3 arc-
 * second) apart from "projected UTM" (USGS's 1-meter DEM -- confirmed
 * directly by parsing a real downloaded tile's GeoKeyDirectoryTag: KeyID
 * 1024=1 (ModelTypeProjected) and KeyID 3072=26911 (EPSG:26911, NAD83 /
 * UTM zone 11N) for a tile covering 38.5N,116.5W in Nevada). */
#define GEOKEY_GTMODELTYPE    1024 /* 1 = Projected, 2 = Geographic */
#define GEOKEY_PROJECTEDCSTYPE 3072 /* EPSG code when GTModelType == 1 */

typedef struct {
    FILE *fp;
    int big;           /* BigTIFF (offsets are 8 bytes) vs classic (4 bytes) */
    int little_endian;
} RawTiff;

static uint64_t RawReadUint(RawTiff *rt, int nbytes) {
    unsigned char b[8] = { 0 };
    if (fread(b, 1, (size_t)nbytes, rt->fp) != (size_t)nbytes) return 0;
    uint64_t v = 0;
    if (rt->little_endian) {
        for (int i = nbytes - 1; i >= 0; i--) v = (v << 8) | b[i];
    } else {
        for (int i = 0; i < nbytes; i++) v = (v << 8) | b[i];
    }
    return v;
}

static double RawReadDouble(RawTiff *rt) {
    unsigned char b[8];
    if (fread(b, 1, 8, rt->fp) != 8) return 0.0;
    if (!rt->little_endian) {
        for (int i = 0; i < 4; i++) { unsigned char t = b[i]; b[i] = b[7 - i]; b[7 - i] = t; }
    }
    double d;
    memcpy(&d, b, 8);
    return d;
}

/* Reads tags 33550/33922/42113/34735 out of IFD 0 (the geo tags are the same
 * across every overview level, so there's no need to look past directory 0).
 * outScale/outTiepoint must each have room for 6 doubles; *outScaleCount/
 * *outTiepointCount report how many were actually present. outNodata is a
 * caller-supplied buffer of outNodataSize bytes for the ASCII value (left
 * untouched, with *outHasNodata set to 0, if the tag isn't present).
 * *outIsProjected is set to 1 if GeoKeyDirectoryTag's GTModelTypeGeoKey
 * (1024) says Projected rather than Geographic, and *outEpsg to whatever
 * ProjectedCSTypeGeoKey (3072) gave in that case (0 otherwise) -- see
 * DtmSetAddFile's UTM handling. */
static int ReadGeoTagsRaw(const char *path, double *outScale, int *outScaleCount,
                           double *outTiepoint, int *outTiepointCount,
                           char *outNodata, size_t outNodataSize, int *outHasNodata,
                           int *outIsProjected, int *outEpsg) {
    *outScaleCount = 0;
    *outTiepointCount = 0;
    *outHasNodata = 0;
    *outIsProjected = 0;
    *outEpsg = 0;

    FILE *fp = fopen(path, "rb");
    if (!fp) return -1;

    unsigned char hdr[16];
    if (fread(hdr, 1, 16, fp) != 16) { fclose(fp); return -1; }

    RawTiff rt;
    rt.fp = fp;
    rt.little_endian = (hdr[0] == 'I' && hdr[1] == 'I');
    uint16_t magic = rt.little_endian ? (uint16_t)(hdr[2] | (hdr[3] << 8)) : (uint16_t)(hdr[3] | (hdr[2] << 8));
    rt.big = (magic == 43);

    uint64_t first_ifd;
    if (rt.big) {
        fseek(fp, 8, SEEK_SET);
        first_ifd = RawReadUint(&rt, 8);
    } else {
        fseek(fp, 4, SEEK_SET);
        first_ifd = RawReadUint(&rt, 4);
    }

    int count_field_bytes = rt.big ? 8 : 4;
    int value_field_bytes = rt.big ? 8 : 4;
    int entry_size = 2 + 2 + count_field_bytes + value_field_bytes;

    fseek(fp, (long)first_ifd, SEEK_SET);
    uint64_t n_entries = RawReadUint(&rt, rt.big ? 8 : 2);
    long entries_start = (long)first_ifd + (rt.big ? 8 : 2);

    for (uint64_t e = 0; e < n_entries; e++) {
        long entry_pos = entries_start + (long)e * entry_size;
        fseek(fp, entry_pos, SEEK_SET);
        uint16_t tag = (uint16_t)RawReadUint(&rt, 2);
        uint16_t type = (uint16_t)RawReadUint(&rt, 2);
        uint64_t count = RawReadUint(&rt, count_field_bytes);

        if (tag != DTM_TAG_GEOPIXELSCALE && tag != DTM_TAG_GEOTIEPOINTS && tag != DTM_TAG_GDAL_NODATA && tag != DTM_TAG_GEOKEYDIRECTORY) continue;

        int type_size = (type == 2) ? 1 : (type == 12) ? 8 : (type == 3) ? 2 : 0; /* ASCII, DOUBLE, SHORT (GeoKeyDirectory) */
        if (type_size == 0) continue;
        uint64_t total_bytes = (uint64_t)type_size * count;

        long value_field_pos = entry_pos + 4 + count_field_bytes;
        uint64_t data_offset;
        if (total_bytes <= (uint64_t)value_field_bytes) {
            data_offset = (uint64_t)value_field_pos; /* value is inline, right here */
        } else {
            fseek(fp, value_field_pos, SEEK_SET);
            data_offset = RawReadUint(&rt, value_field_bytes);
        }

        if (tag == DTM_TAG_GEOPIXELSCALE || tag == DTM_TAG_GEOTIEPOINTS) {
            int n = (int)count;
            if (n > 6) n = 6;
            double *out = (tag == DTM_TAG_GEOPIXELSCALE) ? outScale : outTiepoint;
            fseek(fp, (long)data_offset, SEEK_SET);
            for (int i = 0; i < n; i++) out[i] = RawReadDouble(&rt);
            if (tag == DTM_TAG_GEOPIXELSCALE) *outScaleCount = n;
            else *outTiepointCount = n;
        } else if (tag == DTM_TAG_GEOKEYDIRECTORY) {
            /* Array of SHORTs: a 4-value header (version, keyRevision,
             * minorRevision, numberOfKeys) then numberOfKeys groups of 4
             * (KeyID, TIFFTagLocation, Count, Value/Offset). Only the two
             * keys that matter here (1024, 3072) are ever inline SHORT
             * values (TIFFTagLocation 0) in practice, so that's all this
             * reads -- not a general GeoKey parser. */
            fseek(fp, (long)data_offset, SEEK_SET);
            uint16_t numKeys = 0;
            if (count >= 4) {
                RawReadUint(&rt, 2); /* KeyDirectoryVersion, always 1 */
                RawReadUint(&rt, 2); /* KeyRevision */
                RawReadUint(&rt, 2); /* MinorRevision */
                numKeys = (uint16_t)RawReadUint(&rt, 2);
            }
            for (uint16_t k = 0; k < numKeys && (uint64_t)(4 + (k + 1) * 4) <= count; k++) {
                uint16_t keyId = (uint16_t)RawReadUint(&rt, 2);
                uint16_t loc = (uint16_t)RawReadUint(&rt, 2);
                RawReadUint(&rt, 2); /* keyCount, unused for inline SHORT keys */
                uint16_t val = (uint16_t)RawReadUint(&rt, 2);
                if (loc != 0) continue; /* not an inline value -- not needed for 1024/3072 */
                if (keyId == GEOKEY_GTMODELTYPE) *outIsProjected = (val == 1);
                else if (keyId == GEOKEY_PROJECTEDCSTYPE) *outEpsg = val;
            }
        } else { /* DTM_TAG_GDAL_NODATA, ASCII */
            size_t n = (size_t)count;
            if (n >= outNodataSize) n = outNodataSize - 1;
            fseek(fp, (long)data_offset, SEEK_SET);
            if (fread(outNodata, 1, n, fp) == n) {
                outNodata[n] = '\0';
                *outHasNodata = 1;
            }
        }
    }

    fclose(fp);
    return 0;
}

/* Decodes a UTM EPSG code into zone/hemisphere. NAD83 zones (26901-26923,
 * used by USGS's 1-meter DEM -- CONUS never needs zone > 23) and WGS84
 * zones (32601-32660 north, 32701-32760 south) both map to the same
 * ellipsoid for this project's purposes (see the datum comment on
 * EnsureDtmCoverage -- NAD83 vs WGS84 is a ~1-2m discrepancy, already
 * accepted as negligible elsewhere in this project). Returns 1 and fills
 * outZone (1-60)/outNorth on a recognized UTM code, 0 otherwise. */
static int EpsgToUtmZone(int epsg, int *outZone, int *outNorth) {
    if (epsg >= 26901 && epsg <= 26923) { *outZone = epsg - 26900; *outNorth = 1; return 1; }
    if (epsg >= 32601 && epsg <= 32660) { *outZone = epsg - 32600; *outNorth = 1; return 1; }
    if (epsg >= 32701 && epsg <= 32760) { *outZone = epsg - 32700; *outNorth = 0; return 1; }
    return 0;
}

/* WGS84/GRS80 ellipsoid -- functionally identical for NAD83 too (both
 * derive from the same reference ellipsoid to well under 1mm difference,
 * far tighter than the ~1-2m datum-origin discrepancy already accepted
 * elsewhere). Standard non-iterative Snyder transverse Mercator forward/
 * inverse formulas (used by PROJ/GDAL and every other UTM implementation
 * this project's numbers were cross-checked against), accurate to
 * sub-meter within a single UTM zone -- far tighter than this tool needs
 * at 1m pixel resolution. */
#define UTM_A  6378137.0
#define UTM_F  (1.0 / 298.257223563)
#define UTM_K0 0.9996
#define UTM_E0 500000.0

static void LatLonToUtm(double lat, double lon, int zone, int north, double *outE, double *outN) {
    double e2 = UTM_F * (2.0 - UTM_F);
    double ep2 = e2 / (1.0 - e2);
    double phi = lat * M_PI / 180.0;
    double lambda = lon * M_PI / 180.0;
    double lambda0 = ((zone - 1) * 6 - 180 + 3) * M_PI / 180.0;

    double sinPhi = sin(phi), cosPhi = cos(phi), tanPhi = tan(phi);
    double Nrad = UTM_A / sqrt(1.0 - e2 * sinPhi * sinPhi);
    double T = tanPhi * tanPhi;
    double C = ep2 * cosPhi * cosPhi;
    double A = cosPhi * (lambda - lambda0);
    double M = UTM_A * ((1 - e2/4 - 3*e2*e2/64 - 5*e2*e2*e2/256) * phi
                         - (3*e2/8 + 3*e2*e2/32 + 45*e2*e2*e2/1024) * sin(2*phi)
                         + (15*e2*e2/256 + 45*e2*e2*e2/1024) * sin(4*phi)
                         - (35*e2*e2*e2/3072) * sin(6*phi));

    *outE = UTM_E0 + UTM_K0 * Nrad * (A + (1-T+C)*A*A*A/6.0 + (5-18*T+T*T+72*C-58*ep2)*A*A*A*A*A/120.0);
    *outN = UTM_K0 * (M + Nrad*tanPhi*(A*A/2.0 + (5-T+9*C+4*C*C)*A*A*A*A/24.0 + (61-58*T+T*T+600*C-330*ep2)*A*A*A*A*A*A/720.0));
    if (!north) *outN += 10000000.0;
}

static void UtmToLatLon(double e, double n, int zone, int north, double *outLat, double *outLon) {
    double e2 = UTM_F * (2.0 - UTM_F);
    double ep2 = e2 / (1.0 - e2);
    double e1 = (1.0 - sqrt(1.0 - e2)) / (1.0 + sqrt(1.0 - e2));
    double lambda0 = ((zone - 1) * 6 - 180 + 3) * M_PI / 180.0;

    double x = e - UTM_E0;
    double y = north ? n : n - 10000000.0;
    double M = y / UTM_K0;
    double mu = M / (UTM_A * (1 - e2/4 - 3*e2*e2/64 - 5*e2*e2*e2/256));

    double phi1 = mu + (3*e1/2 - 27*e1*e1*e1/32)*sin(2*mu)
                      + (21*e1*e1/16 - 55*e1*e1*e1*e1/32)*sin(4*mu)
                      + (151*e1*e1*e1/96)*sin(6*mu)
                      + (1097*e1*e1*e1*e1/512)*sin(8*mu);

    double sinPhi1 = sin(phi1), cosPhi1 = cos(phi1), tanPhi1 = tan(phi1);
    double C1 = ep2 * cosPhi1 * cosPhi1;
    double T1 = tanPhi1 * tanPhi1;
    double N1 = UTM_A / sqrt(1.0 - e2 * sinPhi1 * sinPhi1);
    double R1 = UTM_A * (1.0 - e2) / pow(1.0 - e2 * sinPhi1 * sinPhi1, 1.5);
    double D = x / (N1 * UTM_K0);

    double phi = phi1 - (N1 * tanPhi1 / R1) * (D*D/2.0
                 - (5+3*T1+10*C1-4*C1*C1-9*ep2)*D*D*D*D/24.0
                 + (61+90*T1+298*C1+45*T1*T1-252*ep2-3*C1*C1)*D*D*D*D*D*D/720.0);
    double lambda = lambda0 + (D - (1+2*T1+C1)*D*D*D/6.0
                    + (5-2*C1+28*T1-3*C1*C1+8*ep2+24*T1*T1)*D*D*D*D*D/120.0) / cosPhi1;

    *outLat = phi * 180.0 / M_PI;
    *outLon = lambda * 180.0 / M_PI;
}

#define DTM_MAX_LEVELS 16
#define DTM_TILE_CACHE_SLOTS 2048

typedef struct {
    int valid;
    int level;
    uint32_t tile_col, tile_row;
    unsigned char *data; /* level_info[level].tile_w * tile_h samples, bytes_per_sample bytes each -- raw bytes, not a typed pointer, since a DtmSet can mix uint16 (Wyoming) and float32 (USGS) files */
    size_t capacity;     /* bytes currently allocated in `data` */
} TileCacheSlot;

typedef struct {
    int queried;
    uint32_t width, height, tile_w, tile_h;
} LevelInfo;

struct DtmFile {
    TIFF *tiff;
    char *path;
    int current_dir;
    int num_levels;
    /* Native coordinate bounds -- level-independent, every overview covers
     * the same extent. Degrees (lon/lat) for a geographic file, meters
     * (easting/northing) for a projected UTM file -- see isUtm below and
     * FileNativeXY, which converts an incoming query point into whichever
     * system this particular file actually uses before anything touches
     * these bounds. */
    double xMin, xMax, yMin, yMax;
    int isUtm;        /* 0 = geographic lat/lon (Wyoming, USGS 1/3as), 1 = projected UTM (USGS 1m) */
    int utmZone;       /* 1-60, meaningless unless isUtm */
    int utmNorth;      /* 1 = northern hemisphere, meaningless unless isUtm */
    double nodata;
    int sample_is_float;  /* 0 = uint16 (Wyoming wyolidar), 1 = float32 (USGS 3DEP) -- see DtmSetAddFile */
    int bytes_per_sample; /* 2 or 4, matching sample_is_float */
    LevelInfo level_info[DTM_MAX_LEVELS];
    TileCacheSlot cache[DTM_TILE_CACHE_SLOTS];
};

static int EnsureLevelInfo(DtmFile *f, int level) {
    LevelInfo *li = &f->level_info[level];
    if (li->queried) return 1;
    if (f->current_dir != level) {
        if (!TIFFSetDirectory(f->tiff, (tdir_t)level)) {
            fprintf(stderr, "dtm: %s: couldn't seek to directory %d\n", f->path, level);
            return 0;
        }
        f->current_dir = level;
    }
    uint32_t w = 0, h = 0, tw = 0, th = 0;
    TIFFGetField(f->tiff, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(f->tiff, TIFFTAG_IMAGELENGTH, &h);
    if (!TIFFGetField(f->tiff, TIFFTAG_TILEWIDTH, &tw) || !TIFFGetField(f->tiff, TIFFTAG_TILELENGTH, &th)) {
        fprintf(stderr, "dtm: %s: directory %d isn't tiled\n", f->path, level);
        return 0;
    }
    li->width = w;
    li->height = h;
    li->tile_w = tw;
    li->tile_h = th;
    li->queried = 1;
    return 1;
}

static void SilenceLibtiffWarningsOnce(void) {
    /* libtiff logs "Unknown field with tag ..." for every custom GeoTIFF/GDAL
     * tag (33550, 33922, 34735, 42112, 42113, ...) every time it parses a
     * directory -- expected and harmless (ReadGeoTagsRaw above reads the
     * ones we need directly), but noisy enough across 11 pyramid levels
     * times 3 files that it drowns out real warnings. */
    static int done = 0;
    if (done) return;
    done = 1;
    TIFFSetWarningHandler(NULL);
}

int DtmSetAddFile(DtmSet *set, const char *path) {
    SilenceLibtiffWarningsOnce();
    if (set->file_count >= DTM_MAX_FILES) {
        fprintf(stderr, "dtm: %s: DtmSet already has %d files (DTM_MAX_FILES)\n", path, DTM_MAX_FILES);
        return -1;
    }

    double scale[6], tiepoint[6];
    int scaleCount = 0, tiepointCount = 0, hasNodata = 0, isProjected = 0, epsg = 0;
    char nodataStr[64];
    if (ReadGeoTagsRaw(path, scale, &scaleCount, tiepoint, &tiepointCount, nodataStr, sizeof(nodataStr), &hasNodata, &isProjected, &epsg) != 0) {
        fprintf(stderr, "dtm: couldn't open %s to read its GeoTIFF tags\n", path);
        return -1;
    }
    if (scaleCount < 2 || tiepointCount < 6) {
        fprintf(stderr, "dtm: %s: missing ModelPixelScale/ModelTiepoint GeoTIFF tags -- not a recognized DTM layout\n", path);
        return -1;
    }
    int utmZone = 0, utmNorth = 0;
    if (isProjected && !EpsgToUtmZone(epsg, &utmZone, &utmNorth)) {
        fprintf(stderr, "dtm: %s: projected CRS EPSG:%d isn't a recognized UTM zone -- not a recognized DTM layout\n", path, epsg);
        return -1;
    }

    /* No explicit mmap call here on purpose: this "r" mode (as opposed to
     * "rm") already gets one from libtiff itself. On Unix, TIFFOpen's
     * default client mmaps the ENTIRE file read-only the moment it's
     * opened, and TIFFReadEncodedTile (in GetTile below) reads compressed
     * tile bytes straight out of that mapping rather than calling
     * read()/pread() -- confirmed with strace against a live 4.4GB DTM
     * file (one mmap() of the whole file at TIFFOpen time, then zero
     * further read/pread syscalls on that fd for the rest of the run).
     * That mapping is lazy at the OS page-cache level exactly the way
     * kdtree's own mmap'd shard files are: opening reserves address
     * space, but pages only fault in from disk as tiles are actually
     * touched. Do NOT "fix" this by adding a second, redundant explicit
     * mmap layer -- there's nothing for it to improve on. */
    TIFF *tiff = TIFFOpen(path, "r");
    if (!tiff) {
        fprintf(stderr, "dtm: couldn't open %s\n", path);
        return -1;
    }

    /* Two known layouts: Wyoming's wyolidar COGs (16-bit unsigned int
     * meters) and USGS 3DEP's national tiles (32-bit float meters,
     * confirmed against a real downloaded n45w109 1/3-arc-second tile --
     * SampleFormat=3/IEEEFP, BitsPerSample=32). Both single-band. */
    uint16_t bits = 0, sampleFormat = 0, samplesPerPixel = 0;
    TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &bits);
    TIFFGetFieldDefaulted(tiff, TIFFTAG_SAMPLEFORMAT, &sampleFormat);
    TIFFGetFieldDefaulted(tiff, TIFFTAG_SAMPLESPERPIXEL, &samplesPerPixel);
    int isUint16 = (bits == 16 && sampleFormat == SAMPLEFORMAT_UINT);
    int isFloat32 = (bits == 32 && sampleFormat == SAMPLEFORMAT_IEEEFP);
    if (!(isUint16 || isFloat32) || samplesPerPixel != 1) {
        fprintf(stderr, "dtm: %s: expected single-band 16-bit unsigned int or 32-bit float (got %u bits, sampleFormat=%u, samples/pixel=%u) -- not a recognized DTM layout\n",
                path, bits, sampleFormat, samplesPerPixel);
        TIFFClose(tiff);
        return -1;
    }

    uint32_t width = 0, height = 0;
    TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &height);

    /* Tiepoint is (I,J,K, X,Y,Z) for one raster corner -- these files always
     * tie (0,0) to the northwest corner, but compute from whatever I,J is
     * given rather than assuming, in case that ever changes. X/Y here are
     * in whatever coordinate system the file itself uses (degrees for a
     * geographic file, UTM meters for a projected one) -- the math is
     * identical either way, only the units differ. */
    double i0 = tiepoint[0], j0 = tiepoint[1], x0 = tiepoint[3], y0 = tiepoint[4];
    double sx = scale[0], sy = scale[1];
    double xNW = x0 - i0 * sx;
    double yNW = y0 + j0 * sy;
    double xSE = xNW + width * sx;
    double ySE = yNW - height * sy;

    DtmFile *f = calloc(1, sizeof(DtmFile));
    if (!f) { fprintf(stderr, "dtm: out of memory\n"); TIFFClose(tiff); return -1; }
    f->tiff = tiff;
    f->path = strdup(path);
    f->current_dir = 0;
    f->xMin = xNW < xSE ? xNW : xSE;
    f->xMax = xNW < xSE ? xSE : xNW;
    f->yMin = ySE < yNW ? ySE : yNW;
    f->yMax = ySE < yNW ? yNW : ySE;
    f->isUtm = isProjected;
    f->utmZone = utmZone;
    f->utmNorth = utmNorth;
    f->sample_is_float = isFloat32;
    f->bytes_per_sample = isFloat32 ? 4 : 2;

    /* Fallback defaults if a file is somehow missing its own NODATA tag --
     * 65535 for Wyoming's uint16 tiles (observed default there), -999999
     * for float32 (USGS 3DEP's own convention, confirmed on a real tile). */
    f->nodata = hasNodata ? atof(nodataStr) : (isFloat32 ? -999999.0 : 65535.0);

    int ndirs = TIFFNumberOfDirectories(tiff);
    if (ndirs > DTM_MAX_LEVELS) ndirs = DTM_MAX_LEVELS;
    f->num_levels = ndirs;

    if (!EnsureLevelInfo(f, 0)) {
        fprintf(stderr, "dtm: %s: couldn't read directory 0's tile layout\n", path);
        free(f->path);
        free(f);
        TIFFClose(tiff);
        return -1;
    }

    set->files[set->file_count++] = f;
    if (f->isUtm) {
        double lat1, lon1, lat2, lon2;
        UtmToLatLon(f->xMin, f->yMin, f->utmZone, f->utmNorth, &lat1, &lon1);
        UtmToLatLon(f->xMax, f->yMax, f->utmZone, f->utmNorth, &lat2, &lon2);
        printf("dtm: opened %s: %dx%d px, %d pyramid level%s, UTM zone %d%c E[%.0f, %.0f] N[%.0f, %.0f] (~lon [%.4f, %.4f] lat [%.4f, %.4f]), nodata=%.0f\n",
               path, (int)f->level_info[0].width, (int)f->level_info[0].height, f->num_levels,
               f->num_levels == 1 ? "" : "s", f->utmZone, f->utmNorth ? 'N' : 'S',
               f->xMin, f->xMax, f->yMin, f->yMax, lon1, lon2, lat1, lat2, f->nodata);
    } else {
        printf("dtm: opened %s: %dx%d px, %d pyramid level%s, lon [%.4f, %.4f] lat [%.4f, %.4f], nodata=%.0f\n",
               path, (int)f->level_info[0].width, (int)f->level_info[0].height, f->num_levels,
               f->num_levels == 1 ? "" : "s", f->xMin, f->xMax, f->yMin, f->yMax, f->nodata);
    }
    return 0;
}

void DtmSetClose(DtmSet *set) {
    for (int i = 0; i < set->file_count; i++) {
        DtmFile *f = set->files[i];
        for (int s = 0; s < DTM_TILE_CACHE_SLOTS; s++) free(f->cache[s].data);
        TIFFClose(f->tiff);
        free(f->path);
        free(f);
    }
    set->file_count = 0;
}

/* Direct-mapped tile cache: a hash collision just means an eviction (an
 * extra re-decode next time), never incorrectness -- fine at the working
 * set a single walkthrough/profile query touches. */
static unsigned char *GetTile(DtmFile *f, int level, uint32_t tile_col, uint32_t tile_row) {
    LevelInfo *li = &f->level_info[level];
    size_t tile_samples = (size_t)li->tile_w * li->tile_h;
    size_t tile_bytes = tile_samples * (size_t)f->bytes_per_sample;

    uint64_t h = (uint64_t)level * 1000003ull + (uint64_t)tile_col * 2654435761ull + (uint64_t)tile_row * 40503ull;
    TileCacheSlot *slot = &f->cache[h % DTM_TILE_CACHE_SLOTS];

    if (slot->valid && slot->level == level && slot->tile_col == tile_col && slot->tile_row == tile_row) {
        return slot->data;
    }

    if (slot->capacity < tile_bytes) {
        unsigned char *bigger = realloc(slot->data, tile_bytes);
        if (!bigger) { fprintf(stderr, "dtm: out of memory growing a tile cache slot\n"); return NULL; }
        slot->data = bigger;
        slot->capacity = tile_bytes;
    }

    if (f->current_dir != level) {
        if (!TIFFSetDirectory(f->tiff, (tdir_t)level)) {
            fprintf(stderr, "dtm: %s: couldn't seek to directory %d\n", f->path, level);
            return NULL;
        }
        f->current_dir = level;
    }

    uint32_t x = tile_col * li->tile_w;
    uint32_t y = tile_row * li->tile_h;
    tmsize_t bytesExpected = (tmsize_t)tile_bytes;
    tmsize_t bytesRead = TIFFReadTile(f->tiff, slot->data, x, y, 0, 0);
    if (bytesRead != bytesExpected) {
        fprintf(stderr, "dtm: %s: TIFFReadTile(level=%d, x=%u, y=%u) read %ld bytes, expected %ld\n",
                f->path, level, x, y, (long)bytesRead, (long)bytesExpected);
        slot->valid = 0;
        return NULL;
    }

    slot->valid = 1;
    slot->level = level;
    slot->tile_col = tile_col;
    slot->tile_row = tile_row;
    return slot->data;
}

/* Raw sample at full pixel coordinates (col, row) within `level`, or
 * f->nodata's sentinel meaning if the tile can't be read at all. Coordinates
 * are clamped to the level's raster bounds by the caller (SampleFile). */
static double SamplePixel(DtmFile *f, int level, uint32_t col, uint32_t row) {
    LevelInfo *li = &f->level_info[level];
    uint32_t tile_col = col / li->tile_w;
    uint32_t tile_row = row / li->tile_h;
    unsigned char *tile = GetTile(f, level, tile_col, tile_row);
    if (!tile) return f->nodata;
    uint32_t localCol = col % li->tile_w;
    uint32_t localRow = row % li->tile_h;
    size_t sampleIdx = (size_t)localRow * li->tile_w + localCol;
    if (f->sample_is_float) {
        float raw;
        memcpy(&raw, tile + sampleIdx * 4, 4);
        return (double)raw;
    }
    /* uint16 (Wyoming): reinterpret the same 16 bits unsigned (a plain
     * int16_t read would turn values >= 32768, like the 65535 nodata
     * sentinel, negative). */
    uint16_t raw;
    memcpy(&raw, tile + sampleIdx * 2, 2);
    return (double)raw;
}

/* Converts a query point into whichever coordinate system `f`'s bounds
 * are actually expressed in -- identity for a geographic file, a forward
 * UTM projection for a projected one. Doing the conversion here (once per
 * query point) rather than converting the file's own bounds to lat/lon
 * once at load time is what keeps FileCovers/SampleFile's bounds check
 * and pixel math *exact*: a UTM tile's true shape is an axis-aligned box
 * in UTM meters, not in lat/lon (the grid doesn't follow meridians), so
 * checking containment against a lat/lon-approximated box would be the
 * one introducing error, not this. */
static void FileNativeXY(const DtmFile *f, double lat, double lon, double *outX, double *outY) {
    if (f->isUtm) LatLonToUtm(lat, lon, f->utmZone, f->utmNorth, outX, outY);
    else { *outX = lon; *outY = lat; }
}

static int FileCovers(const DtmFile *f, double lat, double lon) {
    double x, y;
    FileNativeXY(f, lat, lon, &x, &y);
    return y >= f->yMin && y <= f->yMax && x >= f->xMin && x <= f->xMax;
}

static double SampleFile(DtmFile *f, double lat, double lon, int level) {
    if (level < 0) level = 0;
    if (level >= f->num_levels) level = f->num_levels - 1;
    if (!EnsureLevelInfo(f, level)) return DTM_NODATA;
    LevelInfo *li = &f->level_info[level];

    double x, y;
    FileNativeXY(f, lat, lon, &x, &y);
    double colf = (x - f->xMin) / (f->xMax - f->xMin) * li->width;
    double rowf = (f->yMax - y) / (f->yMax - f->yMin) * li->height;

    /* Pixel-center convention: pixel index p's center is at colf = p + 0.5. */
    double pc = colf - 0.5;
    double pr = rowf - 0.5;
    long c0 = (long)floor(pc), r0 = (long)floor(pr);
    double tx = pc - (double)c0, ty = pr - (double)r0;
    long c1 = c0 + 1, r1 = r0 + 1;

    if (c0 < 0) c0 = 0;
    if (c0 > (long)li->width - 1) c0 = li->width - 1;
    if (c1 < 0) c1 = 0;
    if (c1 > (long)li->width - 1) c1 = li->width - 1;
    if (r0 < 0) r0 = 0;
    if (r0 > (long)li->height - 1) r0 = li->height - 1;
    if (r1 < 0) r1 = 0;
    if (r1 > (long)li->height - 1) r1 = li->height - 1;

    double v00 = SamplePixel(f, level, (uint32_t)c0, (uint32_t)r0);
    double v10 = SamplePixel(f, level, (uint32_t)c1, (uint32_t)r0);
    double v01 = SamplePixel(f, level, (uint32_t)c0, (uint32_t)r1);
    double v11 = SamplePixel(f, level, (uint32_t)c1, (uint32_t)r1);

    int n00 = (v00 == f->nodata), n10 = (v10 == f->nodata), n01 = (v01 == f->nodata), n11 = (v11 == f->nodata);
    if (n00 && n10 && n01 && n11) return DTM_NODATA;
    /* Substitute any individually-nodata corner (e.g. right at a coverage
     * edge) with the average of whichever corners ARE valid, rather than
     * poisoning the whole bilinear blend with a 65535 outlier. */
    if (n00 || n10 || n01 || n11) {
        double sum = 0; int cnt = 0;
        if (!n00) { sum += v00; cnt++; } if (!n10) { sum += v10; cnt++; }
        if (!n01) { sum += v01; cnt++; } if (!n11) { sum += v11; cnt++; }
        double avg = sum / cnt;
        if (n00) v00 = avg;
        if (n10) v10 = avg;
        if (n01) v01 = avg;
        if (n11) v11 = avg;
    }

    double top = v00 + (v10 - v00) * tx;
    double bot = v01 + (v11 - v01) * tx;
    return top + (bot - top) * ty;
}

double DtmSampleMeters(DtmSet *set, double lat, double lon, int level) {
    for (int i = 0; i < set->file_count; i++) {
        if (FileCovers(set->files[i], lat, lon)) {
            return SampleFile(set->files[i], lat, lon, level);
        }
    }
    return DTM_NODATA;
}

int DtmSetCovers(const DtmSet *set, double lat, double lon) {
    for (int i = 0; i < set->file_count; i++) {
        if (FileCovers(set->files[i], lat, lon)) return 1;
    }
    return 0;
}

void DtmLonLatToLocalMeters(double originLat, double originLon, double lat, double lon, double *outEastM, double *outNorthM) {
    double metersPerDegLat = 111320.0; /* close enough everywhere; the ~0.5% pole-to-equator variation doesn't matter at walkthrough scale */
    double metersPerDegLon = 111320.0 * cos(originLat * M_PI / 180.0);
    *outEastM = (lon - originLon) * metersPerDegLon;
    *outNorthM = (lat - originLat) * metersPerDegLat;
}

void DtmLocalMetersToLonLat(double originLat, double originLon, double eastM, double northM, double *outLat, double *outLon) {
    double metersPerDegLat = 111320.0;
    double metersPerDegLon = 111320.0 * cos(originLat * M_PI / 180.0);
    *outLon = originLon + eastM / metersPerDegLon;
    *outLat = originLat + northM / metersPerDegLat;
}

#ifdef DTM_TEST
int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s file.tif [file.tif ...]\n", argv[0]);
        return 1;
    }
    DtmSet set = { 0 };
    for (int i = 1; i < argc; i++) {
        if (DtmSetAddFile(&set, argv[i]) != 0) return 1;
    }

    struct { const char *name; double lat, lon; } probes[] = {
        { "Cody townsite",       44.5263, -109.0565 },
        { "Meeteetse townsite",  44.1467, -108.8710 },
        { "Cloud Peak (Bighorn)", 44.3853, -107.1683 },
        { "Out of coverage (way off)", 0.0, 0.0 },
    };
    for (size_t i = 0; i < sizeof(probes) / sizeof(probes[0]); i++) {
        double m0 = DtmSampleMeters(&set, probes[i].lat, probes[i].lon, 0);
        printf("%-28s (%.4f, %.4f): %.1f m", probes[i].name, probes[i].lat, probes[i].lon, m0);
        if (m0 == DTM_NODATA) printf(" (NODATA)");
        printf("\n");
    }

    /* Coarse scan of each file's own full-res coverage to compare against
     * the STATISTICS_MINIMUM/MAXIMUM already read out of GDAL_METADATA by
     * hand (Bighorn 1094-4013m, Cody 1231-3811m, Meeteetse 1105-2756m). A
     * coarse (~200x200) scan at a mid overview level won't reproduce those
     * exactly (it can miss the single extreme pixel), but should land in
     * the right neighborhood -- a sanity check, not an exact match. */
    for (int fi = 0; fi < set.file_count; fi++) {
        DtmFile *f = set.files[fi];
        double mn = 1e18, mx = -1e18;
        int level = f->num_levels > 4 ? 4 : f->num_levels - 1;
        for (int yi = 0; yi <= 40; yi++) {
            double y = f->yMin + (f->yMax - f->yMin) * yi / 40.0;
            for (int xi = 0; xi <= 40; xi++) {
                double x = f->xMin + (f->xMax - f->xMin) * xi / 40.0;
                /* Interpolate in the file's own native units (degrees or
                 * UTM meters), then convert to lat/lon -- SampleFile always
                 * takes lat/lon and re-projects forward itself, so a UTM
                 * file needs this round trip rather than treating x/y as
                 * lat/lon directly. */
                double lat, lon;
                if (f->isUtm) UtmToLatLon(x, y, f->utmZone, f->utmNorth, &lat, &lon);
                else { lat = y; lon = x; }
                double m = SampleFile(f, lat, lon, level);
                if (m == DTM_NODATA) continue;
                if (m < mn) mn = m;
                if (m > mx) mx = m;
            }
        }
        printf("%s: coarse scan at level %d -> min=%.0fm max=%.0fm\n", f->path, level, mn, mx);
    }

    DtmSetClose(&set);
    return 0;
}
#endif
