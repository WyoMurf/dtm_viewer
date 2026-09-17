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
#define DTM_TAG_GEOPIXELSCALE 33550
#define DTM_TAG_GEOTIEPOINTS  33922
#define DTM_TAG_GDAL_NODATA   42113

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

/* Reads just tags 33550/33922/42113 out of IFD 0 (the geo tags are the same
 * across every overview level, so there's no need to look past directory 0).
 * outScale/outTiepoint must each have room for 6 doubles; *outScaleCount/
 * *outTiepointCount report how many were actually present. outNodata is a
 * caller-supplied buffer of outNodataSize bytes for the ASCII value (left
 * untouched, with *outHasNodata set to 0, if the tag isn't present). */
static int ReadGeoTagsRaw(const char *path, double *outScale, int *outScaleCount,
                           double *outTiepoint, int *outTiepointCount,
                           char *outNodata, size_t outNodataSize, int *outHasNodata) {
    *outScaleCount = 0;
    *outTiepointCount = 0;
    *outHasNodata = 0;

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

        if (tag != DTM_TAG_GEOPIXELSCALE && tag != DTM_TAG_GEOTIEPOINTS && tag != DTM_TAG_GDAL_NODATA) continue;

        int type_size = (type == 2) ? 1 : (type == 12) ? 8 : 0; /* only ASCII and DOUBLE matter for our 3 tags */
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

#define DTM_MAX_LEVELS 16
#define DTM_TILE_CACHE_SLOTS 2048

typedef struct {
    int valid;
    int level;
    uint32_t tile_col, tile_row;
    int16_t *data;    /* level_info[level].tile_w * tile_h samples */
    size_t capacity;  /* samples currently allocated in `data` */
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
    double lonMin, lonMax, latMin, latMax; /* level-independent: every overview covers the same geographic extent */
    double nodata;
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
    int scaleCount = 0, tiepointCount = 0, hasNodata = 0;
    char nodataStr[64];
    if (ReadGeoTagsRaw(path, scale, &scaleCount, tiepoint, &tiepointCount, nodataStr, sizeof(nodataStr), &hasNodata) != 0) {
        fprintf(stderr, "dtm: couldn't open %s to read its GeoTIFF tags\n", path);
        return -1;
    }
    if (scaleCount < 2 || tiepointCount < 6) {
        fprintf(stderr, "dtm: %s: missing ModelPixelScale/ModelTiepoint GeoTIFF tags -- not a recognized DTM layout\n", path);
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

    uint16_t bits = 0, sampleFormat = 0, samplesPerPixel = 0;
    TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &bits);
    TIFFGetFieldDefaulted(tiff, TIFFTAG_SAMPLEFORMAT, &sampleFormat);
    TIFFGetFieldDefaulted(tiff, TIFFTAG_SAMPLESPERPIXEL, &samplesPerPixel);
    if (bits != 16 || sampleFormat != SAMPLEFORMAT_UINT || samplesPerPixel != 1) {
        fprintf(stderr, "dtm: %s: expected single-band 16-bit unsigned int (got %u bits, sampleFormat=%u, samples/pixel=%u) -- not a recognized DTM layout\n",
                path, bits, sampleFormat, samplesPerPixel);
        TIFFClose(tiff);
        return -1;
    }

    uint32_t width = 0, height = 0;
    TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &height);

    /* Tiepoint is (I,J,K, X,Y,Z) for one raster corner -- these files always
     * tie (0,0) to the northwest corner, but compute from whatever I,J is
     * given rather than assuming, in case that ever changes. */
    double i0 = tiepoint[0], j0 = tiepoint[1], x0 = tiepoint[3], y0 = tiepoint[4];
    double sx = scale[0], sy = scale[1];
    double lonNW = x0 - i0 * sx;
    double latNW = y0 + j0 * sy;
    double lonSE = lonNW + width * sx;
    double latSE = latNW - height * sy;

    DtmFile *f = calloc(1, sizeof(DtmFile));
    if (!f) { fprintf(stderr, "dtm: out of memory\n"); TIFFClose(tiff); return -1; }
    f->tiff = tiff;
    f->path = strdup(path);
    f->current_dir = 0;
    f->lonMin = lonNW < lonSE ? lonNW : lonSE;
    f->lonMax = lonNW < lonSE ? lonSE : lonNW;
    f->latMin = latSE < latNW ? latSE : latNW;
    f->latMax = latSE < latNW ? latNW : latSE;

    f->nodata = hasNodata ? atof(nodataStr) : 65535.0; /* 65535 is the observed default across every DTM tile inspected so far */

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
    printf("dtm: opened %s: %dx%d px, %d pyramid level%s, lon [%.4f, %.4f] lat [%.4f, %.4f], nodata=%.0f\n",
           path, (int)f->level_info[0].width, (int)f->level_info[0].height, f->num_levels,
           f->num_levels == 1 ? "" : "s", f->lonMin, f->lonMax, f->latMin, f->latMax, f->nodata);
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
static int16_t *GetTile(DtmFile *f, int level, uint32_t tile_col, uint32_t tile_row) {
    LevelInfo *li = &f->level_info[level];
    size_t tile_samples = (size_t)li->tile_w * li->tile_h;

    uint64_t h = (uint64_t)level * 1000003ull + (uint64_t)tile_col * 2654435761ull + (uint64_t)tile_row * 40503ull;
    TileCacheSlot *slot = &f->cache[h % DTM_TILE_CACHE_SLOTS];

    if (slot->valid && slot->level == level && slot->tile_col == tile_col && slot->tile_row == tile_row) {
        return slot->data;
    }

    if (slot->capacity < tile_samples) {
        int16_t *bigger = realloc(slot->data, tile_samples * sizeof(int16_t));
        if (!bigger) { fprintf(stderr, "dtm: out of memory growing a tile cache slot\n"); return NULL; }
        slot->data = bigger;
        slot->capacity = tile_samples;
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
    tmsize_t bytesExpected = (tmsize_t)tile_samples * sizeof(int16_t);
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
    int16_t *tile = GetTile(f, level, tile_col, tile_row);
    if (!tile) return f->nodata;
    uint32_t localCol = col % li->tile_w;
    uint32_t localRow = row % li->tile_h;
    /* The raster is SAMPLEFORMAT_UINT: reinterpret the same 16 bits unsigned
     * (a plain int16_t read would turn values >= 32768, like the 65535
     * nodata sentinel, negative). */
    uint16_t raw = (uint16_t)tile[localRow * li->tile_w + localCol];
    return (double)raw;
}

static int FileCovers(const DtmFile *f, double lat, double lon) {
    return lat >= f->latMin && lat <= f->latMax && lon >= f->lonMin && lon <= f->lonMax;
}

static double SampleFile(DtmFile *f, double lat, double lon, int level) {
    if (level < 0) level = 0;
    if (level >= f->num_levels) level = f->num_levels - 1;
    if (!EnsureLevelInfo(f, level)) return DTM_NODATA;
    LevelInfo *li = &f->level_info[level];

    double colf = (lon - f->lonMin) / (f->lonMax - f->lonMin) * li->width;
    double rowf = (f->latMax - lat) / (f->latMax - f->latMin) * li->height;

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
            double lat = f->latMin + (f->latMax - f->latMin) * yi / 40.0;
            for (int xi = 0; xi <= 40; xi++) {
                double lon = f->lonMin + (f->lonMax - f->lonMin) * xi / 40.0;
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
