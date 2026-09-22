/* Standalone validation test for longley_rice.c against the NTIA/itm
 * reference implementation's own published point-to-point test vectors:
 * https://github.com/NTIA/itm -- p2p.csv (scenario parameters + expected
 * A__db) and pfls.csv (the terrain profile for each scenario, one line
 * per p2p.csv row, in the exact `pfl[]` layout ITM_P2P_TLS expects: np,
 * step__meter, then np+1 elevations). Both files are copied verbatim from
 * that public-domain repo into this repo's root as test fixtures.
 *
 * Usage: ./longley_rice_test [p2p.csv] [pfls.csv]  (both default to the
 * same-named files in the current directory.)
 *
 * No raylib dependency -- matches the dtm_test/dtm_fetch_test convention
 * of a small standalone main() with pass/fail printed per case. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "longley_rice.h"

#define TOLERANCE_DB 0.1

/* Reads one line from f into a freshly malloc'd, NUL-terminated buffer
 * (trailing \n and \r stripped). Returns NULL at EOF. Grows to fit lines
 * of any length -- pfls.csv rows can be tens of KB (one very long terrain
 * profile encoded as a single CSV line). */
static char *ReadLine(FILE *f) {
    size_t cap = 4096, len = 0;
    char *buf = malloc(cap);
    if (!buf) return NULL;
    int c;
    int sawAny = 0;
    while ((c = fgetc(f)) != EOF) {
        sawAny = 1;
        if (c == '\n') break;
        if (len + 1 >= cap) {
            cap *= 2;
            char *grown = realloc(buf, cap);
            if (!grown) { free(buf); return NULL; }
            buf = grown;
        }
        buf[len++] = (char)c;
    }
    if (!sawAny) { free(buf); return NULL; }
    if (len > 0 && buf[len - 1] == '\r') len--;
    buf[len] = '\0';
    return buf;
}

/* Splits a comma-separated line of doubles into a freshly malloc'd array.
 * Returns the element count via *outCount; caller frees the array. */
static double *ParseDoubleCsvLine(char *line, int *outCount) {
    int cap = 8;
    for (char *p = line; *p; p++)
        if (*p == ',') cap++;
    cap += 2;
    double *arr = malloc(sizeof(double) * (size_t)cap);
    int n = 0;
    char *tok = strtok(line, ",");
    while (tok) {
        if (n >= cap) { cap *= 2; arr = realloc(arr, sizeof(double) * (size_t)cap); }
        arr[n++] = atof(tok);
        tok = strtok(NULL, ",");
    }
    *outCount = n;
    return arr;
}

/* Splits a comma-separated header line into column name tokens (in place,
 * NUL-splitting `line`); returns the count and fills `names[]` with
 * pointers into `line`. */
static int ParseHeader(char *line, char *names[], int maxNames) {
    int n = 0;
    char *tok = strtok(line, ",");
    while (tok && n < maxNames) {
        names[n++] = tok;
        tok = strtok(NULL, ",");
    }
    return n;
}

static int ColumnIndex(char *names[], int nNames, const char *want) {
    for (int i = 0; i < nNames; i++)
        if (strcmp(names[i], want) == 0)
            return i;
    fprintf(stderr, "p2p.csv: expected column \"%s\" not found in header\n", want);
    exit(1);
}

int main(int argc, char **argv) {
    const char *p2pPath = argc > 1 ? argv[1] : "p2p.csv";
    const char *pflsPath = argc > 2 ? argv[2] : "pfls.csv";

    FILE *p2pFile = fopen(p2pPath, "r");
    if (!p2pFile) { fprintf(stderr, "could not open %s\n", p2pPath); return 1; }
    FILE *pflsFile = fopen(pflsPath, "r");
    if (!pflsFile) { fprintf(stderr, "could not open %s\n", pflsPath); fclose(p2pFile); return 1; }

    char *headerLine = ReadLine(p2pFile);
    if (!headerLine) { fprintf(stderr, "%s is empty\n", p2pPath); return 1; }

    char *names[32];
    int nNames = ParseHeader(headerLine, names, 32);

    int col_h_tx = ColumnIndex(names, nNames, "h_tx__meter");
    int col_h_rx = ColumnIndex(names, nNames, "h_rx__meter");
    int col_epsilon = ColumnIndex(names, nNames, "epsilon");
    int col_sigma = ColumnIndex(names, nNames, "sigma");
    int col_N_0 = ColumnIndex(names, nNames, "N_0");
    int col_f = ColumnIndex(names, nNames, "f__mhz");
    int col_pol = ColumnIndex(names, nNames, "pol");
    int col_climate = ColumnIndex(names, nNames, "climate");
    int col_time = ColumnIndex(names, nNames, "time");
    int col_location = ColumnIndex(names, nNames, "location");
    int col_situation = ColumnIndex(names, nNames, "situation");
    int col_mdvar = ColumnIndex(names, nNames, "mdvar");
    int col_A_db = ColumnIndex(names, nNames, "A__db");
    free(headerLine);

    int total = 0, passed = 0;
    double maxDiscrepancy = 0.0, sumDiscrepancy = 0.0;
    int worstCase = -1;
    double worstExpected = 0.0, worstComputed = 0.0;

    char *row;
    while ((row = ReadLine(p2pFile)) != NULL) {
        if (row[0] == '\0') { free(row); continue; } /* skip trailing blank line */

        char *pflLine = ReadLine(pflsFile);
        if (!pflLine) {
            fprintf(stderr, "%s ran out of rows before %s did (row %d)\n", pflsPath, p2pPath, total + 1);
            free(row);
            break;
        }

        int nFields = 0;
        double *fields = ParseDoubleCsvLine(row, &nFields);

        int nPfl = 0;
        double *pfl = ParseDoubleCsvLine(pflLine, &nPfl);

        double h_tx__meter = fields[col_h_tx];
        double h_rx__meter = fields[col_h_rx];
        double epsilon = fields[col_epsilon];
        double sigma = fields[col_sigma];
        double N_0 = fields[col_N_0];
        double f__mhz = fields[col_f];
        int pol = (int)fields[col_pol];
        int climate = (int)fields[col_climate];
        double time = fields[col_time];
        double location = fields[col_location];
        double situation = fields[col_situation];
        int mdvar = (int)fields[col_mdvar];
        double expected_A__db = fields[col_A_db];

        double computed_A__db = 0.0;
        long warnings = 0;
        int rtn = ITM_P2P_TLS(h_tx__meter, h_rx__meter, pfl, climate, N_0, f__mhz, pol, epsilon, sigma,
            mdvar, time, location, situation, &computed_A__db, &warnings);

        total++;

        if (rtn != SUCCESS && rtn != SUCCESS_WITH_WARNINGS) {
            printf("case %2d: FAIL  ITM_P2P_TLS returned error %d (expected A__db=%.2f)\n", total, rtn, expected_A__db);
        } else {
            double discrepancy = fabs(computed_A__db - expected_A__db);
            sumDiscrepancy += discrepancy;
            if (discrepancy > maxDiscrepancy) {
                maxDiscrepancy = discrepancy;
                worstCase = total;
                worstExpected = expected_A__db;
                worstComputed = computed_A__db;
            }

            int ok = discrepancy <= TOLERANCE_DB;
            if (ok) passed++;

            printf("case %2d: %-4s expected=%.2f dB  computed=%.2f dB  discrepancy=%.4f dB\n",
                total, ok ? "PASS" : "FAIL", expected_A__db, computed_A__db, discrepancy);
            if (rtn == SUCCESS_WITH_WARNINGS)
                printf("            (SUCCESS_WITH_WARNINGS, warnings=0x%04lx)\n", (unsigned long)warnings);
        }

        free(fields);
        free(pfl);
        free(row);
        free(pflLine);
    }

    fclose(p2pFile);
    fclose(pflsFile);

    double meanDiscrepancy = total > 0 ? sumDiscrepancy / total : 0.0;

    printf("\n==== longley_rice_test summary ====\n");
    printf("%d/%d cases within %.2f dB tolerance\n", passed, total, TOLERANCE_DB);
    printf("max discrepancy:  %.4f dB", maxDiscrepancy);
    if (worstCase >= 0)
        printf("  (case %d: expected=%.2f dB, computed=%.2f dB)", worstCase, worstExpected, worstComputed);
    printf("\n");
    printf("mean discrepancy: %.4f dB\n", meanDiscrepancy);

    return (passed == total) ? 0 : 1;
}
